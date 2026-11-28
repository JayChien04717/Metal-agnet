# -*- coding: utf-8 -*-
r"""HFSS Eigenmode verification of Purcell filter + readout resonator frequencies for Design11_Transmon2Q_Manhattan.

Renders the loaded readout fixture (feedline with C_in, TapJunctions, C_p, PurcellReadout modules with PEC straps,
open extensions, and Q1 / Q2 with Manhattan pads & lumped Lj; flux lines as ground, TC omitted to keep mesh focused
on readout channels).

Solver defaults (skill/em_simulation.md):
  - HFSS eigenmode, max_passes = 6, basis_order = 1
  - Explicit 10 um length mesh on all metal objects (mesh_all)
  - Auto wire bonds every 350 um on etched paths + PEC straps at PurcellReadout.strap_points()
  - 50 Ohm ports at feed_cin start (P8 side) and feed_out start (P7 side)

Fast tuning options (all opt-in; defaults preserve the full-chip settings above):
  - --channel Q1|Q2     render only that channel (filter, resonator, its qubit/claw, and local feed section)
  - --local-mesh        apply explicit 10 um mesh only to filter/resonator CPW traces, their gaps (non-model
                        sheets), coupling hammer, C_p, and qubit readout claw; rest of chip adapts freely
  - --no-wire-bonds     disable auto wire bonds on etched paths while tuning (keeps hammer PEC straps)
  - --pass-by-pass      solve one adaptive pass at a time, saving eigenfrequencies/Q/elements after every pass
  - --early-stop-mhz X  with --pass-by-pass, stop early once every channel mode moves < X MHz between passes
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import traceback
from contextlib import chdir
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(ROOT / 'code' / 'resonator'),
    str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
]
import project_paths
project_paths.configure(__file__)

import numpy as np
from qiskit_metal import Dict
import build as D11
import hfss_dirs
import purcell_resonator_pair_sim as S
import save_aedt
from run_record import create_run, finish_run

CHANNELS = ('Q2', 'Q1')  # sorted by frequency: Q2 (6.15 GHz), Q1 (6.35 GHz)
Q_MIN = 10.0
BAND_DEFAULT = (5.6, 6.9)


def build_fixture(channel: str | None = None, wire_bonds: bool = True):
    d = D11.build_design(variant='sim', channels=True, flux_as_ground=True, check_keepout_um=None)
    t = d.qgeometry.tables['path']
    if wire_bonds:
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)
        # No wire bonds on the filter / resonator CPWs (PurcellReadout and the resonator extension res_open): the
        # fabricated chip has no air bridges there (user 2026-11-26). Hammer straps (strap_points) are drawn separately.
        no_wb = [d.components[f'{p}_{q}'].id for q in ('Q1', 'Q2') for p in ('prr', 'res_open')
                 if f'{p}_{q}' in d.components]
        t.loc[t['component'].isin(no_wb), 'hfss_wire_bonds'] = False
    else:
        t['hfss_wire_bonds'] = False

    if channel == 'Q1':
        comps = ['feed_cin', 'cin', 'feed_in', 'n0_Q1', 'feed_mid',
                 'Q1', 'stub_p_Q1', 'cp_Q1', 'prr_Q1', 'res_open_Q1']
        ports = [('feed_cin', 'start', 50), ('feed_mid', 'end', 50)]
    elif channel == 'Q2':
        comps = ['feed_mid', 'n0_Q2', 'feed_out',
                 'Q2', 'stub_p_Q2', 'cp_Q2', 'prr_Q2', 'res_open_Q2']
        ports = [('feed_mid', 'start', 50), ('feed_out', 'start', 50)]
    else:
        comps = ['feed_cin', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2']
        for q in ('Q1', 'Q2'):
            comps += [f'{p}_{q}' for p in ('n0', 'stub_p', 'cp', 'prr', 'res_open')]
        ports = [('feed_cin', 'start', 50), ('feed_out', 'start', 50)]
    opens = []
    return d, comps, ports, opens


def _geom_polygons(geom):
    """Yield exterior coordinate lists for a Shapely Polygon or MultiPolygon."""
    if geom is None or geom.is_empty:
        return
    gtype = geom.geom_type
    if gtype == 'Polygon':
        yield list(geom.exterior.coords)[:-1]
    elif gtype in ('MultiPolygon', 'GeometryCollection'):
        for sub in geom.geoms:
            yield from _geom_polygons(sub)


def apply_local_mesh(renderer, design, target_channels: tuple[str, ...], mesh_len: str = S.MESH) -> dict:
    """Apply explicit length mesh ONLY to filter/resonator CPW traces, their gaps (as non-model sheets),
    the coupling hammer, C_p, and the qubit readout claw. Rest of the chip adapts freely."""
    id2n = {c.id: n for n, c in design.components.items()}
    metal_objs = []
    gap_objs = []

    metal_comps = set()
    gap_comps = set()
    for q in target_channels:
        metal_comps.update([f'prr_{q}', f'res_open_{q}', f'cp_{q}'])
        gap_comps.update([f'prr_{q}', f'res_open_{q}', f'cp_{q}'])

    # 1. Existing metal sheet objects in HFSS
    for kind in ('path', 'poly'):
        for r in design.qgeometry.tables[kind].itertuples():
            cn = id2n[int(r.component)]
            if bool(r.subtract):
                continue
            if cn in metal_comps:
                metal_objs.append(f'{r.name}_{cn}')
            elif cn in target_channels and r.name == 'readout_connector_arm':
                metal_objs.append(f'{r.name}_{cn}')

    # 2. Non-model sheet objects for the CPW/hammer/claw gaps (subtracted from ground plane during render)
    poly_idx = 0
    for kind in ('path', 'poly'):
        for r in design.qgeometry.tables[kind].itertuples():
            cn = id2n[int(r.component)]
            if not bool(r.subtract):
                continue
            want = (cn in gap_comps) or (cn in target_channels and r.name == 'readout_connector_etcher')
            if not want:
                continue
            geom = r.geometry
            if kind == 'path':
                geom = geom.buffer(float(r.width) / 2.0, cap_style=2, join_style=2)
            geom = geom.simplify(0.0005, preserve_topology=True)
            for coords in _geom_polygons(geom):
                if len(coords) < 3:
                    continue
                pts = [[f'{x:.6f}mm', f'{y:.6f}mm', '0mm'] for x, y in coords]
                poly_idx += 1
                nm = f'lm_{cn}_{r.name}_{poly_idx}'
                obj = renderer.modeler.draw_polyline(pts, closed=True, nonmodel=True, name=nm)
                gap_objs.append(str(obj))

    all_objs = metal_objs + gap_objs
    renderer.modeler.mesh_length('local_mesh', all_objs, MaxLength=mesh_len)
    return dict(mesh_length=mesh_len, metal_objects=metal_objs, gap_objects=gap_objs, total_objects=len(all_objs))


def _modes_list(freqs, imags) -> list[dict]:
    modes = []
    for i, (f, im) in enumerate(zip(freqs, imags)):
        f, im = float(f), float(im)
        kappa = 2.0 * abs(im) * 1e3
        qval = f / (2.0 * abs(im)) if im else float('inf')
        modes.append(dict(idx=i, f_GHz=f, im_GHz=im, kappa_MHz=kappa, Q=qval))
    return modes


def analyse_modes(raw_modes: list[dict], target_channels: tuple[str, ...]) -> dict:
    """Filter high-Q in-band modes and invert each channel pair with purcell_resonator_pair_sim.invert."""
    if len(target_channels) == 1:
        ch = target_channels[0]
        f0 = D11.PURCELL_CHANNELS[ch]['f_GHz']
        band = (f0 - 0.35, f0 + 0.35)
    else:
        band = BAND_DEFAULT

    keep = sorted(
        [m for m in raw_modes if band[0] <= m['f_GHz'] <= band[1] and m['Q'] >= Q_MIN],
        key=lambda m: m['f_GHz'],
    )
    dropped = [m for m in raw_modes if m not in keep]
    channels = {}

    if len(target_channels) == 1:
        ch = target_channels[0]
        f0 = D11.PURCELL_CHANNELS[ch]['f_GHz']
        if len(keep) >= 2:
            # Pick the two high-Q modes nearest target frequency, then sort by frequency
            pair = sorted(sorted(keep, key=lambda m: abs(m['f_GHz'] - f0))[:2], key=lambda m: m['f_GHz'])
            channels[ch] = _invert_pair(ch, pair)
    else:
        if len(keep) >= 4:
            for ch, pair in (('Q2', keep[0:2]), ('Q1', keep[2:4])):
                channels[ch] = _invert_pair(ch, pair)

    return dict(
        all_modes=raw_modes,
        used_modes=keep,
        dropped_modes=dropped,
        in_band_high_q_count=len(keep),
        channels=channels,
    )


def _invert_pair(ch: str, pair: list[dict]) -> dict:
    z1 = complex(pair[0]['f_GHz'], abs(pair[0]['im_GHz']))
    z2 = complex(pair[1]['f_GHz'], abs(pair[1]['im_GHz']))
    inv = S.invert(z1, z2)
    kp = inv['kappa_p_MHz']
    J = abs(inv['J_MHz'])
    inv['kappa_r_est_MHz'] = float((4.0 * J * J / kp) if kp > 0 else float('nan'))
    inv['filter_like_mode_GHz'] = max(pair, key=lambda m: abs(m['im_GHz']))['f_GHz']
    inv['resonator_like_mode_GHz'] = min(pair, key=lambda m: abs(m['im_GHz']))['f_GHz']
    inv['modes'] = pair
    inv['target_f_GHz'] = D11.PURCELL_CHANNELS[ch]['f_GHz']
    inv['df_p_MHz'] = (inv['f_p_GHz'] - inv['target_f_GHz']) * 1e3
    inv['df_r_MHz'] = (inv['f_r_GHz'] - inv['target_f_GHz']) * 1e3
    inv['detuning_fr_minus_fp_MHz'] = (inv['f_r_GHz'] - inv['f_p_GHz']) * 1e3
    return inv


def _read_sd_pass(results_dir: Path, pass_no: int) -> tuple[list[dict], int | None, float | None]:
    """Read eigenmodes (f, f'' = Im(f)) for a specific adaptive pass directly from AEDT *.sd and *.EMConvData."""
    modes: list[dict] = []
    # Find the *.sd file whose AdaptivePass == pass_no (prefer Adaptive_<pass_no> over LastAdaptive)
    candidates = sorted(results_dir.glob('*_SOL*_M*_V*.sd'))
    chosen_text = None
    for f in candidates:
        t = f.read_text(encoding='utf-8', errors='replace')
        m_p = re.search(r'AdaptivePass=(\d+)', t)
        if not m_p or int(m_p.group(1)) != pass_no:
            continue
        chosen_text = t
        m_n = re.search(r"Name='([^']+)'", t)
        if m_n and m_n.group(1) == f'Adaptive_{pass_no}':
            break
    if chosen_text is not None:
        m_modes = re.search(r'Modes\(([^)]*)\)', chosen_text)
        if m_modes and m_modes.group(1).strip():
            vals = [float(x) for x in m_modes.group(1).split(',')]
            for idx, i in enumerate(range(0, len(vals), 2)):
                f_ghz = vals[i] / 1e9
                im_ghz = vals[i + 1] / 1e9
                kappa = 2.0 * abs(im_ghz) * 1e3
                qval = f_ghz / (2.0 * abs(im_ghz)) if im_ghz else float('inf')
                modes.append(dict(idx=idx, f_GHz=f_ghz, im_GHz=im_ghz, kappa_MHz=kappa, Q=qval))

    elems = None
    delta_f_pct = None
    for f in results_dir.glob('*.EMConvData'):
        t = f.read_text(encoding='utf-8', errors='replace')
        em = dict((int(p), int(v)) for p, v in
                  re.findall(r'PassNo=(\d+)\s*Intrinsics=\S*\s*SolvedElements\((\d+)', t))
        dm = dict((int(p), float(v)) for p, v in
                  re.findall(r'PassNo=(\d+).*?MaxDeltaFreq\(([\d.eE+-]+)', t, flags=re.S))
        if pass_no in em:
            elems = em[pass_no]
        if pass_no in dm:
            delta_f_pct = dm[pass_no]
    return modes, elems, delta_f_pct


def _edit_eigen_max_passes(renderer, setup_name: str, setup_cfg: dict, max_passes: int):
    renderer.pinfo.design._setup_module.EditSetup(
        setup_name,
        [
            'NAME:' + setup_name,
            'MinimumFrequency:=', f"{setup_cfg['min_freq_ghz']}GHz",
            'NumModes:=', int(setup_cfg['n_modes']),
            'MaxDeltaFreq:=', float(setup_cfg['max_delta_f']),
            'ConvergeOnRealFreq:=', True,
            'MaximumPasses:=', int(max_passes),
            'MinimumPasses:=', int(setup_cfg['min_passes']),
            'MinimumConvergedPasses:=', int(setup_cfg['min_converged']),
            'PercentRefinement:=', int(setup_cfg['pct_refinement']),
            'IsEnabled:=', True,
            'BasisOrder:=', int(setup_cfg['basis_order']),
        ],
    )


def run_eigen_custom(
    folder: Path,
    d,
    comps: list[str],
    ports: list,
    opens: list,
    *,
    minf: float,
    n_modes: int,
    passes: int,
    target_channels: tuple[str, ...],
    local_mesh: bool,
    wire_bonds: bool,
    pass_by_pass: bool,
    early_stop_mhz: float | None,
    tag: str = 'E',
) -> tuple[dict, dict, list[dict]]:
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim

    info: dict = {}
    pass_history: list[dict] = []
    max_delta_f = 0.01 if (pass_by_pass and early_stop_mhz is not None) else 0.1
    min_passes = 1 if pass_by_pass else min(2, passes)
    setup_cfg = dict(
        min_freq_ghz=minf,
        n_modes=n_modes,
        max_delta_f=max_delta_f,
        max_passes=1 if pass_by_pass else passes,
        min_passes=min_passes,
        min_converged=1,
        pct_refinement=20,
        basis_order=1,
        vars=Dict(),
    )
    t_start = time.time()
    hist_path = folder / 'raw' / 'pass_history.json'
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    sum_hist_path = folder / 'summaries' / 'pass_history.json'

    class Eig(EigenmodeSim):
        def _analyze(self):
            info['straps'] = S.draw_straps(self.renderer, d, comps)
            if local_mesh:
                info['mesh'] = apply_local_mesh(self.renderer, d, target_channels, mesh_len=S.MESH)
            else:
                objs = S.mesh_all(self.renderer, d, comps)
                info['mesh'] = dict(mesh_length=S.MESH, mode='mesh_all', total_objects=len(objs))

            self.sim_setup_name = self.renderer.initialize_eigenmode(**self.setup)
            setup = self.renderer.pinfo.setup

            if not pass_by_pass:
                self.renderer.analyze_setup(self.sim_setup_name)
                self.compute_convergences()
                return

            prev_freqs: list[float] | None = None
            prev_channels: dict | None = None
            design_name = self.renderer.pinfo.design.name
            results_dir = ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'aedt' / 'D11.aedtresults' / f'{design_name}.results'
            for p in range(1, passes + 1):
                if p > 1:
                    _edit_eigen_max_passes(self.renderer, self.sim_setup_name, setup_cfg, p)
                t_p0 = time.time()
                self.renderer.analyze_setup(self.sim_setup_name)
                pass_min = (time.time() - t_p0) / 60.0
                elapsed_min = (time.time() - t_start) / 60.0

                modes, elems, delta_f_pct = _read_sd_pass(results_dir, p)
                if not modes:
                    # Fallback to pyEPR COM readout (note: pyEPR returns kappa = 2*Im(f), so divide by 2)
                    f_arr, kp_arr = setup.get_solutions().eigenmodes()
                    im_arr = [0.5 * v for v in kp_arr] if kp_arr is not None else [0.0] * len(f_arr)
                    modes = _modes_list(f_arr, im_arr)
                ana_p = analyse_modes(modes, target_channels)

                # Extract sorted channel mode frequencies for pass-to-pass movement check
                cur_freqs = []
                for ch in sorted(ana_p['channels']):
                    for m in ana_p['channels'][ch]['modes']:
                        cur_freqs.append(m['f_GHz'])
                cur_freqs.sort()

                mode_deltas = None
                max_mode_delta = None
                bare_deltas = {}
                if prev_freqs is not None and len(cur_freqs) == len(prev_freqs) and len(cur_freqs) == 2 * len(target_channels):
                    mode_deltas = [round((c - pr) * 1e3, 3) for c, pr in zip(cur_freqs, prev_freqs)]
                    max_mode_delta = float(max(abs(x) for x in mode_deltas))
                if prev_channels is not None:
                    for ch, inv in ana_p['channels'].items():
                        if ch in prev_channels:
                            bare_deltas[ch] = dict(
                                df_p_MHz=round((inv['f_p_GHz'] - prev_channels[ch]['f_p_GHz']) * 1e3, 3),
                                df_r_MHz=round((inv['f_r_GHz'] - prev_channels[ch]['f_r_GHz']) * 1e3, 3),
                            )

                early_stopped = bool(
                    early_stop_mhz is not None
                    and p >= 2
                    and max_mode_delta is not None
                    and max_mode_delta < early_stop_mhz
                )

                rec = dict(
                    pass_no=p,
                    pass_min=round(pass_min, 3),
                    elapsed_min=round(elapsed_min, 3),
                    solved_elements=elems,
                    max_delta_f_pct=delta_f_pct,
                    analysis=ana_p,
                    mode_deltas_MHz=mode_deltas,
                    max_mode_delta_MHz=max_mode_delta,
                    bare_deltas_MHz=bare_deltas,
                    early_stopped=early_stopped,
                )
                pass_history.append(rec)
                S.write(hist_path, pass_history)
                S.write(sum_hist_path, pass_history)

                ch_str = ' | '.join(
                    f"{ch}: fp={v['f_p_GHz']:.5f} fr={v['f_r_GHz']:.5f} J={v['J_MHz']:.2f} kp={v['kappa_p_MHz']:.2f}"
                    for ch, v in ana_p['channels'].items()
                )
                print(
                    f"[PASS {p}/{passes}] {pass_min:.2f} min (tot {elapsed_min:.2f} min) "
                    f"elems={elems} dF%={delta_f_pct} max_dMHz={max_mode_delta} | {ch_str}",
                    flush=True,
                )

                prev_freqs = cur_freqs if len(cur_freqs) == 2 * len(target_channels) else None
                prev_channels = ana_p['channels'] if ana_p['channels'] else None
                if early_stopped:
                    print(f"[EARLY STOP] max mode movement {max_mode_delta:.3f} MHz < {early_stop_mhz} MHz", flush=True)
                    break

    eig = Eig(d, 'hfss')
    if wire_bonds:
        eig.renderer._options.update(S.WB)
    for k, v in setup_cfg.items():
        eig.setup[k] = v

    with chdir(folder / 'raw'):
        eig.run(name=f'{tag}_{folder.name[-8:]}', components=comps, open_terminations=opens, port_list=ports)

    design_name = eig.renderer.pinfo.design.name
    results_dir = ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'aedt' / 'D11.aedtresults' / f'{design_name}.results'
    last_p = len(pass_history) if pass_by_pass else passes
    modes, elems, delta_f_pct = _read_sd_pass(results_dir, last_p)
    if not modes:
        f_arr, kp_arr = eig.renderer.pinfo.setup.get_solutions().eigenmodes()
        im_arr = [0.5 * v for v in kp_arr] if kp_arr is not None else [0.0] * len(f_arr)
        modes = _modes_list(f_arr, im_arr)
    ana = analyse_modes(modes, target_channels)
    raw = dict(
        frequencies_GHz=[m['f_GHz'] for m in modes],
        imag_GHz=[m['im_GHz'] for m in modes],
        passes=int(len(pass_history)) if pass_by_pass else passes,
        max_passes=passes,
        solved_elements=elems,
        last_delta_f_pct=delta_f_pct,
        minutes=(time.time() - t_start) / 60.0,
        aedt_design=design_name,
        aedt_results=str(results_dir.relative_to(ROOT)),
        n_straps=len(info.get('straps', [])),
        mesh=info.get('mesh'),
        wire_bonds=wire_bonds,
        min_freq_GHz=minf,
        n_modes=n_modes,
        pass_by_pass=pass_by_pass,
        early_stop_mhz=early_stop_mhz,
        early_stopped=bool(pass_history and pass_history[-1].get('early_stopped', False)),
    )
    S.write(folder / 'raw' / 'eigen_modes.json', raw)
    return raw, ana, pass_history


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=6)
    ap.add_argument('--min-freq', type=float, default=5.4)
    ap.add_argument('--n-modes', type=int, default=6)
    ap.add_argument('--tag', default='ro_eig')
    ap.add_argument('--channel', choices=['Q1', 'Q2'], default=None,
                    help='Render only one channel (filter, resonator, its qubit/claw, and feed section)')
    ap.add_argument('--local-mesh', action='store_true',
                    help='Apply 10 um mesh only to filter/resonator traces, gaps, hammer, C_p, and claw')
    ap.add_argument('--no-wire-bonds', action='store_true',
                    help='Disable auto wire bonds while tuning (keeps hammer PEC straps)')
    ap.add_argument('--pass-by-pass', action='store_true',
                    help='Solve pass-by-pass, saving eigenfrequencies and Q after every pass')
    ap.add_argument('--early-stop-mhz', type=float, default=None,
                    help='With --pass-by-pass, stop once every channel mode moves less than this many MHz')
    a = ap.parse_args()

    target_channels = (a.channel,) if a.channel else CHANNELS
    wire_bonds = not a.no_wire_bonds
    d, comps, ports, opens = build_fixture(channel=a.channel, wire_bonds=wire_bonds)

    folder = create_run(
        'Design11_Transmon2Q_Manhattan', 'resonator', f'hfss_eigen_{a.tag}',
        entrypoint=__file__,
        inputs=[ROOT / 'Design11_Transmon2Q_Manhattan/design.json',
                ROOT / 'Design11_Transmon2Q_Manhattan/layout/build.py'],
        parameters=dict(
            channel=a.channel or 'both',
            components=comps, ports=ports, open_terminations=opens,
            max_passes=a.passes, min_freq_ghz=a.min_freq, n_modes=a.n_modes,
            basis_order=1, mesh_length=S.MESH, local_mesh=a.local_mesh,
            wire_bonds=S.WB if wire_bonds else False,
            pass_by_pass=a.pass_by_pass, early_stop_mhz=a.early_stop_mhz,
            purcell_channels=D11.PURCELL_CHANNELS,
            hammer_line_spacing_um=D11.HAMMER_LINE_SPACING_UM,
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        hfss_dirs.use_d_drive()
        project = ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'aedt' / 'D11.aedt'
        project.parent.mkdir(parents=True, exist_ok=True)
        save_aedt.use_project(project)

        use_custom = a.channel is not None or a.local_mesh or a.no_wire_bonds or a.pass_by_pass or (a.early_stop_mhz is not None)
        if not use_custom:
            raw = S.run_eigen(folder, d, comps, ports, opens, minf=a.min_freq, n_modes=a.n_modes, passes=a.passes, tag='E')
            ana = analyse_modes(_modes_list(raw['frequencies_GHz'], raw['imag_GHz']), target_channels)
            pass_history = []
        else:
            raw, ana, pass_history = run_eigen_custom(
                folder, d, comps, ports, opens,
                minf=a.min_freq, n_modes=a.n_modes, passes=a.passes,
                target_channels=target_channels,
                local_mesh=a.local_mesh, wire_bonds=wire_bonds,
                pass_by_pass=a.pass_by_pass, early_stop_mhz=a.early_stop_mhz,
                tag=f'E{a.channel}' if a.channel else 'E',
            )

        summary = dict(raw=raw, analysis=ana, pass_history=pass_history)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, indent=2, default=float), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
