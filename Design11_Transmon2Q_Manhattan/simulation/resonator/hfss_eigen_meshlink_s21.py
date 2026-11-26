# -*- coding: utf-8 -*-
r"""Fast 3-pass HFSS Eigenmode + Mesh-Link Driven S21 for Design11_Transmon2Q_Manhattan.

1. Runs HFSS Eigenmode (max_passes=3, basis_order=1, explicit 10 um length mesh, 350 um auto wire bonds,
   PEC straps at PurcellReadout.strap_points(), 50 Ohm ports on feedline) to verify Purcell & readout
   resonator frequencies (f_p, f_r, J, kappa_p, kappa_r).
2. Creates a DrivenModal design with identical geometry & ports in the same AEDT project, links its mesh
   to the Eigenmode LastAdaptive mesh (MeshLink ImportMesh=True, 1 pass), and sweeps S21 (one wide Fast sweep, 1 MHz step; no fine sweep).

User-requested fast pass budget (3 passes eigenmode + 1 pass Mesh-Link driven S21) recorded in run.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
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
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
import build as D11
import hfss_dirs
import purcell_resonator_pair_sim as S
import save_aedt
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
CHANNELS = ('Q2', 'Q1')


def build_fixture():
    d = D11.build_design(variant='sim', channels=True, flux_as_ground=True, check_keepout_um=None)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    comps = ['feed_cin', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2']
    for q in ('Q1', 'Q2'):
        comps += [f'{p}_{q}' for p in ('n0', 'stub_p', 'cp', 'prr', 'res_open')]
    ports = [('feed_cin', 'start', 50), ('feed_out', 'start', 50)]
    opens = []
    return d, comps, ports, opens


def analyse_eigen(raw: dict) -> dict:
    freqs = np.asarray(raw['frequencies_GHz'], float)
    imags = np.asarray(raw['imag_GHz'], float)
    modes = [dict(idx=i, f_GHz=float(f), im_GHz=float(im), kappa_MHz=float(2.0 * abs(im) * 1e3))
             for i, (f, im) in enumerate(zip(freqs, imags))]
    in_band = sorted([m for m in modes if 5.7 <= m['f_GHz'] <= 6.8], key=lambda m: m['f_GHz'])
    channels = {}
    if len(in_band) >= 4:
        for ch, pair in (('Q2', in_band[0:2]), ('Q1', in_band[2:4])):
            z1 = pair[0]['f_GHz'] + 1j * abs(pair[0]['im_GHz'])
            z2 = pair[1]['f_GHz'] + 1j * abs(pair[1]['im_GHz'])
            inv = S.invert(z1, z2)
            kp = inv['kappa_p_MHz']
            J = abs(inv['J_MHz'])
            inv['kappa_r_est_MHz'] = float(4.0 * J * J / kp) if kp > 0 else float('nan')
            inv['modes'] = pair
            inv['target_f_GHz'] = D11.PURCELL_CHANNELS[ch]['f_GHz']
            inv['df_p_MHz'] = (inv['f_p_GHz'] - inv['target_f_GHz']) * 1e3
            inv['df_r_MHz'] = (inv['f_r_GHz'] - inv['target_f_GHz']) * 1e3
            channels[ch] = inv
    return dict(all_modes=modes, in_band_count=len(in_band), channels=channels)


def run_meshlink_driven(folder: Path, d, comps, ports, opens, eig_design: str, band, tag='S3'):
    """Run 1-pass DrivenModal S21 linked to eig_design:Setup:LastAdaptive mesh."""
    from scipy.signal import find_peaks
    lo, hi = band
    f0 = 0.5 * (lo + hi)
    info = {'mesh_link_applied': False}

    class MeshLinkDrv(ScatteringImpedanceSim):
        def _analyze(self):
            info['straps'] = S.draw_straps(self.renderer, d, comps)
            info['mesh_objects'] = len(S.mesh_all(self.renderer, d, comps))
            try:
                oDes = self.renderer.pinfo.design._design
                oMod = oDes.GetModule('AnalysisSetup')
                oMod.EditSetup('Setup', [
                    'NAME:Setup',
                    'Frequency:=', f'{f0:.4f}GHz',
                    'MaxDeltaS:=', 0.02,
                    'MaximumPasses:=', 1,
                    'MinimumPasses:=', 1,
                    'MinimumConvergedPasses:=', 1,
                    'PercentRefinement:=', 20,
                    'BasisOrder:=', 1,
                    [
                        'NAME:MeshLink',
                        'ImportMesh:=', True,
                        'Project:=', 'This Project*',
                        'Product:=', 'ElectronicsDesktop',
                        'Design:=', eig_design,
                        'Soln:=', 'Setup : LastAdaptive',
                        ['NAME:Params'],
                        'ForceSourceToSolve:=', False,
                        'PreservePartnerSoln:=', False,
                        'PathRelativeTo:=', 'TargetProject',
                        'ApplyMeshOp:=', False,
                    ],
                ])
                info['mesh_link_applied'] = True
                print(f'MeshLink configured -> {eig_design} : Setup : LastAdaptive', flush=True)
            except Exception as exc:
                info['mesh_link_error'] = str(exc)
                print('WARNING: MeshLink failed, falling back to 1-pass seeded mesh:', exc, flush=True)
            super()._analyze()

    sim = MeshLinkDrv(d, 'hfss')
    sim.renderer._options.update(S.WB)
    sim.setup.update(freq_ghz=f0, max_delta_s=0.02, max_passes=1, min_passes=1, min_converged=1,
                     pct_refinement=20, basis_order=1, vars=Dict())
    sim.setup.sweep_setup.update(name='Wide', start_ghz=lo, stop_ghz=hi,
                                 count=int(round((hi - lo) / 0.001)) + 1, type='Fast', save_fields=False)
    t0 = time.time()
    sim.run(name=f'{tag}_{folder.name[-4:]}', components=comps, open_terminations=opens, port_list=ports)
    _, _, Sdf = sim.renderer.get_params(['S11', 'S21'])
    Sdf.to_csv(folder / 'raw' / 'S_wide.csv')
    fw, sw = Sdf.index.to_numpy(float), Sdf['S21'].to_numpy(complex)
    idx, _ = find_peaks(-20 * np.log10(np.maximum(abs(sw), 1e-12)), prominence=1)
    sim.renderer.pinfo.project._project.Save()

    # Plot S21
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ff, sf = S.load_s(folder / 'raw' / 'S_wide.csv')
    db = 20 * np.log10(np.maximum(abs(sf), 1e-12))
    didx, dprops = find_peaks(-db, prominence=2)
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    ax.plot(ff, db, 'k-', lw=1.6, label='HFSS Mesh-Link Driven |S21|')
    for ch in CHANNELS:
        ax.axvline(D11.PURCELL_CHANNELS[ch]['f_GHz'], ls=':', alpha=0.7, label=f"{ch} target ({D11.PURCELL_CHANNELS[ch]['f_GHz']} GHz)")
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('|S21| (dB)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(f"Design11 Option 1 (10/6 CPW, 3um gnd strip) Mesh-Link S21 (linked={info['mesh_link_applied']})")
    fig.tight_layout()
    (folder / 'figures').mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / 'figures' / 's21_meshlink.png', dpi=150)
    plt.close(fig)

    return dict(
        minutes=(time.time() - t0) / 60.0,
        mesh_link_applied=info['mesh_link_applied'],
        mesh_link_error=info.get('mesh_link_error'),
        aedt_design=sim.renderer.pinfo.design.name,
        n_straps=len(info['straps']),
        mesh_objects=info['mesh_objects'],
        wide_dips_GHz=[float(v) for v in fw[idx]],
        wide_dip_details=[dict(f_GHz=float(ff[i]), dB=float(db[i]), prominence_dB=float(p))
                       for i, p in zip(didx, dprops['prominences'])],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=3)
    ap.add_argument('--min-freq', type=float, default=5.4)
    ap.add_argument('--n-modes', type=int, default=6)
    ap.add_argument('--band', default='5.85,6.65')
    ap.add_argument('--tag', default='r8_opt1')
    a = ap.parse_args()
    lo, hi = map(float, a.band.split(','))

    folder = create_run(
        DESIGN, 'resonator', f'hfss_eig3_meshlink_s21_{a.tag}',
        entrypoint=__file__,
        inputs=[ROOT / DESIGN / 'design.json', ROOT / DESIGN / 'layout/build.py'],
        parameters=dict(
            eigen_max_passes=a.passes, driven_max_passes=1, mesh_link=True,
            pass_budget_note='User-requested fast 3-pass eigenmode + 1-pass Mesh-Link driven S21 check',
            min_freq_ghz=a.min_freq, n_modes=a.n_modes, band_GHz=[lo, hi],
            basis_order=1, mesh_length=S.MESH, wire_bonds=S.WB,
            purcell_channels=D11.PURCELL_CHANNELS,
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        hfss_dirs.use_d_drive()
        project = ROOT / DESIGN / 'simulation' / 'aedt' / 'D11.aedt'
        project.parent.mkdir(parents=True, exist_ok=True)
        save_aedt.use_project(project)

        d_eig, comps, ports, opens = build_fixture()
        eig_raw = S.run_eigen(folder, d_eig, comps, ports, opens, minf=a.min_freq, n_modes=a.n_modes,
                              passes=a.passes, tag='E3')
        eig_ana = analyse_eigen(eig_raw)
        print('=== EIGENMODE RESULT ===', flush=True)
        print(json.dumps(eig_ana, indent=2), flush=True)

        d_drv, comps2, ports2, opens2 = build_fixture()
        s21_res = run_meshlink_driven(folder, d_drv, comps2, ports2, opens2, eig_raw['aedt_design'], (lo, hi), tag='S3')
        print('=== MESHLINK S21 RESULT ===', flush=True)
        print(json.dumps(s21_res, indent=2), flush=True)

        summary = dict(eigen_raw=eig_raw, eigen_analysis=eig_ana, meshlink_s21=s21_res)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print('EIG3_MESHLINK_S21_OK', flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
