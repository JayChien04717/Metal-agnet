"""Q3D circuit with distributed filter-readout coupling along the parallel CPW runs.

The earlier circuit treats the Purcell filter and readout as independent TLs joined only
by the lumped hammer. Where the two lines run side by side (centre spacing = hammer
spacing), this script replaces the single TLs with a two-conductor coupled-line section.

Per-unit-length mutual capacitance C'_m comes from Q3D: two parallel straight CPWs, all
four ends open, solved at two lengths; the length difference removes end effects
(extended-CPW reference approach). Self C' keeps the conformal-mapping value of the
existing model, so C'_m = 0 reproduces the old deembedded model exactly (asserted).
Quasi-TEM homogeneous assumption: L' = C'^-1 / v^2 with the single model velocity.
No parameters are fitted; only the pole extraction of each spectrum is a fit.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import pandas as pd
import hammer_deembed as H
from hammer_deembed import P, write, geometry_dump
import hammer_fixed_s21 as F
from hammer_pole_fit import fit_poles
from qiskit_metal import designs, draw, Dict
from qiskit_metal.qlibrary.core import QComponent
from run_record import create_run, finish_run


class ParallelCPW(QComponent):
    """Two straight parallel CPWs along x, centre spacing `spacing`, all ends are pins."""
    default_options = Dict(length='1mm', spacing='200um', trace_width='16um', trace_gap='9um')

    def make(self):
        p = self.p
        for tag, y in (('top', p.spacing/2), ('bot', -p.spacing/2)):
            line = draw.LineString([(-p.length/2, y), (p.length/2, y)])
            self.add_qgeometry('path', {f'line_{tag}': line}, width=p.trace_width, layer=1)
            self.add_qgeometry('path', {f'cut_{tag}': line}, width=p.trace_width+2*p.trace_gap, subtract=True, layer=1)
            for end, s in (('a', -1), ('b', 1)):
                self.add_pin(f'{end}_{tag}', points=np.array([[s*(p.length/2-p.trace_width), y], [s*p.length/2, y]]),
                             width=p.trace_width, gap=p.trace_gap, input_as_norm=True)


def solve_pair(length_mm, spacing_mm, out, name, max_passes=20, percent_error=.1):
    from qiskit_metal.analyses.quantization import LOManalysis
    d = designs.DesignPlanar()
    d.chips.main.material = P.d5.SI_MATERIAL
    d.chips.main.size.update(size_x='5.2mm', size_y='4.3mm', size_z='-500um', center_x='0mm', center_y='0mm')
    ParallelCPW(d, 'pair', options=dict(**P.CPW, length=f'{length_mm}mm', spacing=f'{spacing_mm*1000:g}um'))
    out.mkdir(parents=True, exist_ok=False)
    write(out/'geometry.json', geometry_dump(d))
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer.options['x_buffer_width_mm'] = 0
    lom.sim.renderer.options['y_buffer_width_mm'] = 0
    lom.sim.setup.update(max_passes=max_passes, min_passes=2, percent_error=percent_error, freq_ghz=6.5,
        solution_order='Medium', auto_increase_solution_order=False)
    t0 = time.time()
    lom.sim.run(name=name, components=['pair'],
                open_terminations=[('pair', f'{e}_{t}') for e in 'ab' for t in ('top', 'bot')], box_plus_buffer=False)
    cm = lom.sim.capacitance_matrix
    cm.to_csv(out/'matrix.csv')
    if str(lom.sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise RuntimeError(f'Unexpected units {lom.sim.units}')
    sig = sorted([n for n in cm.index if 'ground' not in n.lower()], key=lambda n: 0 if 'top' in n else 1)
    assert len(sig) == 2, sig
    m = cm.loc[sig, sig].to_numpy(float)
    passes = []
    for k, mat in lom.sim.capacitance_all_passes.items():
        frame = mat if isinstance(mat, pd.DataFrame) else pd.DataFrame(mat, index=cm.index, columns=cm.columns)
        frame.to_csv(out/f'pass_{k}.csv')
        passes.append(frame.loc[sig, sig].to_numpy(float))
    delta = float(np.max(abs(passes[-1]-passes[-2])/np.maximum(abs(passes[-1]), 1e-30))) if len(passes) > 1 else None
    result = dict(units=str(lom.sim.units), signals=sig, length_mm=length_mm, spacing_mm=spacing_mm,
        maxwell_fF=m.tolist(), mutual_fF=float(-m[0, 1]), ground_fF=m.sum(axis=1).tolist(),
        max_passes=max_passes, percent_error=percent_error, last_pass_max_relative_change=delta,
        solver_converged=bool(lom.sim.is_converged), minutes=(time.time()-t0)/60,
        aedt_design=lom.sim.renderer.pinfo.design.name)
    write(out/'caps.json', result)
    lom.sim.renderer.pinfo.project._project.Save()
    print(json.dumps(result), flush=True)
    return result


def coupled_model(f, cap, hammer, lengths, z0, velocity, feed_lengths, spacing, cm_per_mm, sections):
    """Nodal network as hammer_fixed_s21.model (centre-plane hammer), with coupled sections.

    sections: lengths (mm) of the parallel run, ordered from the filter's cp end:
      below = [filter cp-side uncoupled, coupled run below the tee], tee half, tee half, above = coupled run.
    """
    f = np.asarray(f, float)
    w = 2*np.pi*f*1e9
    below_free, below, above = sections['filter_uncoupled_mm'], sections['below_mm'], sections['above_mm']
    th = P.TEE_HALF
    # 0,1 ports; 2,3 cin; 4 n0; 5,6 cp; 7 filter tap; 8 res tap; 9,10 hammer centres;
    # 11 F(y0) 12 F(tee bottom) 13 F(tee top) 14 F(above end); 15 R(y0) 16 R(tee bottom) 17 R(tee top) 18 R(above end)
    n = 19
    Y = np.zeros((len(f), n, n), complex)
    cself = 1/(z0*velocity)                         # F/m, conformal value of the existing model
    cm = cm_per_mm*1e-15/1e-3                       # F/m
    Yc = velocity*np.array([[cself+cm, -cm], [-cm, cself+cm]])

    def tl(a, b, length):
        t = w/velocity*length*1e-3
        dg, off = -1j/z0/np.tan(t), 1j/z0/np.sin(t)
        Y[:, a, a] += dg; Y[:, b, b] += dg; Y[:, a, b] += off; Y[:, b, a] += off

    def ctl(a, b, length):                          # a = (filter, res) start nodes, b = end nodes
        t = (w/velocity*length*1e-3)[:, None, None]
        dg, off = -1j*Yc/np.tan(t), 1j*Yc/np.sin(t)
        ia, ib = np.ix_(a, a), np.ix_(b, b)
        Y[(slice(None),)+ia] += dg; Y[(slice(None),)+ib] += dg
        Y[(slice(None),)+np.ix_(a, b)] += off; Y[(slice(None),)+np.ix_(b, a)] += off

    def pi(a, b, values):
        c, ga, gb = np.asarray(values)*1e-15
        Y[:, a, a] += 1j*w*(c+ga); Y[:, b, b] += 1j*w*(c+gb)
        Y[:, a, b] -= 1j*w*c; Y[:, b, a] -= 1j*w*c

    def stub(node, length, shorted):
        t = w/velocity*length*1e-3
        Y[:, node, node] += -1j/z0/np.tan(t) if shorted else 1j/z0*np.tan(t)

    tl(0, 2, feed_lengths[0]); pi(2, 3, cap); tl(3, 4, feed_lengths[1]); tl(4, 1, feed_lengths[2])
    tl(4, 5, feed_lengths[3]); pi(5, 6, cap)
    tl(6, 11, below_free)
    ctl([11, 15], [12, 16], below)
    ctl([12, 16], [7, 8], th)
    ctl([7, 8], [13, 17], th)
    ctl([13, 17], [14, 18], above)
    stub(14, lengths[0]-P.FILTER_TAP-th-above, shorted=True)        # filter tail to its shorted end
    stub(18, P.RES_TAP-th-above, shorted=False)                       # readout open end
    stub(15, lengths[1]-P.RES_TAP-th-below, shorted=True)             # readout tail to its shorted end
    tl(7, 9, spacing/2); pi(9, 10, hammer); tl(10, 8, spacing/2)
    port = Y[:, :2, :2] - Y[:, :2, 2:] @ np.linalg.solve(Y[:, 2:, 2:], Y[:, 2:, :2])
    eye = np.broadcast_to(np.eye(2), port.shape)
    S = np.linalg.solve(eye+50*port, eye-50*port)
    return S[:, 0, 0], S[:, 1, 0]


def parallel_sections(geometry):
    """Overlap of filter (x=0) and readout (x=spacing) vertical runs, from the solved full geometry."""
    from shapely import wkt
    ids = {n: i for i, n in enumerate(geometry['options'], start=1)}
    lines = {r['component']: wkt.loads(r['wkt']) for r in geometry['path'] if not r['subtract']}
    feed = list(lines[ids['feed_jf']].coords); tail = list(lines[ids['filter_tail']].coords)
    rtail = list(lines[ids['res_tail']].coords); ropen = list(lines[ids['res_open_run']].coords)
    tee_lo, tee_hi = P.TAP_Y-P.TEE_HALF, P.TAP_Y+P.TEE_HALF
    assert np.isclose(feed[-1][1], tee_lo) and np.isclose(rtail[0][1], tee_lo)
    assert np.isclose(tail[0][1], tee_hi) and np.isclose(ropen[0][1], tee_hi)
    below_start = max(feed[0][1], rtail[1][1])                       # readout tail turns at its first vertex
    above_end = min(tail[1][1], ropen[1][1])                         # filter tail turns at its first vertex
    below, above = tee_lo-below_start, above_end-tee_hi
    return dict(below_mm=below, above_mm=above, filter_uncoupled_mm=P.FILTER_TAP-P.TEE_HALF-below,
                note='Straight-segment overlaps between corner vertices; fillet arc length differences ignored')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('hfss_run', type=Path, help='hammer_fixed_s21 run with comparison.json and raw/S_fine.csv')
    ap.add_argument('geometry_run', type=Path, help='candidate run with raw/full_geometry.json')
    ap.add_argument('--coupon-lengths-mm', type=float, nargs=2, default=[1., 2.])
    ap.add_argument('--upper-ghz', type=float, default=6.25)
    ap.add_argument('--reuse-pair', type=Path, help='Existing run whose raw/pair_* coupons are reused')
    args = ap.parse_args()
    hfss, geo_run = args.hfss_run.resolve(), args.geometry_run.resolve()
    metrics = json.loads((hfss/'comparison.json').read_text())
    geometry = json.loads((geo_run/'raw/full_geometry.json').read_text())
    spacing = metrics['spacing_mm']
    sections = parallel_sections(geometry)
    folder = create_run('Design4', 'resonator', 'hammer_coupled_model', entrypoint=__file__,
        inputs=[hfss/'comparison.json', hfss/'raw/S_fine.csv', geo_run/'raw/full_geometry.json'],
        upstream=[str(hfss), str(geo_run)] + ([str(args.reuse_pair.resolve())] if args.reuse_pair else []),
        parameters=dict(EM_rerun_HFSS=False, Q3D_pair_coupons_mm=args.coupon_lengths_mm, spacing_mm=spacing,
            sections=sections, Q3D_max_passes=20, Q3D_requested_percent_error=.1, circuit_parameters_fitted=False,
            assumption='quasi-TEM homogeneous coupled lines, L = C^-1/v^2; self C conformal; mutual C from Q3D difference'))
    print('RUN_FOLDER='+str(folder), flush=True)
    try:
        pairs = []
        for L in args.coupon_lengths_mm:
            out = folder/f'raw/pair_{L:g}mm'
            if args.reuse_pair:
                import shutil
                shutil.copytree(args.reuse_pair.resolve()/f'raw/pair_{L:g}mm', out)
                pairs.append(json.loads((out/'caps.json').read_text()))
            else:
                pairs.append(solve_pair(L, spacing, out, f'D4PAIR_{folder.name[-8:]}_{L:g}mm'))
        dL = args.coupon_lengths_mm[1]-args.coupon_lengths_mm[0]
        cm_per_mm = (pairs[1]['mutual_fF']-pairs[0]['mutual_fF'])/dL
        cg_per_mm = [(b-a)/dL for a, b in zip(pairs[0]['ground_fF'], pairs[1]['ground_fF'])]
        z0, v = metrics['Z0_ohm'], metrics['velocity_m_s']
        cap, feeds, lengths = metrics['cin_cp_fF'], metrics['feed_lengths_mm'], metrics['line_lengths_mm']
        hammer = metrics['models']['deembedded']['hammer_fF']
        mult = metrics['port_polarity']['S21_multiplier']
        data = pd.read_csv(hfss/'raw/S_fine.csv', index_col=0)
        data = data.loc[data.index <= args.upper_ghz]
        f = data.index.to_numpy(float)
        hv = np.column_stack([data[k].map(complex).to_numpy() for k in ('S11', 'S21')])
        # Regression guard: no mutual -> identical to the existing deembedded model
        a0, b0 = F.model(f, cap, hammer, lengths, z0, v, feeds, spacing=spacing, cap_half_span=0., omit_hammer_stubs=False)
        a1, b1 = coupled_model(f, cap, hammer, lengths, z0, v, feeds, spacing, 0., sections)
        assert np.allclose(b0, b1, atol=1e-9) and np.allclose(a0, a1, atol=1e-9), 'Coupled model with Cm=0 must equal old model'
        observed, observed_curve = fit_poles(f, hv, 3)
        cases = {'deembedded': (hammer, 0.), 'coupled_lines': (hammer, cm_per_mm), 'coupled_lines_no_hammer': ([0., 0., 0.], cm_per_mm)}
        export = pd.DataFrame({'frequency_GHz': f, 'HFSS_S21': hv[:, 1], 'HFSS_pole_fit_S21': observed_curve[:, 1]})
        models = {}
        for name, (ham, cmv) in cases.items():
            a, b = coupled_model(f, cap, ham, lengths, z0, v, feeds, spacing, cmv, sections)
            b = b*mult
            export[name+'_S21'] = b
            entry = dict(hammer_fF=list(map(float, ham)), mutual_per_mm_fF=cmv,
                         rms_complex_S21=float(np.sqrt(np.mean(abs(b-hv[:, 1])**2))))
            try:
                pf, _ = fit_poles(f, np.column_stack([a, b]), 3)
                entry.update(poles=pf, J_relative_error=pf['J_MHz']/observed['J_MHz']-1 if pf['J_MHz'] else None,
                    bare_fp_error_MHz=(pf['bare_fp_GHz']-observed['bare_fp_GHz'])*1000,
                    bare_fr_error_MHz=(pf['bare_fr_GHz']-observed['bare_fr_GHz'])*1000)
            except Exception as e:
                entry['pole_fit_error'] = repr(e)
            models[name] = entry
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        style = {'deembedded': dict(label='Q3D circuit - deembedded hammer only', color='tab:orange', linestyle='-.'),
                 'coupled_lines': dict(label="Q3D circuit - hammer + coupled lines (C'm)", color='tab:blue', linestyle='--'),
                 'coupled_lines_no_hammer': dict(label="Q3D circuit - coupled lines only, no hammer", color='tab:green', linestyle=':')}
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        db = lambda z: 20*np.log10(np.maximum(abs(z), 1e-12))
        ax[0].plot(f, db(hv[:, 1]), color='black', lw=2, label='HFSS Driven S21 - full-wave')
        ax[0].plot(f, db(observed_curve[:, 1]), color='tab:red', lw=1, label=f"HFSS two-pole fit (J {observed['J_MHz']:.3f} MHz)")
        for name in cases:
            ax[0].plot(f, db(export[name+'_S21']), lw=1.8, **style[name])
            ax[1].plot(f, abs(export[name+'_S21']-hv[:, 1]), lw=1.8, **style[name])
        for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
        ax[0].set_ylabel('S21 (dB)'); ax[1].set_ylabel('|Q3D circuit - HFSS Driven|'); ax[1].set_xlabel('Frequency (GHz)')
        fig.suptitle(f"Design4 h30 / 200 um: Q3D circuit with distributed coupling, C'm = {cm_per_mm:.3f} fF/mm\n"
                     'Q3D curves are circuit predictions from extracted capacitances (no fitting)')
        fig.tight_layout(); fig.savefig(folder/'comparison.png', dpi=170); plt.close(fig)
        export.to_csv(folder/'model_comparison.csv', index=False)
        result = dict(mutual_per_mm_fF=cm_per_mm, ground_per_mm_fF=cg_per_mm,
            conformal_self_per_mm_fF=1e12/(z0*v), pair_coupons=pairs, sections=sections,
            HFSS_poles=observed, models=models,
            limitations=['Homogeneous quasi-TEM coupled lines; even/odd velocity difference neglected.',
                'Ground-strip interruption by the hammer pocket and bridge loading are not modelled; tee region uses the straight-pair C\'m.',
                'Cin/Cp inputs are the baseline values whose Q3D convergence target is unmet.',
                'Pole inversion is conditional on a two-mode model; HFSS Delta S below 0.01 on one pass only.'])
        write(folder/'comparison.json', result)
        finish_run(folder, status='completed', summary=result)
        print(json.dumps({k: result[k] for k in ('mutual_per_mm_fF', 'ground_per_mm_fF', 'sections')}), flush=True)
        for k, m in models.items():
            print(k, {x: m.get(x) for x in ('rms_complex_S21', 'J_relative_error', 'bare_fp_error_MHz', 'bare_fr_error_MHz')},
                  m.get('poles', {}).get('J_MHz'), m.get('pole_fit_error'))
    except Exception:
        import traceback
        err = traceback.format_exc(); (folder/'error.log').write_text(err)
        finish_run(folder, status='failed', summary={'error': err}); raise


if __name__ == '__main__':
    main()
