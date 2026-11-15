"""Design4 hammer CPW-end sensitivity experiment; all outputs belong to new runs.

Reference plane: hammer centre. Subtract half the capacitance of a continuous
CPW with the same two outer endpoints. This defines an EXCESS shunt capacitance,
not the absolute pad-to-ground capacitance. The channel model must retain TLs
from each tee to the hammer centre. No line-length or capacitor fitting.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import os
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / p) for p in ('code/general', 'Design4/layout',
    'Design4/simulation/resonator', 'Design5/simulation/resonator', 'pattern design')]
import numpy as np
import pandas as pd
import purcell_6500 as P
from qiskit_metal import designs, draw, Dict
from qiskit_metal.qlibrary.core import QComponent
from hammer_cap import HammerCap
from run_record import create_run, finish_run


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str) + '\n', encoding='utf-8')


class ReferenceCPW(QComponent):
    default_options = Dict(half_length='100um', trace_width='16um', trace_gap='9um')
    def make(self):
        p = self.p
        line = draw.LineString([(-p.half_length, 0), (p.half_length, 0)])
        self.add_qgeometry('path', {'trace': line}, width=p.trace_width, layer=1)
        self.add_qgeometry('path', {'cut': line}, width=p.trace_width+2*p.trace_gap, subtract=True, layer=1)
        for name, sign in [('a', -1), ('b', 1)]:
            self.add_pin(name, points=np.array([[sign*(p.half_length-p.trace_width), 0],
                [sign*p.half_length, 0]]), width=p.trace_width, gap=p.trace_gap, input_as_norm=True)


def coupon(length_um, reference=False):
    d = designs.DesignPlanar()
    d.chips.main.material = P.d5.SI_MATERIAL
    # Fixed domain for every sweep point, rather than auto-growing with the coupon.
    d.chips.main.size.update(size_x='5.2mm', size_y='4.3mm', size_z='-500um', center_x='0mm', center_y='0mm')
    h = P.PAPER_HAMMER
    if reference:
        ReferenceCPW(d, 'ref', options=dict(**P.CPW, half_length=f"{h['gap_um']/2+h['width_um']+length_um}um"))
        return d, 'ref'
    HammerCap(d, 'cj', options=dict(**P.CPW, head_width=f"{h['width_um']}um",
        head_height=f"{h['height_um']}um", head_gap=f"{h['gap_um']}um",
        ground_gap=f"{h['ground_um']}um", lead_length=f'{length_um}um'))
    return d, 'cj'


def geometry_dump(d):
    data = {}
    for kind in ('poly', 'path'):
        data[kind] = [dict(component=int(r.component), name=r.name, subtract=bool(r.subtract),
            wkt=r.geometry.wkt, width=getattr(r, 'width', None)) for r in d.qgeometry.tables[kind].itertuples()]
    data['options'] = {n: dict(c.options) for n, c in d.components.items()}
    return data


def solve_cap(d, component, out, name):
    from qiskit_metal.analyses.quantization import LOManalysis
    out.mkdir(parents=True, exist_ok=False)
    write(out/'geometry.json', geometry_dump(d))
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer.options['x_buffer_width_mm'] = 0
    lom.sim.renderer.options['y_buffer_width_mm'] = 0
    lom.sim.setup.update(max_passes=8, min_passes=2, percent_error=0.3, freq_ghz=6.5,
        solution_order='Medium', auto_increase_solution_order=False)
    started = time.time()
    lom.sim.run(name=name, components=[component], open_terminations=[(component, 'a'), (component, 'b')], box_plus_buffer=False)
    cm = lom.sim.capacitance_matrix
    cm.to_csv(out/'matrix.csv')
    units = lom.sim.units
    if str(units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise RuntimeError(f'Unexpected capacitance units: {units}')
    signals = [x for x in cm.index if 'ground' not in x.lower() and not x.startswith('g_wb')]
    if len(signals) != (1 if component == 'ref' else 2):
        raise RuntimeError(f'Unexpected nets: {signals}')
    # Keep electrode order explicit, not the arbitrary matrix order.
    if component != 'ref':
        signals.sort(key=lambda n: 0 if ('head_a' in n or 'lead_a' in n or 'frame' in n) else 1)
    sub = cm.loc[signals, signals].to_numpy(float)
    mutual = -sub[0, 1] if len(signals) == 2 else 0.
    grounds = sub.sum(axis=1)
    pass_values = []
    for k, matrix in lom.sim.capacitance_all_passes.items():
        frame = matrix if isinstance(matrix, pd.DataFrame) else pd.DataFrame(matrix, index=cm.index, columns=cm.columns)
        frame.to_csv(out/f'pass_{k}.csv')
        pass_values.append(frame.loc[signals, signals].to_numpy(float))
    delta = None
    if len(pass_values) > 1:
        # pyEPR all-pass exports use F even when the final matrix uses fF.
        delta = float(np.max(np.abs(pass_values[-1]-pass_values[-2])/np.maximum(np.abs(pass_values[-1]), 1e-30)))
    result = dict(units=str(units), signals=signals, mutual_fF=float(mutual), ground_fF=grounds.tolist(),
        last_pass_max_relative_change=delta, solver_converged=bool(lom.sim.is_converged), minutes=(time.time()-started)/60,
        aedt_design=lom.sim.renderer.pinfo.design.name, ground_bonds='none: finite local coupon, common electrostatic ground')
    write(out/'caps.json', result)
    lom.sim.renderer.pinfo.project._project.Save()
    print(json.dumps(result), flush=True)
    return result


def plot_sweep(rows, folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    table = pd.DataFrame(rows)
    table.to_csv(folder/'cpw_sweep.csv', index=False)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(table.length_um, table.mutual_fF, 'o-')
    ax[0].set_ylabel('Mutual capacitance (fF)')
    for key in ('raw_a_fF', 'raw_b_fF', 'excess_a_fF', 'excess_b_fF'):
        ax[1].plot(table.length_um, table[key], 'o-', label=key)
    for a in ax:
        a.set_xlabel('CPW lead length per side (um)'); a.grid(alpha=.3)
    ax[1].set_ylabel('Shunt capacitance (fF)'); ax[1].legend()
    fig.tight_layout(); fig.savefig(folder/'cpw_sweep.png', dpi=160); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lengths', default='50,100,200,400,800')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    lengths = [float(x) for x in args.lengths.split(',')]
    parameters = dict(lengths_um=lengths, hammer=P.PAPER_HAMMER, cpw=P.CPW,
        channel_lengths_mm=P.LENGTHS, spacing_mm=.4, trim_um=0, length_fitting=False,
        reference='continuous CPW between identical outer endpoints; half Cref removed per side',
        q3d_max_passes=8, q3d_percent_error=.3, domain_mm=[5.2,4.3,.5],
        substrate=P.d5.SI_MATERIAL, geometry_role='simulation', final_selection=False)
    folder = create_run('Design4', 'resonator', 'hammer_cpw_deembed', entrypoint=__file__,
        inputs=[ROOT/'Design4/design.json'], parameters=parameters)
    print('RUN_FOLDER='+str(folder), flush=True)
    try:
        for length in lengths:
            for ref in (False, True):
                d, c = coupon(length, ref)
                for table in d.qgeometry.tables.values():
                    if 'geometry' in table:
                        assert all(g.is_valid and not g.is_empty for g in table.geometry)
        if args.dry_run:
            finish_run(folder, status='completed', summary={'dry_run':True, 'geometry_checks':'passed', 'solver_run':False})
            return
        import save_aedt, hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT/'Design4/simulation'))
        manifest = json.loads((folder/'run.json').read_text(encoding='utf-8'))
        manifest['environment']['solver_version'] = save_aedt._desktop().get_version()
        manifest['geometry_revision'] = hashlib.sha256(json.dumps(parameters, sort_keys=True).encode()).hexdigest()
        write(folder/'run.json', manifest)
        rows = []
        for length in lengths:
            results = []
            for ref in (False, True):
                d, c = coupon(length, ref)
                results.append(solve_cap(d, c, folder/'raw'/f'{c}_{length:g}um',
                    f"D4DE_{folder.name[-8:]}_{c}_{length:g}"))
            cap, reference = results
            baseline = reference['ground_fF'][0]/2
            rows.append(dict(length_um=length, mutual_fF=cap['mutual_fF'],
                raw_a_fF=cap['ground_fF'][0], raw_b_fF=cap['ground_fF'][1], reference_half_fF=baseline,
                excess_a_fF=cap['ground_fF'][0]-baseline, excess_b_fF=cap['ground_fF'][1]-baseline))
            plot_sweep(rows, folder)
        finish_run(folder, status='completed', summary=dict(stage='cpw_length_sweep', rows=rows,
            remaining='fixed-length channel HFSS and before/after circuit comparison'))
    except Exception:
        (folder/'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error':traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
