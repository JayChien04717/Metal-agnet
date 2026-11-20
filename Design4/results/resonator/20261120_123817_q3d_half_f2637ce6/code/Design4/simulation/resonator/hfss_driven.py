"""Design4 HFSS driven S21 of one stub variant (open fixture, no qubit).

Geometry: Design4/layout/build.py build(variant) with the filter / resonator lengths that the frequency model
(purcell_readout_model.length_for_f) gives for f_target from the Q3D loads of the variant's q3d_module run
(C_j, head excess, C_p mutual + plate ground). Nothing is fitted.

    python -I -B Design4/simulation/resonator/hfss_driven.py --variant ref|half --q3d <q3d run folder> [--passes 6]

Solver (skill.md sections 12-13, shared purcell_resonator_pair_sim.run_driven): first order, 10 um length mesh on
every metal object, qiskit auto bonds every 350 um plus PEC straps at PurcellReadout.strap_points(), two 50 Ohm
lumped feed ports, Fast wide sweep f0 +- 0.5 GHz then a 50 kHz Fast fine sweep around the dips. The pass budget
(default 6, as Design8's matching run) is an execution setting, not convergence evidence.
Analysis: simulation/resonator/compare.py.
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design4' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)

import build as B  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402


def q3d_loads(q3d_run, variant):
    s = json.loads((Path(q3d_run) / 'summary.json').read_text(encoding='utf-8'))
    if s.get('variant') != variant:
        raise ValueError(f'{q3d_run} is variant {s.get("variant")!r}, not {variant!r}')
    d = s['derived']
    return {k: float(d[k]) for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground', 'cp_frame_ground')}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variant', required=True, choices=sorted(B.GEOMETRY['variants']))
    ap.add_argument('--q3d', type=Path, required=True)
    ap.add_argument('--passes', type=int, default=6)
    args = ap.parse_args()
    v, q3d = args.variant, args.q3d.resolve()
    loads = q3d_loads(q3d, v)
    lengths = B.model_lengths(v, caps={k: loads[k] for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground')})
    f0 = B.GEOMETRY['f_target_GHz']
    params = dict(variant=v, line_spacing_um=B.spacing_um(v), stub_metal_um=B.stub_metal_um(v), lengths=lengths,
                  q3d_loads=loads, passes=args.passes, max_delta_s=0.02, band_GHz=[f0 - 0.5, f0 + 0.5],
                  mesh=S.MESH, wire_bonds=S.WB, straps='PurcellReadout.strap_points, PEC, 10 um wire', ports=B.PORTS,
                  geometry='Design4 build(variant) open fixture', geometry_revision=B.GEOMETRY['revision'])
    folder = create_run('Design4', 'resonator', f'hfss_driven_{v}', entrypoint=__file__,
                        inputs=[ROOT / 'Design4/parameters/geometry.json', q3d / 'summary.json'], upstream=[str(q3d)],
                        parameters=params)
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import save_aedt
        import hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design4' / 'simulation'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        S.write(folder / 'run.json', rec)
        d, info = B.build(v, lengths=lengths)
        (folder / 'raw').mkdir(exist_ok=True)
        S.write(folder / 'raw' / 'geometry.json', S.geometry_dump(d))
        S.write(folder / 'build_info.json', info)
        comps = [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']
        raw = S.run_driven(folder, d, comps, B.PORTS, [], f0, (f0 - 0.5, f0 + 0.5), passes=args.passes,
                           tag=f'S{v}')
        print(json.dumps(raw, default=float), flush=True)
        finish_run(folder, status='completed', summary=dict(variant=v, driven=raw, lengths=lengths,
                   validation='fixed pass budget; last Delta S recorded, not claimed converged'))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
