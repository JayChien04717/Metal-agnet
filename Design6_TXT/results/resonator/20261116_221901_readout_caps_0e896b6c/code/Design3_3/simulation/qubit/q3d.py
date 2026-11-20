"""Q3D of Design3_3 (straight ㄈ brackets, layout.build_design), with the Design3e procedure.

Same solver settings, rendered components and post-processing as Design3e/run_q3d.py (so the
numbers compare with Design3e e4); only the TC electrode geometry differs.

    python -B Design3_3/q3d.py arms6p7_18p9 --length 6.7 18.9
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import json
from pathlib import Path
import sys

HERE = _P('Design3_3')
sys.path[:0] = [_I('Design3_3')]
import run_q3d as e      # copy of the Design3e runner (verify_layout, Design3_2 base)  # noqa: E402
from layout import build_design, ARM_LENGTHS  # noqa: E402

base = e.base


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('--project', type=Path, default=Path.home() / 'AppData/Local/Temp/D33Tune/D33Tune.aedt')
    parser.add_argument('--length', type=float, nargs=2, default=ARM_LENGTHS,
                        metavar=('FX1', 'FX2'), help='Upper and lower arm extension per side, in µm.')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--passes', type=int, default=12)
    parser.add_argument('--min-passes', type=int, default=4)
    parser.add_argument('--error', type=float, default=.65)
    parser.add_argument('--order', choices=('Low', 'Medium', 'High'), default='Medium')
    args = parser.parse_args()

    def build(**kw):
        # the Q3D runner passes fx1_options=None etc.; build_design fills its defaults into dicts
        return build_design(rounded=True, arm_lengths=args.length,
                            **{k: v for k, v in kw.items() if v is not None})

    e.verify_layout(build(for_process=True, fx_options={'array_end_metal': 'False'}))
    if args.check:
        return
    base.RESULTS = base.q3d.RESULTS = _P('Design3_3/results/general/legacy')
    base.RESULTS.mkdir(exist_ok=True)
    base.save_aedt.use_project(args.project)
    base.q3d.build_design = build

    class Runner(base.q3d.Q3DRunner):
        def _prepare(self, design):
            lom = super()._prepare(design)
            base.join_squid_ground(lom.sim.renderer, merge_readouts=not args.local)
            renderer = lom.sim.renderer

            def analyze(name):
                setup = renderer.pinfo.get_setup(name)
                setup.parent._setup_module.EditSetup(name, [
                    f'NAME:{name}', 'AdaptiveFreq:=', '5GHz', 'SaveFields:=', False, 'Enabled:=', True,
                    ['NAME:Cap', 'MaxPass:=', args.passes, 'MinPass:=', args.min_passes,
                     'MinConvPass:=', 2, 'PerError:=', args.error, 'PerRefine:=', 25,
                     'AutoIncreaseSolutionOrder:=', False, 'SolutionOrder:=', args.order,
                     'Solver Type:=', 'Iterative']])
                try:
                    setup.analyze(name)
                except Exception:
                    desktop = base.save_aedt._desktop()._desktop
                    print(desktop.GetMessages(args.project.stem, f'{args.tag}_q3d', 0), flush=True)
                    raise
            renderer.analyze_setup = analyze
            run = lom.sim.run
            lom.sim.run = lambda **kw: run(name=args.tag, **kw)
            return lom

    base.q3d.OPEN_PINS = [p for p in base.q3d.OPEN_PINS if p != ('TC', 'readout_line')]
    if not args.local:
        base.q3d.COMPONENTS = [*base.q3d.COMPONENTS, 'res_FX1', 'res_FX2']
        base.q3d.OPEN_PINS = [p for p in base.q3d.OPEN_PINS if p[1] != 'readout_line']
        base.q3d.OPEN_PINS += [('res_FX1', 'start'), ('res_FX2', 'start')]
    runner = Runner(setup=dict(max_passes=args.passes, min_passes=args.min_passes, percent_error=args.error))
    try:
        runner.run(args.tag, with_resonators=False)
        info = dict(vars(args), arm_style='bracket', gap=[4, 4])
        (base.RESULTS / f'{args.tag}_geometry.json').write_text(json.dumps(info, indent=2, default=str))
        conv, _ = runner.lom.sim.renderer.pinfo.setup.get_convergence()
        conv.to_csv(base.RESULTS / f'{args.tag}_convergence.csv')
    finally:
        base.save()


if __name__ == '__main__':
    main()
