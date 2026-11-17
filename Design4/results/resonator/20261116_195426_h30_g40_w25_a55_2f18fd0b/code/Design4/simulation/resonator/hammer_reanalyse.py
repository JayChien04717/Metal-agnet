"""Reanalyse a recorded baseline with the explicit omitted-stub diagnostic; no EM solve."""
from pathlib import Path
import argparse
import json
import shutil
import hammer_fixed_s21 as F
from run_record import create_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('capacitance_run', type=Path)
    args = parser.parse_args()
    baseline, upstream = args.baseline.resolve(), args.capacitance_run.resolve()
    source = json.loads((baseline/'run.json').read_text(encoding='utf-8'))
    F.P.PAPER_HAMMER.update(source['parameters']['hammer'])
    folder = create_run('Design4', 'resonator', 'hammer_model_audit', entrypoint=__file__,
        inputs=[baseline/'run.json', baseline/'raw/S_fine.csv', upstream/'cpw_sweep.csv'],
        upstream=[str(baseline), str(upstream)],
        parameters=dict(EM_rerun=False, diagnostic='omit hammer arms only', fit_parameters=[]))
    print('RUN_FOLDER='+str(folder), flush=True)
    shutil.copy2(baseline/'raw/S_fine.csv', folder/'raw/S_fine.csv')
    import hfss_hammer_straps as E
    E.SPACING_MM=.4
    E.GAP_UM=F.P.PAPER_HAMMER['gap_um']
    E.LENGTHS=F.P.LENGTHS
    E.NO_HAMMER=False
    design, _ = E.build(0)
    current=F.geometry_dump(design)
    previous=json.loads((baseline/'raw/channel_geometry.json').read_text())
    assert all(current[k]==previous[k] for k in ('poly','path')), 'Baseline geometry differs'
    def values(name):
        caps=json.loads((baseline/'raw'/name/'caps.json').read_text())
        return [caps['mutual_fF'], *caps['ground_fF']]
    metrics=F.compare(folder, upstream, values('cp_coupon'), values('original_hammer'), design)
    metrics['EM_rerun']=False
    metrics['interpretation']='Omitted-stub curve is a diagnostic, not an alternative physical layout'
    finish_run(folder, status='completed', summary=metrics)


if __name__ == '__main__':
    main()
