"""Hybrid EM-circuit model, EM part: HFSS driven 3-port of the C_j branch block only.

Block = feed_jf (filter, from C_p pin b down to jf) + jf + stub_f + cj + stub_r + jr + res_a + open_res
(filter / resonator straight segments around C_j, the resonator open end, both T junctions and C_j).
Ports (50 Ohm, reference planes at the pins):
  1  feed_jf.start   filter top    = C_p pin b
  2  jf.prime_start  filter bottom -> 'filt' route to the shorted end
  3  jr.prime_end    resonator top -> 'res' route to the shorted end
Everything else (C_in, C_p, feed lines, long filter / resonator lines) stays in the circuit model (theory_abcd).

Usage: python hfss_block6.py <tag> [f_start f_stop count sweep_type]   (env SOL_F, MAXP, BASIS, RUN_NOTE)
Output: results/blocks/<tag>_S3.csv (complex S11..S33), <tag>_meta.json (block geometry, remaining lengths), convergence.
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
import design6_layout as d6, ftf_layout as f3

BLOCK = ['feed_jf', 'jf', 'stub_f', 'cj', 'stub_r', 'jr', 'res_a', 'open_res']
PORTS = [('feed_jf', 'start', 50), ('jf', 'prime_start', 50), ('jr', 'prime_end', 50)]
OUT = os.path.join(HERE, 'results', 'blocks')


def block_meta(d):
    """Remaining (circuit) line lengths outside the block, from the drawn routes."""
    return dict(Lf_rest=f3.path_length(d, 'filt'), Lr_rest=f3.path_length(d, 'res'),
                CJ=dict(head_w_um=d6.CJ_HW, head_h_um=d6.CJ_HH, gap_um=d6.CJ_HG), D_CJ=d6.D_CJ, DX=0.20,
                L_FILTER=d6.L_FILTER, L_RES_TOTAL=d6.L_RES_TOTAL)


if __name__ == '__main__':
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    tag = sys.argv[1] if len(sys.argv) > 1 else 'blk_T5'
    f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (5.8, 7.6)
    count = int(sys.argv[4]) if len(sys.argv) > 4 else int(round((f1 - f0) * 1000)) + 1
    stype = sys.argv[5] if len(sys.argv) > 5 else 'Fast'
    os.makedirs(OUT, exist_ok=True)
    d = d6.build_design()
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in dict(freq_ghz=float(os.environ.get('SOL_F', 6.7)), max_delta_s=0.002, max_passes=int(os.environ.get('MAXP', 20)),
                     min_passes=6, min_converged=2, pct_refinement=30, basis_order=int(os.environ.get('BASIS', 2))).items():
        sim.setup[k] = v
    sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=count, type=stype, save_fields=False)
    t0 = time.time()
    # ResEnd is pin-only: an end not listed here is SHORTED by the renderer -> the resonator open end must be listed
    sim.run(name=tag, components=BLOCK, open_terminations=[('res_a', 'start')], port_list=PORTS)
    try:
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(os.path.join(OUT, f'{tag}_convergence.csv')); print('CONVERGENCE:'); print(conv.tail(4).to_string())
    except Exception as e:
        print('convergence not read:', repr(e)[:150])
    names = [f'S{i}{j}' for i in (1, 2, 3) for j in (1, 2, 3)]
    _, _, S = sim.renderer.get_params(names)
    S.index.name = 'f_GHz'; S.to_csv(os.path.join(OUT, f'{tag}_S3.csv'))
    meta = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), note=os.environ.get('RUN_NOTE', ''), sweep=[f0, f1, count, stype],
                sol_f=float(os.environ.get('SOL_F', 6.7)), ports=[p[:2] for p in PORTS], **block_meta(d))
    json.dump(meta, open(os.path.join(OUT, f'{tag}_meta.json'), 'w'), indent=1)
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) === saved', os.path.join(OUT, f'{tag}_S3.csv'))
    try: sim.close()
    except Exception: pass
