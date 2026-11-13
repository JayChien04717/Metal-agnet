"""HFSS driven modal: mutual inductance flux line -> fluxonium loop.

Design2: renders Q1 alone (FluxoniumPocket with the simulation loop, u_array='True', same as
Design1 / 3): the two U sides are lumped-L boundaries (l_inductance / 2 each), the JJ gap in the
middle of the lead line is Port 2 (port only) and the flux bias line pin is Port 1 (lumped port
across its cpw_gap endcap; the other leg is shorted to the pocket wall). Readout, charge and fake
flux line pins are open-terminated.

    M = Im Z21 / (2 pi f)        (quasi-static, far below any resonance)

Also printed as checks: L_fl = Im Z11 / w (flux U-loop self inductance) and
L_loop = Im Z22 / w (array + geometric loop inductance, ~ l_inductance).
Method: ../Design1/notes/mutual_inductance_method.md.

Usage:  python hfss_mutual_fx.py [tag]
or:     from hfss_mutual_fx import run
        run('fl10', fx_options=dict(flux_line_options=dict(fl_sep='10um')))
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal import designs, Dict  # noqa: E402
from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim  # noqa: E402
from ftf_layout import FX1_OPTIONS, Q1, CPW_WIDTH, CPW_GAP  # noqa: E402
from fluxonium_withJJ import FluxoniumPocket  # noqa: E402
from q3d_ftf import RESULTS  # noqa: E402

PHI0 = 2.067833848e-15  # Wb

QUBIT = 'Q1'
PORT_LIST = [(QUBIT, 'flux_bias_line', 50)]
JJ_TO_PORT = [(QUBIT, 'rect_jj', 50, False)]     # port only: open-loop EMF
OPEN_PINS = [(QUBIT, 'readout_line'), (QUBIT, 'charge_line'), (QUBIT, 'fake_flux_bias_line')]

# adaptive mesh at 100 MHz. Delta S is a weak criterion here: max_delta_s 0.002 stopped
# after 4 passes (7k elements) with M 7 % low, so at least 12 passes are forced
# (12 -> 15 passes: M changes 0.02 %, results/runs/mutual_fx_conv vs mutual_fx_p15).
# Discrete sweep 10-100 MHz: M must be flat over it (quasi-static check).
SETUP = dict(freq_ghz=0.1, max_delta_s=0.0002, max_passes=15, min_passes=12, min_converged=3,
             pct_refinement=30, basis_order=1,
             sweep_setup=Dict(name='Sweep', start_ghz=0.01, stop_ghz=0.1, count=10,
                              step_ghz=None, type='Discrete', save_fields=False))


def build_fx_design(fx_options=None):
    """Q1 of Design2 alone (FluxoniumPocket, simulation loop u_array='True') at its chip position."""
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.chips.main.size.size_z = '-525um'
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    opt = dict(FX1_OPTIONS, u_array='True')
    for k, v in (fx_options or {}).items():
        base = opt.get(k, FluxoniumPocket.default_options.get(k))
        opt[k] = Dict(base, **v) if isinstance(v, dict) and isinstance(base, dict) else v
    FluxoniumPocket(design, QUBIT, options=dict(Q1, **opt))
    return design


def _port_index(sim, z):
    """0-based index of the JJ port in the Z matrix: the port with the larger Im Zii at f_min
    (the JJ port sees the array inductance, ~200 nH; the flux port only the short
    U-loop, ~pH). The GetExcitations order does not match the exported network-data
    order (first run: excitations [rect_jj, flux_line], matrix [flux_line, rect_jj]),
    so the names are only stored for reference."""
    try:
        exc = list(sim.renderer.pinfo.design._boundaries.GetExcitations())
        names = [e for e in exc if isinstance(e, str) and e.startswith('LumpPort')]
    except Exception:
        names = None
    # lowest frequency only: above the loop resonance (~2.7 GHz) Im Zjj turns negative
    return (0 if z['Z11'][0].imag > z['Z22'][0].imag else 1), names


def run(tag='mutual_fx', fx_options=None, setup=None):
    """Design2 Q1 (FluxoniumPocket u_array, FX1_OPTIONS + fx_options)."""
    return run_design(tag, build_fx_design(fx_options), QUBIT, 'flux_bias_line', 'rect_jj',
                      OPEN_PINS, setup=setup, info=dict(fx_options=fx_options))


def run_design(tag, design, qubit, flux_pin, jj, open_pins=(), setup=None, info=None):
    """M between pin `flux_pin` (Port 1) and junction `jj` (Port 2, port only) of
    component `qubit` in `design`; other junctions keep their hfss_inductance."""
    sim = ScatteringImpedanceSim(design, 'hfss')
    for k, v in dict(SETUP, **(setup or {})).items():
        sim.setup[k] = v
    t0 = time.time()
    sim.run_sim(name=tag, components=[qubit], open_terminations=list(open_pins),
                port_list=[(qubit, flux_pin, 50)], jj_to_port=[(qubit, jj, 50, False)],
                box_plus_buffer=True)

    names = ['Z11', 'Z21', 'Z12', 'Z22']
    freqs, curves, _ = sim.renderer.get_params(names)
    z = dict(zip(names, (np.asarray(c) for c in curves)))
    j, ports = _port_index(sim, z)
    fl = 1 - j
    w = 2 * np.pi * np.asarray(freqs)
    zf = z[f'Z{fl + 1}{fl + 1}']
    zj = z[f'Z{j + 1}{j + 1}']
    zm = (z['Z21'] + z['Z12']) / 2               # reciprocal; average the two
    df = pd.DataFrame({
        'f_MHz': np.asarray(freqs) / 1e6,
        'M_pH': zm.imag / w * 1e12,
        'L_fl_pH': zf.imag / w * 1e12,
        'L_loop_nH': zj.imag / w * 1e9,
        'ReZ21_ohm': zm.real, 'ImZ21_ohm': zm.imag,
        'Z21_Z12_mismatch': np.abs(z['Z21'] - z['Z12']) / np.abs(zm),
    })
    df.to_csv(os.path.join(RESULTS, f'{tag}_Z.csv'), index=False)
    try:
        # adaptive-pass table (delta S per pass); get_convergences() also builds an
        # eigenmode frequency report, so read the setup directly
        conv_t, _ = sim.renderer.pinfo.setup.get_convergence()
        conv_t.to_csv(os.path.join(RESULTS, f'{tag}_convergence.csv'))
    except Exception as e:
        print('convergence table not saved:', e)

    M = df['M_pH'].to_numpy()
    summary = dict(tag=tag, M_pH_fmin=float(M[0]), M_pH_mean=float(np.mean(M)),
                   M_spread_pct=float(np.ptp(M) / abs(np.mean(M)) * 100),
                   I_phi0_mA=float(PHI0 / abs(M[0] * 1e-12) * 1e3),
                   L_fl_pH=float(df['L_fl_pH'].mean()), L_loop_nH=float(df['L_loop_nH'].mean()),
                   excitations=ports, jj_port_in_matrix=j + 1, **(info or {}),
                   setup={k: v for k, v in sim.setup.items() if k != 'sweep_setup'},
                   minutes=round((time.time() - t0) / 60, 1))
    with open(os.path.join(RESULTS, f'{tag}_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(df.to_string(index=False))
    print(f"\nM(f_min) = {M[0]:.4f} pH, mean M = {summary['M_pH_mean']:.3f} pH (spread {summary['M_spread_pct']:.1f} % over the sweep), "
          f"I_phi0 = {summary['I_phi0_mA']:.3f} mA, L_loop = {summary['L_loop_nH']:.1f} nH, "
          f"{summary['minutes']} min")
    try:
        sim.close()
    except Exception:
        pass
    return df, summary


if __name__ == '__main__':
    run(*(sys.argv[1:2] or ['mutual_fx']))
