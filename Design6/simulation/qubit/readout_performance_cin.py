"""Estimate Design6 r10 dispersive shifts, relaxation, and conditional assignment with a series Cin."""
from argparse import ArgumentParser
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from scipy.constants import e, h, hbar, k as k_B
from scipy.integrate import cumulative_trapezoid
from scipy.special import erfc, erfcinv

import readout_qubit_r9 as Q

DESIGN = Q.ROOT / 'Design6'
sys.path.insert(0, str(DESIGN / 'simulation/resonator'))
from compare_loaded_s21_r9 import BASE_MODEL, read
from readout_network import port_Y, RES_TAP, TEE_HALF
from run_record import create_run, finish_run

CHANNELS = ('Q1', 'Q2')
QCAP = 10e5
XQP = 1e-6
TEFF_K = .060
AL_GAP_eV = 180e-6
ETA = (.1, .2, .4)
PHOTONS = (1, 5, 10)


def variant_g(baseline):
    saved = read(BASE_MODEL / 'compare.json')
    source = Path(saved['model_run']) / 'run.json'
    params = read(source)['parameters']
    settings = read(baseline / 'run.json')['parameters']
    geom = read(baseline / 'raw/geometry.json')
    indices = {name: index for index, name in enumerate(geom['options'], 1)}
    lines = {row['component']: row['wkt'] for row in geom['path'] if not row['subtract']}

    def y_coordinates(component):
        line = lines[indices[component]]
        assert line.startswith('LINESTRING (')
        return [float(vertex.split()[1]) for vertex in line[12:-1].split(',')]

    lengths = {q: (settings['lengths_mm'][q]['L_filter'], settings['lengths_mm'][q]['L_res'] + params['res_extra_mm'])
               for q in CHANNELS}
    sections = {}
    for q in CHANNELS:
        ro = y_coordinates('res_open_' + q)
        ft = y_coordinates('filter_tail_' + q)
        below = ro[0] - max(ft[1], ro[1])
        sections[q] = dict(saved['sections_mm'][q], below_mm=below,
                           res_open_free_mm=RES_TAP - TEE_HALF - below)
    c = saved['corrections']
    geo = dict(params['geometry'])
    geo.update({key: geo[key] * c['line_scale_F'] for key in ('z0_ohm', 'velocity_m_s')})
    _, cpa, cpb = params['cp_fF']
    cp = (c['cp_mutual_full_q3d_fF'], cpa, cpb - c['delta_used_fF'])
    correction_run = BASE_MODEL / 'run.json'
    return (lengths, geo, cp, params['hammer_deembedded_fF'], sections), dict(
        cm=read(correction_run)['parameters']['cm_pair_fF_per_mm'], feed_res=c['feed_res_fF'],
        extra_cj=c['extra_lumped_cj_fF']), source


def qubit_levels(p, q, fr, cap, levels, ncut):
    qubit = Q.C.scq.Transmon(EJ=p['EJsum_GHz'][q], EC=p['Ec_GHz'][q], ng=0, ncut=ncut)
    vals, vectors = qubit.eigensys(evals_count=levels)
    energy = vals - vals[0]
    charge = vectors.conj().T @ qubit.n_operator() @ vectors
    gn = cap['g/2pi [MHz]'] * fr / cap['f_r [GHz]'] / 1000
    n01 = abs(charge[0, 1])
    return dict(f01_GHz=energy[1], anharmonicity_MHz=1000 * (energy[2] - 2 * energy[1]), n01=n01,
                g01_MHz=1000 * gn * n01, g_charge_MHz=1000 * gn,
                critical_photons_01=((fr - energy[1]) / (2 * gn * n01))**2,
                critical_photons_12=((fr - energy[2] + energy[1]) / (2 * gn * abs(charge[1, 2])))**2,
                **Q.dispersive(energy, charge, fr, gn))


def incident_reflection(f, arm_mm, cin_fF, geo):
    if not cin_fF:
        return 0j
    w = 2 * np.pi * f * 1e9
    return np.exp(-2j * w * arm_mm * 1e-3 / geo['velocity_m_s']) / (1 + 2j * w * geo['z0_ohm'] * cin_fF * 1e-15)


def output_share(f, q, cin, geo):
    distance = cin['arm_from_Q2_tee_mm'] + (geo['feeds_mm']['feed_mid'] if q == 'Q1' else 0)
    gamma = incident_reflection(f, distance, cin['Cin_fF'], geo)
    return abs(1 + gamma)**2 / (2 * (1 + gamma.real))


def purcell(q, levels, beta, cin, args, kwargs):
    """Q3D charge lever arm times the vacuum voltage noise of the terminated readout network."""
    f = levels['f01_GHz']
    Y = port_Y(f, *args, **kwargs, cin_fF=cin['Cin_fF'], cin_arm_mm=cin['arm_from_Q2_tee_mm'], readout_probe=q)
    voltage = np.linalg.solve(Y + np.diag([1 / 50, 1 / 50, 0]), [0, 0, 1])
    z = voltage[2]
    port_resistances = abs(voltage[:2])**2 / 50
    np.testing.assert_allclose(z.real, sum(port_resistances), rtol=2e-6, atol=1e-12)
    omega = 2 * np.pi * f * 1e9
    gamma = 8 * e**2 * beta**2 * levels['n01']**2 * omega * z.real / hbar
    return dict(T1_us=1e6 / gamma, Re_Z_readout_ohm=z.real,
                Re_Z_per_port_ohm=dict(zip(('P8', 'P7'), port_resistances)))


def modal_matrix(fit, two_chi_MHz, state, drive_GHz):
    return np.array([[1000 * (fit['bare_fp_GHz'] - drive_GHz) - .5j * fit['kappa_p_MHz'], fit['J_MHz']],
                     [fit['J_MHz'], 1000 * (fit['bare_fr_GHz'] - drive_GHz) + state * two_chi_MHz]])


def driven_response(matrix, coupling, time_us):
    """Linear step response to one incident photon per microsecond."""
    generator = -2j * np.pi * matrix
    poles, vectors = np.linalg.eig(generator)
    weights = np.linalg.solve(vectors, [np.sqrt(coupling), 0.])
    return vectors @ (weights[:, None] * np.expm1(poles[:, None] * time_us) / poles[:, None])


def assignment(fit, q, two_chi_MHz, cin, geo, down_per_us, up_per_us, thermal):
    # Compare each constant-tone pulse at the same peak resonator occupancy across both states.
    time = np.linspace(0, 1, 1001)
    use = np.arange(20, len(time), 5)
    p_jump = -0.5 * (np.expm1(-down_per_us * time[use]) + np.expm1(-up_per_us * time[use]))
    rows = []
    for f in np.linspace(min(fit['centers_GHz']) - .012, max(fit['centers_GHz']) + .012, 151):
        share = output_share(f, q, cin, geo)
        kappa = 2 * np.pi * fit['kappa_p_MHz']
        a0, a1 = [driven_response(modal_matrix(fit, two_chi_MHz, state, f), kappa * (1 - share), time)
                  for state in (-.5, .5)]
        peak = np.maximum.accumulate(np.maximum(abs(a0[1])**2, abs(a1[1])**2))[use]
        information = cumulative_trapezoid(kappa * share * abs(a1[0] - a0[0])**2, time, initial=0)[use] / peak
        for photons in PHOTONS:
            for eta in ETA:
                noise_error = .5 * erfc(np.sqrt(eta * photons * information / 2))
                fidelity = (1 - noise_error) * (1 - p_jump)
                index = np.argmax(fidelity)
                rows.append(dict(n_peak=photons, eta_after_P8=eta, drive_GHz=f, integration_ns=1000 * time[use[index]],
                                 conditional_lower_bound=fidelity[index], noiseless_relaxation_error=p_jump[index],
                                 gaussian_error=noise_error[index],
                                 unreset_label_model=thermal + (1 - 2 * thermal) * fidelity[index],
                                 chip_input_dBm=10 * np.log10(photons / peak[index] * 1e6 * h * f * 1e9 / 1e-3)))
            valid = p_jump < .01
            allowed = (.01 - p_jump[valid]) / (1 - p_jump[valid])
            needed = 2 * erfcinv(2 * allowed)**2 / (photons * information[valid])
            index = np.argmin(needed)
            rows.append(dict(n_peak=photons, eta_for_99_lower_bound=needed[index], drive_GHz=f,
                             integration_ns=1000 * time[use[valid][index]]))
    examples = [max((r for r in rows if r['n_peak'] == n and r.get('eta_after_P8') == eta),
                    key=lambda r: r['conditional_lower_bound']) for n in PHOTONS for eta in ETA]
    thresholds = [min((r for r in rows if r['n_peak'] == n and 'eta_for_99_lower_bound' in r),
                      key=lambda r: r['eta_for_99_lower_bound']) for n in PHOTONS]
    return dict(examples=examples, thresholds=thresholds)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('comparison', type=Path)
    args = parser.parse_args()
    compare = args.comparison.resolve()
    manifest = read(compare / 'run.json')
    if manifest['status'] != 'completed':
        parser.error('The S21 comparison must be complete')
    baseline, candidate = map(Path, manifest['upstream'][:2])
    cfile = candidate / 'raw/capacitor.json'
    cin = read(cfile)
    fits = read(compare / 'summary.json')['variants']
    g_args, g_kw, model_source = variant_g(baseline)
    params = Q.CZ / 'params_r9_q3d.json'
    cmat = Q.Q3D / 'raw/cmat_fF.csv'
    core_summary = Q.Q3D / 'summary_out.json'
    folder = create_run('Design6', 'qubit', 'readout_performance_cin', entrypoint=__file__,
                        upstream=[str(p) for p in (compare, Q.Q3D, Q.CZ, BASE_MODEL)],
                        inputs=[compare / 'summary.json', cfile, baseline / 'run.json', baseline / 'raw/geometry.json',
                                params, cmat, core_summary, BASE_MODEL / 'compare.json', BASE_MODEL / 'run.json', model_source],
                        parameters=dict(candidate_only=True, Qcap=QCAP, xqp=XQP, T_eff_K=TEFF_K, Al_gap_eV=AL_GAP_eV,
                                        parent_geometry=manifest['parameters']['parent_geometry'], EJ='fixed layout L_J',
                                        charge_cutoffs=[15, 25], local_levels_and_photon_cutoffs=[6, 8],
                                        bath='T_eff detailed balance applied to the combined modeled downward rates',
                                        pulse='constant tone, integration 20–1000 ns, calibrated peak resonator occupancy',
                                        eta_after_P8=ETA, peak_photons=PHOTONS, Purcell='distributed variant G port impedance'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        p, core = Q.C.load_params(params), read(core_summary)
        matrix = pd.read_csv(cmat, index_col=0)
        report = {}
        for kind in ('baseline', 'candidate'):
            cap = dict(cin, Cin_fF=0.) if kind == 'baseline' else cin
            report[kind] = {}
            for q, fit in fits[kind]['channels'].items():
                if not fit['resolved'] or fit['J_MHz'] is None:
                    raise RuntimeError(f'No resolved two-mode inversion for {kind} {q}')
                ro = core['final']['readout'][q]
                levels = qubit_levels(p, q, fit['bare_fr_GHz'], ro, 8, 25)
                coarse = qubit_levels(p, q, fit['bare_fr_GHz'], ro, 6, 15)
                chi = levels['two_chi_MHz']
                modes = [np.sort_complex(np.linalg.eigvals(modal_matrix(fit, chi, s, 0))) for s in (-.5, .5)]
                shifts = (modes[1] - modes[0]).real
                nominal = np.sort_complex(np.linalg.eigvals(modal_matrix(fit, chi, 0, 0)))
                np.testing.assert_allclose(-2 * nominal.imag, fit['FWHM_MHz'], rtol=1e-5)
                ambient = purcell(q, levels, ro['beta = Cc/Csigma'], cap, g_args, g_kw)
                sensitivity = purcell(q, levels, ro['beta = Cc/Csigma'], cap, g_args, dict(g_kw, feed_res=0.))
                nmap = core['net_mapping']
                xy = Q.cm.charge_line(matrix, (nmap[q],), nmap['XY1' if q == 'Q1' else 'XY2'],
                                      levels['f01_GHz'] * 1e9, levels['n01'])
                f = levels['f01_GHz'] * 1e9
                omega = 2 * np.pi * f
                thermal_ratio = np.exp(-h * f / (k_B * TEFF_K))
                thermal = thermal_ratio / (1 + thermal_ratio)
                rates = dict(cap=omega / QCAP, qp=XQP * np.sqrt(2 * e * AL_GAP_eV * omega / hbar) / np.pi,
                             Purcell=1e6 / ambient['T1_us'], XY=1e6 / xy['T1 limit [us]'])
                down = sum(rates.values())
                up = thermal_ratio * down
                report[kind][q] = dict(levels=levels, chi_truncation_delta_MHz=chi - coarse['two_chi_MHz'],
                                       hybrid_state_separation_MHz=shifts.tolist(),
                                       separation_over_FWHM=(abs(shifts) / fit['FWHM_MHz']).tolist(),
                                       estimated_output_share=[output_share(f, q, cap, g_args[1]) for f in fit['centers_GHz']],
                                       t1_limits_us={key: 1e6 / rate for key, rate in rates.items()},
                                       total_T1_us=1e6 / (down + up), purcell_impedance=ambient,
                                       purcell_T1_us_without_direct_feed_res=sensitivity['T1_us'],
                                       thermal_population=thermal, thermal_no_reset_fidelity_ceiling=1 - thermal,
                                       pulse_model=assignment(fit, q, chi, cap, g_args[1], down / 1e6, up / 1e6, thermal))
        finish_run(folder, status='completed', summary=report)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
