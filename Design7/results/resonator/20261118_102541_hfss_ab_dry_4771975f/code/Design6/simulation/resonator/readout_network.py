"""Numerical Design6 variant G admittance with optional series Cin and readout probe."""
import numpy as np

NAMES = ("Q1", "Q2")
FILTER_TAP, RES_TAP, TEE_HALF = 1.20, 1.30, 0.15  # mm, Design4/layout/purcell_6500.py


def port_Y(f, lengths, geo, cp, hammer, sec, cm=0.0, feed_res=0.0, extra_cj=0.0, coupled=True,
           *, cin_fF=0., cin_arm_mm=None, readout_probe=None):
    """Feed port admittance, with an optional third port at a readout voltage antinode."""
    w = 2 * np.pi * np.asarray(f, complex) * 1e9
    z0, v = geo['z0_ohm'], geo['velocity_m_s']
    cself = 1 / (z0 * v); cmm = cm * 1e-12             # F/m
    Yc = v * np.array([[cself + cmm, -cmm], [-cmm, cself + cmm]])
    per = 14; core_n = 4 + per * 2
    n = core_n + 2 * bool(cin_fF) + bool(readout_probe)
    Y = np.zeros((n, n), complex)

    def tl(a, b, mm):
        t = w / v * mm * 1e-3; dg, off = -1j / z0 / np.tan(t), 1j / z0 / np.sin(t)
        Y[a, a] += dg; Y[b, b] += dg; Y[a, b] += off; Y[b, a] += off

    def ctl(a, b, mm):
        t = w / v * mm * 1e-3; dg, off = -1j * Yc / np.tan(t), 1j * Yc / np.sin(t)
        Y[np.ix_(a, a)] += dg; Y[np.ix_(b, b)] += dg; Y[np.ix_(a, b)] += off; Y[np.ix_(b, a)] += off

    def pi(a, b, c):
        cc, ga, gb = np.asarray(c) * 1e-15
        Y[a, a] += 1j * w * (cc + ga); Y[b, b] += 1j * w * (cc + gb); Y[a, b] -= 1j * w * cc; Y[b, a] -= 1j * w * cc

    def stub(a, mm, shorted):
        t = w / v * mm * 1e-3
        Y[a, a] += -1j / z0 / np.tan(t) if shorted else 1j / z0 * np.tan(t)

    F = geo['feeds_mm']
    tl(0, 2, F['feed_in']); tl(2, 3, F['feed_mid'])
    if cin_fF:
        if not 0 < cin_arm_mm < F['feed_out']:
            raise ValueError('Input capacitor must be between the Q2 tee and P7')
        tl(3, core_n, cin_arm_mm)
        pi(core_n, core_n + 1, (cin_fF, 0, 0))
        tl(core_n + 1, 1, F['feed_out'] - cin_arm_mm)
    else:
        tl(3, 1, F['feed_out'])
    for i, q in enumerate(NAMES):
        a, b, f1, ft, tf, fb, f0, r1, rt, tr, rb, r0, hl, hr = range(4 + per * i, 4 + per * (i + 1))
        Lf, Lr = lengths[q]; s = sec[q]; th = TEE_HALF
        tl(2 + i, a, F[f'stub_p_{q}']); pi(a, b, cp)
        if coupled:
            tl(b, f1, s['filter_free_mm'])
            ctl([f1, r1], [ft, rt], s['above_mm']); ctl([ft, rt], [tf, tr], th)
            ctl([tf, tr], [fb, rb], th); ctl([fb, rb], [f0, r0], s['below_mm'])
            stub(f0, Lf - FILTER_TAP - th - s['below_mm'], True)
            if readout_probe == q:
                tl(r0, n - 1, s['res_open_free_mm'])
            else:
                stub(r0, s['res_open_free_mm'], False)
            stub(r1, Lr - RES_TAP - th - s['above_mm'], True)
            if feed_res:
                pi(a, r1, (feed_res, 0, 0))
        else:
            tl(b, tf, FILTER_TAP); stub(tf, Lf - FILTER_TAP, True)
            stub(tr, RES_TAP, False); stub(tr, Lr - RES_TAP, True)
            for x in (f1, ft, fb, f0, r1, rt, rb, r0):
                Y[x, x] += 1.0                                       # unused nodes, isolated
        sp = geo['spacing_mm'][q] / 2
        tl(tf, hl, sp); pi(hl, hr, (hammer[0] + extra_cj, hammer[1], hammer[2])); tl(hr, tr, sp)
    if readout_probe and (readout_probe not in NAMES or not coupled):
        raise ValueError('Readout probe requires a named channel in the coupled model')
    ports = [0, 1] + ([n - 1] if readout_probe else [])
    internal = list(range(2, n - bool(readout_probe)))
    return (Y[np.ix_(ports, ports)] - Y[np.ix_(ports, internal)]
            @ np.linalg.solve(Y[np.ix_(internal, internal)], Y[np.ix_(internal, ports)]))
