"""Length <-> frequency and hammer C_j <-> coupling helpers for the Purcell + readout module.

Circuit model (not a QComponent) for the layout ``pattern design/resonator/purcell_readout.PurcellReadout``.
Quarter-wave filter and quarter-wave readout resonator (open at C_p / at the qubit pad, shorted at the meander
end), joined at one tap each by the hammer C_j. Methods: note/purcell_filter_method.md (quarter-wave frequency
and short-tail trim, closed-form J; layout and parasitic numbers in section 6), skill/purcell_resonator_design.md
(order of work, CLAW_DL). Everything here is a circuit/closed-form estimate: take final f_p, f_r and J
from the matching HFSS fixture (skill/purcell_resonator_design.md section 2).

Units: lengths in mm, capacitances in fF, frequencies in GHz, J / kappa in MHz unless the name says otherwise.

Contents
* ``CPW``: line presets (Z0, v, C').
* ``lambda4_f`` / ``lambda2_f``: open-end-loaded line         * ``tapped_lambda4_f``: exact TL root with tap stub
* ``length_for_f``: inverse (total length for a target f)    * ``hammer_cj`` / ``hammer_excess_ground``: Q3D fit
* ``coupling``: J (closed form, hammer only), kappa_p, kappa_r, plus the parallel-run parasitic J
* ``parallel_j``: distributed coupling of side-by-side runs (2-D c'_m table ``CPW_PAIR_2D`` below)
Route lengths after fillets: ``purcell_readout.filleted_length`` (layout side).
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, least_squares

# ---- CPW presets ------------------------------------------------------------------------------------------
# z0 [Ohm], v [m/s]; c_line [fF/mm] = 1 / (Z0 v). Sources:
#  10/6: skill/chip_drc.md A.1 (conformal mapping; eps_eff 6.225 = (1 + 11.45) / 2, i.e. Si at 10 mK, eps_r 11.45)
#  16/9: skill/purcell_resonator_design.md section 2 / Design4_Single_purcell-res quasi-static model (eps_r 11.5, eps_eff 6.25)
#  *_eps11.9: same cross-sections with eps_eff (1 + 11.9) / 2 = 6.45 (HFSS library 'silicon'); Z0 scaled by
#  sqrt(6.225 / 6.45) resp. sqrt(6.25 / 6.45). Pick the preset that matches the solver's substrate.
C0 = 299792458.0
CPW = {
    '10/6': dict(width_um=10, gap_um=6, z0=50.9, v=1.2016e8, source='skill/chip_drc.md A.1 (eps_r 11.45)'),
    '16/9': dict(width_um=16, gap_um=9, z0=49.89486, v=1.199169832e8,
                 source='Design4_Single_purcell-res/notes/hammer_fixed_length_validation_20261116.md (eps_r 11.5)'),
}
CPW['10/6_eps11.9'] = dict(CPW['10/6'], z0=50.9 * np.sqrt(6.225 / 6.45), v=C0 / np.sqrt(6.45),
                           source='10/6 scaled to eps_eff 6.45 (eps_r 11.9)')
CPW['16/9_eps11.9'] = dict(CPW['16/9'], z0=49.89486 * np.sqrt(6.25 / 6.45), v=C0 / np.sqrt(6.45),
                           source='16/9 scaled to eps_eff 6.45 (eps_r 11.9)')
for _c in CPW.values():
    _c['c_line'] = 1e15 / (_c['z0'] * _c['v']) / 1e3          # fF/mm


def line(cpw='10/6'):
    return CPW[cpw] if isinstance(cpw, str) else cpw


def open_end_c(cpw='10/6'):
    """Fringe C of an etched CPW open end (end gap = trace gap), fF: the usual extension
    Delta_l ~ (W + 2S) / 4 (quasi-static rule of thumb, not simulated here)."""
    c = line(cpw)
    return c['c_line'] * (c['width_um'] + 2 * c['gap_um']) / 4 / 1000


# ---- frequency from length ---------------------------------------------------------------------------------
def lambda4_f(L_mm, c_end_fF=0.0, cpw='10/6'):
    """Open-short lambda/4: f = v / (4 (L + C_end / C')), GHz (note/purcell_filter_method.md step 1)."""
    c = line(cpw)
    return c['v'] / (4 * (L_mm + c_end_fF / c['c_line']) * 1e-3) / 1e9


def lambda2_f(L_mm, c_end1_fF=0.0, c_end2_fF=0.0, cpw='10/6'):
    """Open-open lambda/2: f = v / (2 (L + (C1 + C2) / C')), GHz."""
    c = line(cpw)
    return c['v'] / (2 * (L_mm + (c_end1_fF + c_end2_fF) / c['c_line']) * 1e-3) / 1e9


def _y_open(f, length, c_load, c):
    """Input admittance (S, imaginary part) of a line of `length` mm ending in a shunt C [fF] (open end)."""
    y0 = 1 / c['z0']
    th = 2 * np.pi * f * 1e9 * length * 1e-3 / c['v']
    yl = 2 * np.pi * f * 1e9 * c_load * 1e-15
    t = np.tan(th)
    return y0 * (yl + y0 * t) / (y0 - yl * t)


def _y_short(f, length, c):
    th = 2 * np.pi * f * 1e9 * length * 1e-3 / c['v']
    return -1 / (c['z0'] * np.tan(th))


def tap_susceptance(f, L_mm, tap_mm, c_end_fF, stub_mm, c_stub_fF, cpw='10/6'):
    """Total susceptance (S) at the tap of a shorted lambda/4 with a tapped open stub.
    tap_mm: tap -> open end (along the line); L_mm - tap_mm: tap -> short. The stub (hammer arm) of stub_mm
    ends in c_stub_fF (head excess ground + C_j for the bare mode: the other line is held at ground)."""
    c = line(cpw)
    return (_y_open(f, tap_mm, c_end_fF, c) + _y_short(f, L_mm - tap_mm, c)
            + (_y_open(f, stub_mm, c_stub_fF, c) if stub_mm > 0 or c_stub_fF else 0.0))


def tapped_lambda4_f(L_mm, tap_mm, c_end_fF=0.0, stub_mm=0.0, c_stub_fF=0.0, cpw='10/6'):
    """Bare fundamental (GHz) of the tapped quarter-wave line: zero of tap_susceptance (Foster: B rises
    through 0) below the unloaded v / 4L. Method: note/purcell_filter_method.md step 2."""
    f_hi = lambda4_f(L_mm, 0, cpw) * 0.9999
    f_lo = lambda4_f(L_mm, c_end_fF + c_stub_fF, cpw) * 0.6
    fs = np.linspace(f_lo, f_hi, 600)
    b = np.array([tap_susceptance(f, L_mm, tap_mm, c_end_fF, stub_mm, c_stub_fF, cpw) for f in fs])
    roots = [i for i in range(len(fs) - 1) if b[i] < 0 < b[i + 1]]
    if not roots:
        raise ValueError('no quarter-wave root found; check lengths / loads')
    i = roots[-1]
    return brentq(lambda f: tap_susceptance(f, L_mm, tap_mm, c_end_fF, stub_mm, c_stub_fF, cpw), fs[i], fs[i + 1],
                  xtol=1e-9)


def length_for_f(f_ghz, tap_mm, c_end_fF=0.0, stub_mm=0.0, c_stub_fF=0.0, cpw='10/6'):
    """Total filleted length (mm, open end -> short, through the tap) that puts the bare mode at f_ghz.
    The tap position stays fixed; the shorted arm absorbs the change (short-tail trim, method note step 4)."""
    L0 = lambda4_f(1.0, 0, cpw) / f_ghz                               # unloaded lambda/4
    lo, hi = max(tap_mm + 0.05, 0.5 * L0), 1.05 * L0
    return brentq(lambda L: tapped_lambda4_f(L, tap_mm, c_end_fF, stub_mm, c_stub_fF, cpw) - f_ghz, lo, hi,
                  xtol=1e-9)


def trim_slope(L_mm, tap_mm, c_end_fF=0.0, stub_mm=0.0, c_stub_fF=0.0, cpw='10/6', h_mm=0.005):
    """s = df / dL_short in MHz/um (central difference), method note step 4."""
    f = lambda L: tapped_lambda4_f(L, tap_mm, c_end_fF, stub_mm, c_stub_fF, cpw)
    return (f(L_mm + h_mm) - f(L_mm - h_mm)) / (2 * h_mm) * 1e3 / 1e3


# ---- hammer C_j (Q3D data and fit) ---------------------------------------------------------------------------
# Local Q3D mutual capacitance of the hammer (heads + arms, CPW 16/9 um, eps_r 11.5, head_width 25 um,
# ground_gap 16 um, corner 2 um). 'C_m' = -Maxwell off-diagonal (fF). All rows are Q3D results; convergence
# as stated in the source.
HAMMER_Q3D = [
    dict(gap=40, height=30, width=25, C_m=0.49278, src='Design4_Single_purcell-res run 20261116_195508_h30_g40_w25_a55_q3d_cae2f56a (20 passes, converged)'),
    dict(gap=40, height=45, width=25, C_m=0.84207, src='Design4_Single_purcell-res run 20261116_200033_h45_g40_w25_a55_q3d_73df346f (converged)'),
    dict(gap=40, height=60, width=25, C_m=1.25486, src='Design4_Single_purcell-res run 20261116_200156_h60_g40_w25_a55_q3d_33be5e01 (converged)'),
    dict(gap=40, height=80, width=25, C_m=1.84797, src='Design4_Single_purcell-res run 20261116_200519_h80_g40_w25_a55_q3d_58545fab (converged)'),
    dict(gap=42, height=65, width=25, C_m=1.45200, src='Design4_Single_purcell-res hammer_cpw_refine_2d2775df, 800 um leads (extension not converged)'),
    dict(gap=20, height=65, width=25, C_m=2.82917, src='Design4_Single_purcell-res hammer_fixed_length_validation note, 800 um leads'),
    dict(gap=42, height=100, width=25, C_m=2.52294, src='Design4_Single_purcell-res hammer_fixed_length_validation note, 800 um leads'),
]
# Deembedded excess ground per head (fF) relative to a straight CPW of the same outer length (same sources).
HAMMER_EXCESS_Q3D = [
    dict(gap=42, height=65, C_ex=2.40), dict(gap=20, height=65, C_ex=3.36), dict(gap=42, height=100, C_ex=6.06),
    dict(gap=20, height=30, C_ex=-0.356),
]
# Validity of the fit: width 25 um, gap 20-42 um, height 30-100 um, 16/9 um handles. Design5_6_purcell-res r3 (width 50, gap 50,
# height 65: Q3D 1.476 fF) is 15 % above the width-25 fit: wider heads couple more, refit before using them.


@lru_cache(maxsize=None)
def _cj_fit():
    d = np.array([[r['gap'], r['height'], r['C_m']] for r in HAMMER_Q3D])
    res = least_squares(lambda q: (d[:, 1] - q[2]) * (q[0] / d[:, 0] + q[1]) / d[:, 2] - 1, [1.0, 0.004, 12.0])
    a, b, h0 = res.x
    rel = (d[:, 1] - h0) * (a / d[:, 0] + b) / d[:, 2] - 1
    return dict(a_fF=a, b_fF_per_um=b, h0_um=h0, max_rel_err=float(abs(rel).max()), rel_err=rel.tolist())


def hammer_cj(gap_um, height_um, width_um=25.0):
    """Hammer mutual C_j (fF) = (height - h0) (a / gap + b), least-squares fit to HAMMER_Q3D (max error printed by
    hammer_fit_report()). Estimate only: run the Q3D coupon (purcell_readout.hammer_coupon) for a final value."""
    if not 15 <= gap_um <= 50 or not 25 <= height_um <= 110:
        raise ValueError('hammer_cj: outside the fitted range (gap 20-42, height 30-100 um)')
    if abs(width_um - 25) > 5:
        raise ValueError('hammer_cj: fit only for head_width ~25 um; run Q3D for other widths')
    q = _cj_fit()
    return (height_um - q['h0_um']) * (q['a_fF'] / gap_um + q['b_fF_per_um'])


def gap_for_cj(cj_fF, height_um, width_um=25.0):
    """Inverse of hammer_cj at fixed height: head gap (um) for a target C_j."""
    return brentq(lambda g: hammer_cj(g, height_um, width_um) - cj_fF, 15, 50)


def height_for_cj(cj_fF, gap_um, width_um=25.0):
    """Inverse of hammer_cj at fixed gap: head height (um) for a target C_j."""
    return brentq(lambda h: hammer_cj(gap_um, h, width_um) - cj_fF, 25, 110)


@lru_cache(maxsize=None)
def _excess_fit():
    d = np.array([[r['gap'], r['height'], r['C_ex']] for r in HAMMER_EXCESS_Q3D])
    A = np.c_[np.ones(len(d)), d[:, 1], d[:, 1] / d[:, 0]]
    coef, *_ = np.linalg.lstsq(A, d[:, 2], rcond=None)
    return coef


def hammer_excess_ground(gap_um, height_um):
    """Excess ground C of one head over a straight CPW (fF), linear fit c0 + c1 h + c2 h / g to
    HAMMER_EXCESS_Q3D (4 points, about +-1 fF). Enters the frequency only through the tap (~1 MHz per fF)."""
    c0, c1, c2 = _excess_fit()
    return float(c0 + c1 * height_um + c2 * height_um / gap_um)


def hammer_fit_report():
    q = _cj_fit()
    rows = [dict(r, fit=hammer_cj(r['gap'], r['height'])) for r in HAMMER_Q3D]
    return dict(q, rows=rows, excess_coef=_excess_fit().tolist())


# ---- coupling -------------------------------------------------------------------------------------------------
def j_hammer(f_ghz, cj_fF, L_p, tap_p, c_end_p, L_r, tap_r, c_end_r, cpw='10/6'):
    """Closed-form J (MHz) of the hammer alone, code/general/purcell_j.py formula:
    J = (f/2) C_j cos(theta_p) cos(theta_r) / sqrt(C_eff,p C_eff,r), open-end loads as line extensions."""
    cl = line(cpw)['c_line']
    xp, xr = c_end_p / cl, c_end_r / cl
    cp = np.cos(np.pi / 2 * (tap_p + xp) / (L_p + xp))
    cr = np.cos(np.pi / 2 * (tap_r + xr) / (L_r + xr))
    ce_p, ce_r = cl * (L_p + xp) / 2, cl * (L_r + xr) / 2
    return float(f_ghz * 1e3 / 2 * cj_fF * cp * cr / np.sqrt(ce_p * ce_r))


def kappa_hanger(f_ghz, c_p_fF, cpw='10/6', z_env=None):
    """kappa_p / 2pi (MHz) of a lambda/4 filter whose open end couples through C_p to a matched feed line
    (environment Z0/2): kappa = 2 w^3 C_p^2 Z_env Zr ... with Z_env = Z0 / 2 -> kappa = 4 w^3 C_p^2 Z_env Zr / pi.
    First estimate only; HFSS / Q3D values for C_p in skill/purcell_resonator_design.md section 1."""
    c = line(cpw)
    w = 2 * np.pi * f_ghz * 1e9
    ze = c['z0'] / 2 if z_env is None else z_env
    return float(4 * w ** 3 * (c_p_fF * 1e-15) ** 2 * ze * c['z0'] / np.pi / (2 * np.pi) / 1e6)


def kappa_readout(j_mhz, kappa_p_mhz, delta_mhz=0.0):
    """Effective readout linewidth (MHz) through the filter, weak-J limit: 4 J^2 kappa_p / (kappa_p^2 + 4 Delta^2)."""
    return float(4 * j_mhz ** 2 * kappa_p_mhz / (kappa_p_mhz ** 2 + 4 * delta_mhz ** 2))


# ---- parallel-run parasitic coupling -----------------------------------------------------------------------------
# 2-D cross-section extraction of two parallel CPWs (code/resonator/cpw_pair_cap2d.py, eps_r 11.9; run
# 2026-11-18T23:12:15, script sha256 3dbe8160... = git 5964ee0, Python 3.12.13 / numpy 1.26.4 / scipy 1.17.1).
# Rows: (centre spacing um, c'_m fF/mm, k_C, k_L). Grid check: c'_m changes < 1.3 % from h0 1 um to 0.25 um;
# 16/9 at 200 um, eps_r 11.5 gives 0.546 fF/mm against 0.563 fF/mm from the Design4_Single_purcell-res Q3D pair coupons.
# Regenerate with `python -I -B code/resonator/cpw_pair_cap2d.py --out <run dir>` and paste the rows here.
CPW_PAIR_2D = {
    '10_6_eps11.9': [
        (150, 0.40008, 0.0023224, 0.0022779),
        (200, 0.22874, 0.0013278, 0.0012865),
        (250, 0.14922, 0.0008662, 0.00082871),
        (300, 0.10573, 0.00061378, 0.00058064),
        (400, 0.061525, 0.00035714, 0.00033348),
        (500, 0.04004, 0.00023243, 0.00021822),
        (600, 0.027698, 0.00016079, 0.00015507),
        (800, 0.014516, 8.4266e-05, 9.0948e-05),
    ],
    '16_9_eps11.9': [
        (150, 0.99001, 0.0056709, 0.0055624),
        (200, 0.56354, 0.003228, 0.0031275),
        (250, 0.36733, 0.0021041, 0.0020131),
        (300, 0.25968, 0.0014874, 0.001407),
        (400, 0.15093, 0.00086451, 0.00080705),
        (500, 0.098215, 0.00056258, 0.00052803),
        (600, 0.067938, 0.00038915, 0.00037517),
        (800, 0.035603, 0.00020393, 0.00021999),
    ],
}


def cm_per_length(spacing_um, cpw='10/6'):
    """2-D c'_m (fF/mm), k_C and k_L of two parallel CPWs at a centre spacing (log-interpolated in
    CPW_PAIR_2D, eps_r 11.9). Returns dict(c_m, k_C, k_L)."""
    key = {'10/6': '10_6_eps11.9', '16/9': '16_9_eps11.9'}[cpw.split('_')[0] if isinstance(cpw, str) else '10/6']
    t = np.array(CPW_PAIR_2D[key], float)
    return {k: float(np.exp(np.interp(np.log(spacing_um), np.log(t[:, 0]), np.log(t[:, i]))))
            for i, k in ((1, 'c_m'), (2, 'k_C'), (3, 'k_L'))}


def parallel_j(f_ghz, runs, L_p, c_end_p, L_r, c_end_r, cpw='10/6'):
    """Distributed coupling J (MHz) of side-by-side runs between the filter and the resonator.

    runs: list of dict(spacing_um, zp=(a, b), zr=(a, b)) with the run's positions along the filter and the
    resonator, measured from each line's open end (mm); zr may run backwards (b < a) when the two lines run in
    opposite directions. Quarter-wave profiles V ~ cos(pi z / 2 L_eff), I ~ sin(...), open-end loads as line
    extensions. Returns dict(J_C, J_L, J_max = |J_C| + |J_L|, J_min = ||J_C| - |J_L||): the sign between the
    capacitive and inductive parts depends on the relative current direction and is not resolved here, so use
    J_max as the bound (quasi-TEM, k_L ~ k_C: they cancel for co-directional runs near the short ends).
    """
    cl = line(cpw)['c_line']
    xp, xr = c_end_p / cl, c_end_r / cl
    Lp, Lr = L_p + xp, L_r + xr
    jc = jl = 0.0
    for run in runs:
        k = cm_per_length(run['spacing_um'], cpw)
        n = 200
        t = (np.arange(n) + 0.5) / n
        zp = run['zp'][0] + (run['zp'][1] - run['zp'][0]) * t + xp
        zr = run['zr'][0] + (run['zr'][1] - run['zr'][0]) * t + xr
        dz = abs(run['zp'][1] - run['zp'][0]) / n
        vp, vr = np.cos(np.pi / 2 * zp / Lp), np.cos(np.pi / 2 * zr / Lr)
        ip, ir = np.sin(np.pi / 2 * zp / Lp), np.sin(np.pi / 2 * zr / Lr)
        norm = 2 / np.sqrt(Lp * Lr)          # mode normalisation: int_0^L cos^2 = L / 2
        jc += f_ghz * 1e3 / 2 * k['k_C'] * norm * np.sum(vp * vr) * dz
        jl += f_ghz * 1e3 / 2 * k['k_L'] * norm * np.sum(ip * ir) * dz
    return dict(J_C=float(jc), J_L=float(jl), J_max=float(abs(jc) + abs(jl)), J_min=float(abs(abs(jc) - abs(jl))))


def coupling(f_ghz, cj_fF, L_p, tap_p, c_end_p, L_r, tap_r, c_end_r, c_p_fF=None, delta_mhz=0.0,
             runs=(), cpw='10/6'):
    """J from the hammer, the parallel-run parasitic bound and (with c_p_fF) kappa_p and kappa_r. MHz."""
    out = dict(J_hammer=j_hammer(f_ghz, cj_fF, L_p, tap_p, c_end_p, L_r, tap_r, c_end_r, cpw))
    par = parallel_j(f_ghz, list(runs), L_p, c_end_p, L_r, c_end_r, cpw) if runs else dict(J_max=0.0)
    out['J_parallel_max'] = par['J_max']
    out['parasitic_ratio'] = par['J_max'] / out['J_hammer'] if out['J_hammer'] else np.inf
    if c_p_fF is not None:
        out['kappa_p'] = kappa_hanger(f_ghz, c_p_fF, cpw)
        out['kappa_r'] = kappa_readout(out['J_hammer'], out['kappa_p'], delta_mhz)
    return out
