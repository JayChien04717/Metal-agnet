"""Closed-form filter-resonator coupling J of two shorted lambda/4 lines joined by a lumped capacitor C_j at the taps.

J = (f / 2) * C_j * cos(theta_p) * cos(theta_r) / sqrt(C_eff,p * C_eff,r)

* Each open-end load C_end is replaced by an equivalent line extension x = C_end / C' (C' = 1 / (Z0 v)).
* theta = (pi / 2) * (tap + x) / (L + x), with tap the distance from the physical open end and L the physical length.
* C_eff = C' (L + x) / 2 (lambda/4 mode capacitance).

Assumptions: lossless TEM lines, C_j small (no frequency pull), no other coupling path (parallel lines, ground return).
Accuracy: within ~5 % of the full Wallraff circuit for thesis Table 2.4 (see note/purcell_filter_method.md section 2).
"""
import numpy as np


def j_closed_form(f_mhz, cj_fF, c_line_fF_per_mm, tap_p_mm, L_p_mm, c_end_p_fF, tap_r_mm, L_r_mm, c_end_r_fF):
    """Returns dict(J_MHz, cos_p, cos_r, ext_p_mm, ext_r_mm, c_eff_p_fF, c_eff_r_fF). f in MHz -> J in MHz."""
    xp, xr = c_end_p_fF / c_line_fF_per_mm, c_end_r_fF / c_line_fF_per_mm
    tp, tr = L_p_mm + xp, L_r_mm + xr
    cp, cr = np.cos(np.pi / 2 * (tap_p_mm + xp) / tp), np.cos(np.pi / 2 * (tap_r_mm + xr) / tr)
    ce_p, ce_r = c_line_fF_per_mm * tp / 2, c_line_fF_per_mm * tr / 2
    return dict(J_MHz=float(f_mhz / 2 * cj_fF * cp * cr / np.sqrt(ce_p * ce_r)), cos_p=float(cp), cos_r=float(cr),
                ext_p_mm=float(xp), ext_r_mm=float(xr), c_eff_p_fF=float(ce_p), c_eff_r_fF=float(ce_r))
