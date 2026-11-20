"""Small-JJ EJ <-> normal-state resistance <-> overlap area <-> e-beam mask widths.

Rules and tables: skill/ebeam_junction.md. Every quantity is per single JJ,
including each branch of a SQUID. EJ is a design input, never an HFSS result.

    EJ/h [GHz]  --Ambegaokar-Baratoff (Al, Delta = 180 ueV, T -> 0)-->  Rn [kOhm]
                   Rn [kOhm] = 140.433954 / (EJ/h [GHz])
    Rn          --RA = 625 Ohm um^2-->  actual overlap area A = RA / Rn
    A           = (jj_finger_width + 30 nm) * (jj_cross_electrode_width + 70 nm)

* ``jj_finger_width``: mask x width of the vertical finger (default 80 nm, fixed).
* ``jj_cross_electrode_width``: mask y width of the horizontal electrode; this is
  the parameter that sets the area. Rounded to the 1 nm GDS grid.
* "actual" values are predicted fabricated sizes; never use them as mask sizes.

Calibration status: RA comes from one measured junction (0.100 x 0.250 um actual,
about 25 kOhm); the +30 / +70 nm growth values are process assumptions. Results
are estimates. After fabrication, fit them to test-key resistances and pass the
new values (``ra``, ``finger_growth_nm``, ``cross_growth_nm``, ``delta_uev``).

Usage (standard library only)::

    python -I code/qubit/jj_resistance_model.py --ej 4 8.435                 # widths for EJ
    python -I code/qubit/jj_resistance_model.py --sweep                      # width -> area -> Rn -> EJ
    python -I code/qubit/jj_resistance_model.py --squid 23.056 0.4227        # EJ_sum, asymmetry d
    python -I code/qubit/jj_resistance_model.py --squid-flux 23.056 9.746    # EJ(0), EJ(Phi0/2)
    python -I code/qubit/jj_resistance_model.py --check                      # self-consistency
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from math import isfinite

E_CHARGE = 1.602176634e-19          # C (exact, SI 2019)
DELTA_UEV = 180.0                   # Al gap, low-temperature limit
RA_OHM_UM2 = 625.0                  # 25 kOhm * 0.100 um * 0.250 um (one measured calibration JJ)
FINGER_WIDTH_NM = 80.0              # mask finger width (fixed for every small JJ)
FINGER_GROWTH_NM = 30.0             # 80 nm mask -> 110 nm actual (process assumption)
CROSS_GROWTH_NM = 70.0              # horizontal fabrication offset (process assumption)
GRID_NM = 1.0                       # GDS grid used for the applied width


def rn_ej_product(delta_uev: float = DELTA_UEV) -> float:
    """Rn [kOhm] * EJ/h [GHz] from Ic Rn = pi Delta / 2e and EJ = hbar Ic / 2e.

    Rn = Delta / (8 e^2 f) -> 140.433954 kOhm GHz at Delta = 180 ueV.
    """
    return delta_uev * 1e-6 / (8 * E_CHARGE * 1e9) / 1e3


def rn_kohm(ej_ghz: float, delta_uev: float = DELTA_UEV) -> float:
    """Normal-state resistance [kOhm] of one JJ with EJ/h = ``ej_ghz``."""
    _positive(ej_ghz=ej_ghz)
    return rn_ej_product(delta_uev) / ej_ghz


def ej_ghz(rn: float, delta_uev: float = DELTA_UEV) -> float:
    """EJ/h [GHz] of one JJ with normal-state resistance ``rn`` [kOhm]."""
    _positive(rn=rn)
    return rn_ej_product(delta_uev) / rn


def area_um2(rn: float, ra: float = RA_OHM_UM2) -> float:
    """Actual (fabricated) overlap area [um^2] for ``rn`` [kOhm]."""
    _positive(rn=rn)
    return ra / (rn * 1e3)


def rn_from_area(area: float, ra: float = RA_OHM_UM2) -> float:
    """Normal-state resistance [kOhm] of an actual overlap ``area`` [um^2]."""
    _positive(area=area)
    return ra / area / 1e3


@dataclass(frozen=True)
class JunctionSize:
    """One small JJ: design inputs, estimate chain and mask widths (nm)."""
    ej_ghz: float
    rn_kohm: float
    actual_area_um2: float
    jj_finger_width_nm: float              # mask vertical finger
    actual_finger_width_nm: float          # predicted fabricated
    calculated_cross_electrode_width_nm: float   # unrounded mask width
    actual_cross_electrode_width_nm: float       # predicted fabricated
    jj_cross_electrode_width_nm: float     # applied mask width, rounded to grid

    def options(self) -> dict:
        """Component options for FluxoniumJunction / TransmonJunction (one JJ)."""
        return dict(jj_finger_width=f'{self.jj_finger_width_nm:g}nm',
                    jj_cross_electrode_width=f'{self.jj_cross_electrode_width_nm:g}nm')


def size_from_ej(ej: float, finger_width_nm: float = FINGER_WIDTH_NM, *,
                 delta_uev: float = DELTA_UEV, ra: float = RA_OHM_UM2,
                 finger_growth_nm: float = FINGER_GROWTH_NM,
                 cross_growth_nm: float = CROSS_GROWTH_NM,
                 grid_nm: float = GRID_NM) -> JunctionSize:
    """EJ/h [GHz] -> Rn -> area -> mask widths, with ``jj_finger_width`` fixed.

    Uses the unrounded Rn and rounds only the final
    mask width to ``grid_nm``.
    """
    _positive(finger_width_nm=finger_width_nm, grid_nm=grid_nm)
    rn = rn_kohm(ej, delta_uev)
    area = area_um2(rn, ra)
    actual_finger = finger_width_nm + finger_growth_nm
    actual_cross = area * 1e6 / actual_finger
    calculated = actual_cross - cross_growth_nm
    if calculated <= 0:
        raise ValueError(f'EJ/h = {ej} GHz needs a non-positive cross electrode '
                         f'({calculated:.1f} nm); use a narrower jj_finger_width')
    applied = round(calculated / grid_nm) * grid_nm
    return JunctionSize(ej, rn, area, finger_width_nm, actual_finger,
                        calculated, actual_cross, applied)


def finger_from_ej(ej: float, cross_electrode_width_nm: float, **kw) -> float:
    """Unrounded mask finger width [nm] for EJ with a fixed cross electrode."""
    rn = rn_kohm(ej, kw.get('delta_uev', DELTA_UEV))
    area = area_um2(rn, kw.get('ra', RA_OHM_UM2))
    actual_cross = cross_electrode_width_nm + kw.get('cross_growth_nm', CROSS_GROWTH_NM)
    return area * 1e6 / actual_cross - kw.get('finger_growth_nm', FINGER_GROWTH_NM)


def ej_from_size(finger_width_nm: float, cross_electrode_width_nm: float, *,
                 delta_uev: float = DELTA_UEV, ra: float = RA_OHM_UM2,
                 finger_growth_nm: float = FINGER_GROWTH_NM,
                 cross_growth_nm: float = CROSS_GROWTH_NM) -> dict:
    """Mask widths [nm] -> predicted actual area, Rn [kOhm] and EJ/h [GHz]."""
    area = ((finger_width_nm + finger_growth_nm)
            * (cross_electrode_width_nm + cross_growth_nm) * 1e-6)
    rn = rn_from_area(area, ra)
    return dict(actual_area_um2=area, rn_kohm=rn, ej_ghz=ej_ghz(rn, delta_uev))


def squid_options(ej1: float, ej2: float | None = None, **kw) -> dict:
    """TransmonJunction options for a SQUID from the per-JJ EJ1 / EJ2 [GHz].

    Use the individual-branch energies; a symmetric SQUID with a zero-flux
    total EJ_sum has EJ1 = EJ2 = EJ_sum / 2 (loop inductance neglected).
    """
    a = size_from_ej(ej1, **kw)
    b = size_from_ej(ej1 if ej2 is None else ej2, **kw)
    if a.jj_finger_width_nm != b.jj_finger_width_nm:
        raise ValueError('both SQUID branches share jj_finger_width')
    return dict(jj_finger_width=f'{a.jj_finger_width_nm:g}nm',
                jj1_cross_electrode_width=f'{a.jj_cross_electrode_width_nm:g}nm',
                jj2_cross_electrode_width=f'{b.jj_cross_electrode_width_nm:g}nm')


def squid_branches(ej_sum: float, asymmetry: float) -> tuple[float, float]:
    """Branch energies (EJ1 >= EJ2, GHz) from EJ_sum = EJ1 + EJ2 and d = (EJ1 - EJ2) / EJ_sum.

    Ideal SQUID (loop inductance neglected):
    EJ(Phi) = EJ_sum * sqrt(cos^2(pi Phi/Phi0) + d^2 sin^2(pi Phi/Phi0)),
    so EJ(0) = EJ_sum and EJ(Phi0/2) = d * EJ_sum.
    """
    _positive(ej_sum=ej_sum)
    if not 0 <= asymmetry < 1:
        raise ValueError('asymmetry d must satisfy 0 <= d < 1')
    return ej_sum * (1 + asymmetry) / 2, ej_sum * (1 - asymmetry) / 2


def asymmetry_from_flux_points(ej_zero_flux: float, ej_half_flux: float) -> float:
    """d = EJ(Phi0/2) / EJ(0) for an ideal SQUID."""
    _positive(ej_zero_flux=ej_zero_flux)
    return ej_half_flux / ej_zero_flux


def squid_sizes(ej_sum: float, asymmetry: float, **kw) -> dict:
    """Both branch sizes and TransmonJunction options from EJ_sum and asymmetry d.

    jj_finger_width stays fixed; only the cross electrode
    width differs between the branches. Areas scale
    with EJ: A_actual = EJ * RA / (Rn * EJ product).
    """
    ej1, ej2 = squid_branches(ej_sum, asymmetry)
    a, b = size_from_ej(ej1, **kw), size_from_ej(ej2, **kw)
    return dict(jj1=a, jj2=b, area_ratio=a.actual_area_um2 / b.actual_area_um2,
                options=squid_options(ej1, ej2, **kw))


def self_check(tol=1e-9):
    """Internal consistency of the model; raises AssertionError on failure.

    * Rn * EJ = 140.433954 kOhm GHz at Delta = 180 ueV.
    * Calibration junction: 0.100 x 0.250 um actual -> 25 kOhm at RA = 625 Ohm um^2.
    * EJ -> mask widths -> EJ round trip (before grid rounding).
    * SQUID branches sum to EJ_sum and reproduce the asymmetry.
    """
    assert abs(rn_ej_product() - 140.433954) < 1e-6, rn_ej_product()
    assert abs(rn_from_area(0.100 * 0.250) - 25.0) < tol
    for ej in (2.0, 4.0, 8.0, 12.0, 16.0):
        s = size_from_ej(ej)
        back = ej_from_size(s.jj_finger_width_nm, s.calculated_cross_electrode_width_nm)
        assert abs(back['ej_ghz'] - ej) / ej < tol, (ej, back)
        assert abs(back['actual_area_um2'] - s.actual_area_um2) < tol
        assert abs(finger_from_ej(ej, s.calculated_cross_electrode_width_nm) - s.jj_finger_width_nm) < 1e-6
    e1, e2 = squid_branches(20.0, 0.4)
    assert abs(e1 + e2 - 20.0) < tol and abs((e1 - e2) / 20.0 - 0.4) < tol
    assert abs(asymmetry_from_flux_points(20.0, 8.0) - 0.4) < tol
    return True


def markdown_sweep(widths=(80, 100, 120, 150, 200, 250, 300, 350, 400, 450, 500, 600),
                   finger_width_nm=FINGER_WIDTH_NM):
    lines = [f'| `jj_cross_electrode_width` (nm), finger {finger_width_nm:g} nm | '
             'Actual area (um^2) | Rn (kOhm) | EJ/h (GHz) |', '|---:|---:|---:|---:|']
    for w in widths:
        r = ej_from_size(finger_width_nm, w)
        lines.append(f'| {w} | {r["actual_area_um2"]:.5f} | {r["rn_kohm"]:.3f} | '
                     f'{r["ej_ghz"]:.3f} |')
    return '\n'.join(lines)


def _positive(**values):
    for name, v in values.items():
        if not (isfinite(v) and v > 0):
            raise ValueError(f'{name} must be a positive finite number, got {v}')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--ej', type=float, nargs='+', help='EJ/h per JJ in GHz')
    parser.add_argument('--finger', type=float, default=FINGER_WIDTH_NM,
                        help='jj_finger_width in nm (default 80)')
    parser.add_argument('--sweep', action='store_true', help='area / Rn / EJ vs width')
    parser.add_argument('--check', action='store_true', help='internal consistency check')
    parser.add_argument('--squid', type=float, nargs=2, metavar=('EJ_SUM', 'D'),
                        help='SQUID from EJ_sum = EJ1 + EJ2 (GHz) and asymmetry d')
    parser.add_argument('--squid-flux', type=float, nargs=2, metavar=('EJ_0', 'EJ_HALF'),
                        help='SQUID from EJ at zero and half flux (GHz)')
    args = parser.parse_args()
    if args.check:
        self_check()
        print('model self-check passed')
    if args.sweep:
        print(markdown_sweep(finger_width_nm=args.finger))
    squid = args.squid
    if args.squid_flux:
        squid = (args.squid_flux[0], asymmetry_from_flux_points(*args.squid_flux))
    if squid:
        r = squid_sizes(squid[0], squid[1], finger_width_nm=args.finger)
        print(f'EJ_sum {squid[0]:g} GHz, d = {squid[1]:.6f}, area ratio JJ1/JJ2 = {r["area_ratio"]:.4f}')
        for tag in ('jj1', 'jj2'):
            s = r[tag]
            print(f'  {tag}: EJ/h {s.ej_ghz:.4f} GHz, Rn {s.rn_kohm:.3f} kOhm, actual area '
                  f'{s.actual_area_um2:.5f} um^2, jj_finger_width {s.jj_finger_width_nm:g} nm, '
                  f'{tag}_cross_electrode_width {s.jj_cross_electrode_width_nm:g} nm '
                  f'({s.calculated_cross_electrode_width_nm:.2f} unrounded)')
        print(' ', r['options'])
    for ej in args.ej or ():
        s = size_from_ej(ej, args.finger)
        print({k: round(v, 6) for k, v in asdict(s).items()}, s.options())


if __name__ == '__main__':
    main()
