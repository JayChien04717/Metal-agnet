# -*- coding: utf-8 -*-
r"""Filleted CPW path-length measurement, segment checks, and target-length solver.

Shared across Designs per ``skill/layout_structure.md`` §3 and ``skill/chip_drc.md``
rules B2 (segment >= 2 * R) and B8 (filleted centre-line length == target).
"""
import math
import re

import numpy as np


def filleted_length(coords, fillet_mm):
    """Filleted centre-line length [mm] of a polyline with bend radius ``fillet_mm``.

    Each interior corner of turning angle ``t`` is shortened by
    ``2 * r * tan(t / 2) - r * t``.
    """
    pts = np.asarray(coords, dtype=float)
    if len(pts) < 2:
        return 0.0
    segs = np.diff(pts, axis=0)
    lengths = np.linalg.norm(segs, axis=1)
    # Drop zero-length duplicate consecutive vertices before angle calculation.
    keep = np.concatenate(([True], lengths > 1e-12))
    pts = pts[keep]
    if len(pts) < 2:
        return 0.0
    total = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
    if not fillet_mm or fillet_mm <= 0.0:
        return total
    r = float(fillet_mm)
    for p0, p1, p2 in zip(pts[:-2], pts[1:-1], pts[2:]):
        a, b = p1 - p0, p2 - p1
        na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
        if na <= 1e-12 or nb <= 1e-12:
            continue
        cos_t = float(np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0))
        t = math.acos(cos_t)
        if t > 1e-9:
            total -= 2.0 * r * math.tan(t / 2.0) - r * t
    return total


def _route_row(design, comp, name=None):
    """Return the non-subtract path QGeometry row for ``comp`` (and optional ``name``)."""
    t = design.qgeometry.tables['path']
    rows = t[(t['component'] == design.components[comp].id) & (~t['subtract'])]
    if name is not None:
        rows = rows[rows['name'] == name]
    if len(rows) == 0:
        raise KeyError(f'No centre-line path row found for component {comp!r} (name={name!r})')
    return rows.iloc[0]


def path_length(design, comp, name=None):
    """Filleted centre-line length [mm] of a route or tee line in ``design``."""
    row = _route_row(design, comp, name=name)
    fil = row['fillet']
    r = float(fil) if fil is not None and not math.isnan(float(fil)) and float(fil) > 0 else 0.0
    return filleted_length(row.geometry.coords, r)


def check_min_segments(design, comp, name=None, min_factor=2.0, raise_on_violation=True):
    """Check that every interior straight segment of ``comp`` is ``>= min_factor * fillet``.

    Returns ``(min_segment_mm, required_mm)``. End segments (adjacent to only one
    bend) require ``>= fillet``; interior segments (between two bends) require
    ``>= min_factor * fillet`` (DRC B2: ``min_factor = 2.0``).
    """
    row = _route_row(design, comp, name=name)
    fil = row['fillet']
    r = float(fil) if fil is not None and not math.isnan(float(fil)) and float(fil) > 0 else 0.0
    pts = np.asarray(row.geometry.coords, dtype=float)
    segs = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    segs = segs[segs > 1e-12]
    if len(segs) == 0:
        return 0.0, min_factor * r
    min_seg = float(segs.min())
    required = float(min_factor * r)
    if r > 0.0 and len(segs) > 2:
        inner_min = float(segs[1:-1].min())
        if raise_on_violation and inner_min + 1e-9 < required:
            raise ValueError(
                f'{comp}: interior straight segment {inner_min * 1e3:.2f} um '
                f'< {min_factor:g} * fillet ({required * 1e3:.2f} um)'
            )
    return min_seg, required


def _parse_mm(value):
    """Parse a Qiskit Metal length string ('4.25mm', '4250um') or float (mm) to float mm."""
    if isinstance(value, (int, float)):
        return float(value)
    s = str(value).strip()
    m = re.fullmatch(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*(mm|um|µm|m)?', s)
    if not m:
        raise ValueError(f'Cannot parse length {value!r} as mm')
    num = float(m.group(1))
    unit = m.group(2) or 'mm'
    if unit == 'mm':
        return num
    if unit in ('um', 'µm'):
        return num * 1e-3
    if unit == 'm':
        return num * 1e3
    raise ValueError(f'Unsupported length unit {unit!r}')


def check_meander_spacing(design, comp, min_trench_um=300.0, raise_on_violation=True):
    """Check DRC B7: adjacent parallel meander sections >= ``min_trench_um`` trench-to-trench.

    For a route with ``options.meander.spacing``, measures the ground strip width
    ``spacing_um - cut_width_um`` (where ``cut_width_um = W + 2*S`` comes from the
    route's ``subtract=True`` path QGeometry row). Returns ``None`` when the route
    has no meander spacing option.
    """
    route = design.components[comp]
    meander = route.options.get('meander') if hasattr(route.options, 'get') else None
    spacing_raw = meander.get('spacing') if isinstance(meander, dict) else None
    if spacing_raw is None:
        return None
    spacing_um = _parse_mm(spacing_raw) * 1e3

    t = design.qgeometry.tables['path']
    cut_rows = t[(t['component'] == route.id) & (t['subtract'])]
    if len(cut_rows) > 0 and cut_rows.iloc[0]['width'] == cut_rows.iloc[0]['width']:
        cut_width_um = float(cut_rows.iloc[0]['width']) * 1e3
    else:
        w = _parse_mm(route.options.get('trace_width', design.variables.get('cpw_width', '10um'))) * 1e3
        g = _parse_mm(route.options.get('trace_gap', design.variables.get('cpw_gap', '6um'))) * 1e3
        cut_width_um = w + 2.0 * g

    trench_gap_um = spacing_um - cut_width_um
    if raise_on_violation and trench_gap_um + 1e-6 < float(min_trench_um):
        raise ValueError(
            f'{comp}: meander trench-to-trench gap {trench_gap_um:.2f} um '
            f'(spacing {spacing_um:.2f} um - cut {cut_width_um:.2f} um) '
            f'< {min_trench_um:g} um (DRC B7)'
        )
    return round(trench_gap_um, 2)


def solve_route_length(design, comp, target_mm, extra_length_fn=None, tol_mm=0.001,
                       max_iter=12, check_segments=True, check_meander=True,
                       min_trench_um=300.0):
    """Iterate ``comp.options.total_length`` until the filleted length matches ``target_mm``.

    Parameters
    ----------
    design : qiskit_metal.designs.DesignPlanar
        Design containing the already-instantiated route component ``comp``
        (e.g. ``RouteMixed`` or ``RouteMeander``).
    comp : str
        Route component name in ``design.components``.
    target_mm : float
        Desired total filleted centre-line length in mm (including any
        ``extra_length_fn(design)`` contribution, such as a coupled-line tee's
        second line).
    extra_length_fn : callable | None
        Optional ``fn(design) -> float`` returning additional length [mm] added
        to ``path_length(design, comp)`` when comparing against ``target_mm``.
    tol_mm : float
        Convergence tolerance in mm (default: 0.001 mm = 1 um).
    max_iter : int
        Maximum rebuild iterations (default: 12).
    check_segments : bool
        If True (default), verify DRC B2 (interior segments >= 2 * fillet) after
        convergence.
    check_meander : bool
        If True (default), verify DRC B7 (meander trench-to-trench >= ``min_trench_um``).
    min_trench_um : float
        Minimum allowed ground strip width between parallel meander trenches (default: 300 um).

    Returns
    -------
    dict
        ``{'component', 'target_mm', 'actual_mm', 'error_mm',
           'option_total_length', 'iterations', 'min_segment_mm',
           'meander_trench_gap_um'}``.
    """
    route = design.components[comp]
    target = float(target_mm)

    def measure():
        extra = float(extra_length_fn(design)) if extra_length_fn else 0.0
        return path_length(design, comp) + extra

    actual = measure()
    iters = 0
    for iters in range(1, max_iter + 1):
        err = target - actual
        if abs(err) <= tol_mm:
            break
        cur = _parse_mm(route.options.total_length)
        route.options.total_length = f'{cur + err:.4f}mm'
        route.rebuild()
        actual = measure()

    err = target - actual
    if abs(err) > tol_mm:
        raise RuntimeError(
            f'{comp}: failed to converge to {target:.4f} mm within {tol_mm:.4f} mm '
            f'after {max_iter} iterations (actual={actual:.4f} mm, err={err:+.4f} mm)'
        )

    min_seg = None
    if check_segments:
        min_seg, _ = check_min_segments(design, comp, raise_on_violation=True)

    trench_gap = None
    if check_meander:
        trench_gap = check_meander_spacing(design, comp, min_trench_um=min_trench_um, raise_on_violation=True)

    return {
        'component': comp,
        'target_mm': round(target, 6),
        'actual_mm': round(actual, 6),
        'error_mm': round(actual - target, 6),
        'option_total_length': str(route.options.total_length),
        'iterations': iters,
        'min_segment_mm': None if min_seg is None else round(min_seg, 6),
        'meander_trench_gap_um': trench_gap,
    }
