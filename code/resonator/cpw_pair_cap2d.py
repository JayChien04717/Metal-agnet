"""2-D quasi-static capacitance of one CPW or two parallel CPWs on a substrate (finite-volume Laplace solve).

Used to estimate the distributed filter-resonator parasitic c'_m (fF/mm) of two side-by-side CPW runs as a
function of centre spacing, for the Purcell + readout module (pattern design/resonator/purcell_readout.py,
note/purcell_filter_method.md section 6). It is a cross-section (infinitely long lines) solve: no ends, tees, hammer
pocket, bonds or bridges.

    python -I -B code/resonator/cpw_pair_cap2d.py --out <run dir>

The rows the model uses are pasted into purcell_readout_model.CPW_PAIR_2D.

Model
* Zero-thickness metal on the substrate top face (y = 0): signal strips W, gaps S, ground everywhere else up to
  the domain edge (`outer_ground` beyond the outermost gap).
* Substrate eps_r (default 11.9) of thickness h (default 500 um) between y = -h and 0; vacuum above and below.
  Outer box boundaries are natural (zero normal field); every field line ends on the grounds.
* Tensor-product grid graded geometrically from h0 at every metal edge and at y = 0 / y = -h.
* Charge on a conductor = sum of the finite-volume flux residual over its nodes. Maxwell matrix in fF/mm.

Validation printed by the script: single 10/6 um CPW vs the conformal-mapping value in skill/chip_drc.md A.1
(163.4 fF/mm, eps_eff 6.225) and the 16/9 um pair at 200 um vs the Design4 Q3D pair coupons
(C'_m = 0.563 fF/mm, ground C' 165.8 fF/mm at eps_r 11.5; Design4/notes/hammer_fixed_length_validation_20261116.md).
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import platform
import sys

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

EPS0_FF_PER_MM = 8.8541878128e-12 * 1e15 / 1e3     # F/m -> fF/mm


def graded(a, b, h0, ratio, hmax):
    """Points from a to b, step h0 at both ends growing by `ratio` up to hmax."""
    if b - a <= 2 * h0:
        return np.array([a, b])
    left, right, h = [a], [b], h0
    while right[-1] - left[-1] > 2 * h:
        left.append(left[-1] + h)
        right.append(right[-1] - h)
        h = min(h * ratio, hmax)
    pts = np.array(left + right[::-1])
    return np.unique(np.round(pts, 9))


def axis(keys, h0, ratio, hmax):
    keys = np.unique(np.round(np.asarray(keys, float), 9))
    pts = [graded(a, b, h0, ratio, hmax) for a, b in zip(keys[:-1], keys[1:])]
    return np.unique(np.concatenate(pts))


def solve(centres_um, w_um, s_um, *, eps_r=11.9, h_um=500.0, outer_ground_um=1500.0, air_um=1000.0,
          h0_um=0.25, ratio=1.15, hmax_um=25.0):
    """Maxwell capacitance matrix (fF/mm) of the signal strips centred at `centres_um` (ground is the reference)."""
    centres = sorted(centres_um)
    xmin = centres[0] - w_um / 2 - s_um - outer_ground_um
    xmax = centres[-1] + w_um / 2 + s_um + outer_ground_um
    edges = []
    for c in centres:
        edges += [c - w_um / 2 - s_um, c - w_um / 2, c + w_um / 2, c + w_um / 2 + s_um]
    x = axis([xmin, xmax, *edges], h0_um, ratio, hmax_um)
    y = axis([-h_um - air_um, -h_um, 0.0, air_um], h0_um, ratio, hmax_um)
    nx, ny = len(x), len(y)
    j0 = int(np.argmin(abs(y)))
    assert abs(y[j0]) < 1e-9
    dx, dy = np.diff(x), np.diff(y)
    yc = (y[:-1] + y[1:]) / 2
    eps_cell_row = np.where((yc < 0) & (yc > -h_um), eps_r, 1.0)        # per cell row, same for all columns

    idx = np.arange(nx * ny).reshape(nx, ny)
    rows, cols, vals = [], [], []

    def add(a, b, wgt):
        rows.extend([a, b]); cols.extend([b, a]); vals.extend([wgt, wgt])

    def flat(parts):
        return np.concatenate([np.ravel(p) for p in parts])

    # horizontal edges (i,j)-(i+1,j): half cells above (row j) and below (row j-1)
    up = np.zeros(ny); dn = np.zeros(ny)
    up[:-1] = eps_cell_row * dy / 2
    dn[1:] = eps_cell_row * dy / 2
    hw = (up + dn)[None, :] / dx[:, None]                               # (nx-1, ny)
    a = idx[:-1, :].ravel(); b = idx[1:, :].ravel()
    add(a, b, hw.ravel())
    # vertical edges (i,j)-(i,j+1): half columns left and right, eps of row j
    lf = np.zeros(nx); rt = np.zeros(nx)
    rt[:-1] = dx / 2
    lf[1:] = dx / 2
    vw = (lf + rt)[:, None] * eps_cell_row[None, :] / dy[None, :]       # (nx, ny-1)
    a = idx[:, :-1].ravel(); b = idx[:, 1:].ravel()
    add(a, b, vw.ravel())
    W = coo_matrix((flat(vals), (flat(rows), flat(cols))), shape=(nx * ny, nx * ny)).tocsr()
    A = scipy.sparse.diags(np.asarray(W.sum(axis=1)).ravel()) - W

    # conductors on row j0
    label = np.full(nx, -1)            # -1 dielectric gap, 0 ground, k>=1 signal k
    tol = 1e-9
    label[:] = 0
    for k, c in enumerate(centres, 1):
        gap = (x > c - w_um / 2 - s_um + tol) & (x < c + w_um / 2 + s_um - tol)
        sig = (x >= c - w_um / 2 - tol) & (x <= c + w_um / 2 + tol)
        label[gap] = -1
        label[sig] = k
    fixed = np.zeros(nx * ny, bool)
    fixed[idx[label >= 0, j0]] = True
    free = ~fixed
    A_ff = A[free][:, free].tocsc()
    A_fd = A[free][:, fixed]
    n = len(centres)
    C = np.zeros((n, n))
    for k in range(1, n + 1):
        phi_d = np.zeros(fixed.sum())
        fixed_labels = np.full(nx * ny, -9)
        fixed_labels[idx[:, j0]] = label
        phi_d[fixed_labels[fixed] == k] = 1.0
        phi = np.zeros(nx * ny)
        phi[fixed] = phi_d
        phi[free] = spsolve(A_ff, -A_fd @ phi_d)
        q = A @ phi
        for m in range(1, n + 1):
            C[m - 1, k - 1] = q[idx[label == m, j0]].sum() * EPS0_FF_PER_MM
    return C, dict(nx=nx, ny=ny, h0_um=h0_um, ratio=ratio, hmax_um=hmax_um)


def single(w, s, **kw):
    c, grid = solve([0.0], w, s, **kw)
    c_air, _ = solve([0.0], w, s, **dict(kw, eps_r=1.0))
    return c[0, 0], c[0, 0] / c_air[0, 0], grid


def pair(spacing, w, s, **kw):
    """Two CPWs at `spacing`: Maxwell self / mutual C, and the quasi-TEM coupling factors
    k_C = C'_m / C'_self and k_L = L'_m / L'_self (L = C_air^-1 / c^2, so k_L from the vacuum matrix)."""
    c, grid = solve([-spacing / 2, spacing / 2], w, s, **kw)
    c_air, _ = solve([-spacing / 2, spacing / 2], w, s, **dict(kw, eps_r=1.0))
    l = np.linalg.inv(c_air)
    return dict(spacing_um=spacing, c_self_fF_per_mm=float(c[0, 0]), c_mutual_fF_per_mm=float(-c[1, 0]),
                c_ground_fF_per_mm=float(c[0, 0] + c[1, 0]), k_C=float(-c[1, 0] / c[0, 0]),
                k_L=float(l[1, 0] / l[0, 0]), grid=grid)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True, help='output directory (a new run directory)')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    c0 = 299792458.0
    out = dict(tool='code/resonator/cpw_pair_cap2d.py', created=datetime.datetime.now().isoformat(timespec='seconds'),
               code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               environment=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                                platform=platform.platform()),
               model='2-D quasi-static finite-volume, zero-thickness metal, 500 um substrate, vacuum outside, '
                     'natural outer boundary, 1500 um outer ground', validation={}, sweeps={})

    # grid check and validation
    conv = {}
    for h0 in (1.0, 0.5, 0.25):
        c, ee, _ = single(10, 6, eps_r=11.9, h0_um=h0)
        conv[str(h0)] = dict(c_fF_per_mm=c, eps_eff=ee)
    c, ee, grid = single(10, 6, eps_r=11.9)
    z0 = 1 / (c * 1e-15 / 1e-3 * c0 / np.sqrt(ee))
    out['validation']['cpw_10_6_eps11.9'] = dict(
        c_fF_per_mm=c, eps_eff=ee, v_m_per_s=c0 / np.sqrt(ee), z0_ohm=z0, grid=grid, h0_sweep=conv,
        reference='conformal mapping, thick substrate, eps_r 11.9: 164.4 fF/mm, eps_eff 6.45. skill/chip_drc.md A.1 '
                  'quotes 163.4 fF/mm, eps_eff 6.225 (= (1 + 11.45) / 2, i.e. eps_r 11.45)')
    cm_conv = {str(h0): pair(400, 10, 6, eps_r=11.9, h0_um=h0)['c_mutual_fF_per_mm'] for h0 in (1.0, 0.5, 0.25)}
    out['validation']['cm_10_6_400um_h0_sweep'] = cm_conv
    print('c_m(400 um, 10/6) vs h0:', cm_conv)
    p = pair(200, 16, 9, eps_r=11.5)
    out['validation']['pair_16_9_eps11.5_200um'] = dict(
        p, reference='Design4 Q3D pair coupons (1 mm / 2 mm difference): C\'m 0.563 fF/mm, ground C\' 165.8 / 165.9 fF/mm')
    print('10/6 single CPW (eps 11.9):', json.dumps(out['validation']['cpw_10_6_eps11.9'], default=float)[:400])
    print('16/9 pair 200 um (eps 11.5):', p)

    for name, (w, s, eps) in {'10_6_eps11.9': (10, 6, 11.9), '16_9_eps11.9': (16, 9, 11.9),
                              '16_9_eps11.5': (16, 9, 11.5)}.items():
        rows = []
        for spacing in (150, 200, 250, 300, 400, 500, 600, 800):
            r = pair(spacing, w, s, eps_r=eps)
            rows.append(r)
            print(name, spacing, 'um: c_m', round(r['c_mutual_fF_per_mm'], 4), 'fF/mm')
        out['sweeps'][name] = rows
    path = args.out / 'cpw_pair_cap2d.json'
    path.write_text(json.dumps(out, indent=1, default=float))
    print('wrote', path)


if __name__ == '__main__':
    main()
