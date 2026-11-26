# -*- coding: utf-8 -*-
r"""Full-feed kappa_p estimate for the Design11 readout feed with a series C_in (skill/purcell_resonator_design.md section 4).

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/feed_kappa_model.py \
        --coupons <results/resonator/<q3d_coupons run>> --cp Q1=180,Q2=220 --cin 140 [--overrides JSON]

Model (circuit only, no EM solve): 50 Ohm port P8 - feed_cin - C_in pi network (Q3D coupon: series mutual, both shunt grounds)
- feed_in - tee n0_Q1 - feed_mid - tee n0_Q2 - feed_out - 50 Ohm port P7. Line lengths are the filleted route lengths of
build.py; each tee carries a shunt (stub line + C_p frame-side ground). For channel k:
    kappa_p,k ~ kappa_hanger(f_k, C_p,k) * Re Z_node,k / (Z0 / 2)
with Z_node,k the feed impedance at its tee. Robustness corners: f x (0.99, 1, 1.01), v x (0.98, 1.02),
C_frame x (0.8, 1.2), C_in x (0.9, 1.1); the worst-case ratio Re Z_node / (Z0/2) is reported.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout')]
import project_paths
project_paths.configure(__file__)

import numpy as np
import build as D11
import purcell_readout_model as PRM

LN = PRM.line('10/6')
Z0, V = LN['z0'], LN['v']


def route_len_mm(d, name):
    t = d.qgeometry.tables['path']
    cid = d.components[name].id
    rows = t[(t['component'] == cid) & (~t['subtract'].astype(bool))]
    return float(sum(g.length for g in rows['geometry']))


def feed_lengths(d):
    half_tee = 0.150
    stub = route_len_mm(d, 'stub_p_Q1') + 0.100          # stub_p route + TapJunction stub_length
    return dict(a=route_len_mm(d, 'feed_cin'),
                b=route_len_mm(d, 'feed_in') + half_tee,
                c=half_tee + route_len_mm(d, 'feed_mid') + half_tee,
                dd=half_tee + route_len_mm(d, 'feed_out'),
                stub=stub)


def abcd_line(L_mm, f, v):
    b = 2 * np.pi * f / v * L_mm * 1e-3
    return np.array([[np.cos(b), 1j * Z0 * np.sin(b)], [1j * np.sin(b) / Z0, np.cos(b)]])


def abcd_shunt(Y):
    return np.array([[1, 0], [Y, 1]], complex)


def abcd_series(Z):
    return np.array([[1, Z], [0, 1]], complex)


def z_looking(chain, z_term):
    """Impedance looking into the input of `chain` (list of ABCD, input first) terminated by z_term."""
    M = np.eye(2, dtype=complex)
    for m in chain:
        M = M @ m
    A, B, C, D = M.ravel()
    return (A * z_term + B) / (C * z_term + D)


def node_z(f, L, c_in, c_frame, v, stub_mm):
    """Feed impedance at tee 1 (Q1) and tee 2 (Q2). c_in: dict(series, ga, gb) fF; c_frame fF per tee."""
    w = 2 * np.pi * f
    Y = lambda c_fF: 1j * w * c_fF * 1e-15
    # tee shunt: open stub line + frame ground C at its end -> admittance of the stub branch
    zl_stub = 1 / Y(c_frame)
    y_tee = 1 / z_looking([abcd_line(stub_mm, f, v)], zl_stub)
    # left side of tee 1: C_in (b side shunt, series, a side shunt), line a, port
    left1 = z_looking([abcd_line(L['b'], f, v), abcd_shunt(Y(c_in['gb'])), abcd_series(1 / Y(c_in['series'])),
                       abcd_shunt(Y(c_in['ga'])), abcd_line(L['a'], f, v)], 50.0)
    right2 = z_looking([abcd_line(L['dd'], f, v)], 50.0)
    # right of tee 1: line c to tee 2 (with tee 2 shunt) then right2
    z2_right_with_shunt = 1 / (1 / right2 + y_tee)
    right1 = z_looking([abcd_line(L['c'], f, v)], z2_right_with_shunt)
    z1 = 1 / (1 / left1 + 1 / right1)
    left1_with_shunt = 1 / (1 / left1 + y_tee)
    left2 = z_looking([abcd_line(L['c'], f, v)], left1_with_shunt)
    z2 = 1 / (1 / left2 + 1 / right2)
    return z1, z2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--coupons', required=True, help='q3d_coupons run folder (summaries/summary.json)')
    ap.add_argument('--cp', default='Q1=140,Q2=140', help='plate width per channel (um), must exist in coupons')
    ap.add_argument('--cin', default='140', help='C_in plate width (um)')
    ap.add_argument('--overrides', default='{}')
    a = ap.parse_args()
    sp = Path(a.coupons) / 'summaries' / 'summary.json'
    if sp.exists():
        cs = json.loads(sp.read_text(encoding='utf-8'))['cp']
    else:   # interrupted coupon run: read the per-case raw results
        cs = {}
        for f in sorted((Path(a.coupons) / 'raw').glob('cp_w*/caps.json')):
            r = json.loads(f.read_text(encoding='utf-8'))
            cs[f'{float(f.parent.name[4:].replace("p", ".")):g}'] = dict(mutual_fF=r['mutual_fF'], ground_fF=r['ground_fF'])
    d = D11.build_design(variant='sim', overrides=json.loads(a.overrides), channels=True, check_keepout_um=None)
    L = feed_lengths(d)
    cin_c = cs[a.cin]
    c_in = dict(series=cin_c['mutual_fF'], ga=cin_c['ground_fF']['frame'], gb=cin_c['ground_fF']['plate'])
    cp_w = dict(kv.split('=') for kv in a.cp.split(','))
    out = dict(lengths_mm=L, c_in_fF=c_in, channels={})
    for i, q in enumerate(('Q1', 'Q2')):
        f0 = D11.PURCELL_CHANNELS[q]['f_GHz'] * 1e9
        cpq = cs[cp_w[q]]
        c_frame = cpq['ground_fF']['frame']
        k_h = PRM.kappa_hanger(f0 / 1e9, cpq['mutual_fF'])
        ratios = []
        for fs, vs, cfs, cis in itertools.product((0.99, 1.0, 1.01), (0.98, 1.0, 1.02), (0.8, 1.0, 1.2), (0.9, 1.0, 1.1)):
            ci = {k: v * cis for k, v in c_in.items()}
            z = node_z(f0 * fs, L, ci, c_frame * cfs, V * vs, L['stub'])[i]
            ratios.append(float(z.real / (Z0 / 2)))
        z_nom = node_z(f0, L, c_in, c_frame, V, L['stub'])[i]
        r_nom = float(z_nom.real / (Z0 / 2))
        out['channels'][q] = dict(plate_width_um=cp_w[q], C_p_fF=cpq['mutual_fF'], C_frame_fF=c_frame,
                                  kappa_hanger_MHz=k_h, ReZ_ratio_nominal=r_nom, kappa_p_model_MHz=k_h * r_nom,
                                  ReZ_ratio_min=min(ratios), ReZ_ratio_max=max(ratios),
                                  kappa_p_min_MHz=k_h * min(ratios), kappa_p_max_MHz=k_h * max(ratios))
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
