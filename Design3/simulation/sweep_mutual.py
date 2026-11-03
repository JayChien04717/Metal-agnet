"""Flux-line -> loop mutual inductance vs flux-line geometry (HFSS driven modal, hfss_mutual_fx.py).

Usage:  python sweep_mutual.py          # fl_sep 10 / 8 / 6 um, fl_sep 8 + fl_height 80 um
"""
import pandas as pd
import hfss_mutual_fx as m

ALL = {
    'mutual_fl10': dict(fl_sep='10um'),
    'mutual_fl8': dict(fl_sep='8um'),
    'mutual_fl6': dict(fl_sep='6um'),
    'mutual_fl8_h80': dict(fl_sep='8um', fl_height='80um'),
    'mutual_fl10_h80': dict(fl_sep='10um', fl_height='80um'),
    'mutual_fl12_h80': dict(fl_sep='12um', fl_height='80um'),
    'mutual_fl12_h100': dict(fl_sep='12um', fl_height='100um'),
}
import sys
CASES = {k: v for k, v in ALL.items() if k in sys.argv[1:]} or ALL
rows = []
for tag, fl in CASES.items():
    df, s = m.run(tag, fx_options=dict(flux_line_options=fl))
    rows.append(dict(tag=tag, **fl, M_pH=s['M_pH_fmin'], I_phi0_mA=s['I_phi0_mA'], L_fl_pH=s['L_fl_pH']))
    pd.DataFrame(rows).to_csv(m.os.path.join(m.RESULTS, 'sweep_mutual_%s.csv' % ('_'.join(CASES) if len(sys.argv) > 1 else 'all')), index=False)
print(pd.DataFrame(rows).round(4).to_string())
