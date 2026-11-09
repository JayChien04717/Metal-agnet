"""Use the already-solved eigenmode project open in HFSS: integrate |E| along short probe lines over the filter and
resonator vertical segments, for each mode, to tell which mode lives on which line (no re-solve).
Usage: python eigen_field_probe.py   (HFSS must still hold the solved e6_T3 design)"""
import os, sys, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
import numpy as np, pandas as pd
import design6_layout as d6
import pyEPR.ansys as a

d = d6.build_design()
y_c = float(d.components['cp'].pins[d6.CP_FPIN]['middle'][1]) if hasattr(d6, 'CP_FPIN') else float(d.components['cp'].pins['b']['middle'][1])
DX = 0.20
y_f = y_c - d6.VERT                                        # filter vertical segment: y_f .. y_c at x = 0
y_ro = y_c - d6.VERT - d6.OPEN_DY                          # resonator vertical segment: y_ro .. y_c-0.44 at x = DX
y_top = y_c - 0.44
lines = {'filt': (0.0 + 0.008, y_f, y_c), 'res': (DX + 0.008, y_ro, y_top)}
app = a.HfssApp(); desk = app.get_app_desktop(); proj = desk.get_active_project(); des = proj.get_active_design()
setup = des.get_setup(); sols = setup.get_solutions(); fields = setup.get_fields(); modeler = des.modeler
Z = 0.02
existing = set(modeler.get_objects_in_group('Lines')) if hasattr(modeler, 'get_objects_in_group') else set()
print('existing probe lines:', sorted(existing))
rows = []
for key, (x, ya, yb) in lines.items():
    ys = np.linspace(ya, yb, 5)
    for i in range(4):
        nm = f'probe_{key}_{i}'
        if nm in existing: continue
        modeler.draw_polyline([(f'{x}mm', f'{ys[i]}mm', f'{Z}mm'), (f'{x}mm', f'{ys[i + 1]}mm', f'{Z}mm')], closed=False, nonmodel=True, name=nm)
for n in range(1, 4):
    sols.set_mode(n, 0)
    for key in lines:
        vals = [fields.Mag_E.integrate_line(f'probe_{key}_{i}').evaluate(phase=0) for i in range(4)]
        rows.append(dict(mode=n, line=key, **{f'seg{i}': v for i, v in enumerate(vals)}, total=sum(vals)))
df = pd.DataFrame(rows); print(df.to_string(index=False))
df.to_csv(os.path.join(HERE, 'results', f"{sys.argv[1] if len(sys.argv) > 1 else 'e6_T3'}_field_probe.csv"), index=False)
