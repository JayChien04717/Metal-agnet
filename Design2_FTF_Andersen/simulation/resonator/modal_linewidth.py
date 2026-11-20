"""Calculate the modal capacitance referred to the coupling end."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import sys
import math
import pandas as pd
from pyEPR.ansys import CalcObject
from linewidth import cpw_impedance, read_matrix

design_dir = _P('Design2_FTF_Andersen')
root = _P('.')
sys.path.insert(0, _I('Design5_6_purcell-res/simulation'))
import save_aedt

out = _P('Design2_FTF_Andersen/results/general/legacy/layout3_full/readout_eigenmode')
out.mkdir(parents=True, exist_ok=True)
project = next(p for p in save_aedt._desktop().get_projects() if p.name == "D2_Layout3")
matrix = read_matrix(_P('Design2_FTF_Andersen/results/general/legacy/layout3_full/capacitance_fF.csv'))
rows = []
for side, sign, mode in (("left", -1, 2), ("right", 1, 1)):
    design = project.get_design(f"D2_Readout_Local_{side}_hfss")
    setup = design.get_setup("Setup")
    modeler = design.modeler
    x0 = -1740 if sign < 0 else 140
    sol = setup.get_solutions()
    frequencies, _ = sol.eigenmodes()
    sol.set_mode(mode)
    f = frequencies[mode-1]
    omega = 2*math.pi*f*1e9
    domain = "StoredEnergyDomain"
    if domain not in modeler._modeler.GetObjectsInGroup("Solids"):
        modeler.draw_box_corner([x0*1e-6,-1280e-6,-1000e-6],
                            [1600e-6,3940e-6,1550e-6],
                            name=domain,nonmodel=True,wireframe=True)
    calc = design._fields_calc
    calc.CalcStack("Clear")
    electric = CalcObject([],setup).getQty("E")
    stored = float(electric.times_eps().dot(electric.conj()).real().integrate_vol(domain).evaluate())
    h = CalcObject([],setup).getQty("H")
    calc.CalcStack("Clear")
    magnetic = float(h.times_mu().dot(h.conj()).real().integrate_vol(domain).evaluate())
    lumped_magnetic = 0.0
    x = sign*485
    arm = sign*539
    for label,start,end,inductance in (
        ("JJ",(x,-999),(x,-1001),34.38e-9),
        ("ArrayTop",(x+sign,-996),(arm-sign,-996),100e-9),
        ("ArrayBottom",(x+sign,-1004),(arm-sign,-1004),100e-9),
    ):
        name = "Energy"+label
        modeler.draw_polyline([[f"{xx}um",f"{yy}um","0um"] for xx,yy in (start,end)],
                               closed=False,name=name,nonmodel=True,wireframe=True)
        calc.CalcStack("Clear")
        vr = float(CalcObject([],setup).getQty("E").real().integrate_line_tangent(name).evaluate())
        calc.CalcStack("Clear")
        vi = float(CalcObject([],setup).getQty("E").imag().integrate_line_tangent(name).evaluate())
        lumped_magnetic += (vr**2+vi**2)/(omega**2*inductance)
    voltages = []
    for x in (520, 600, 750, 825, 920):
        name = f"CouplerVoltage{x}"
        if name not in modeler._modeler.GetObjectsInGroup("Lines"):
            modeler.draw_polyline([[f"{sign*x}um","2470.5um","0um"],
                               [f"{sign*x}um","2480.5um","0um"]],
                             closed=False,name=name,nonmodel=True,wireframe=True)
        calc.CalcStack("Clear")
        vr = float(CalcObject([],setup).getQty("E").real().integrate_line_tangent(name).evaluate())
        calc.CalcStack("Clear")
        vi = float(CalcObject([],setup).getQty("E").imag().integrate_line_tangent(name).evaluate())
        voltages.append(complex(vr,vi))
    voltage = sum(voltages)/len(voltages)
    c_mode = stored / abs(voltage)**2
    c_coupling = -float(matrix[f"metal_res_{side}"])*1e-15
    r = cpw_impedance()/2
    g = (omega*c_coupling)**2*r/(1+(omega*c_coupling*r)**2)
    rows.append(dict(side=side,f_GHz=f,Cmode_energy_fF=c_mode*1e15,
                     kappa_ext_mode_MHz=g/c_mode/(2*math.pi)/1e6,
                     field_E_total_normalized=stored,field_H_total_normalized=magnetic,
                     mean_gap_voltage_V=abs(voltage),
                     magnetic_electric_ratio=(magnetic+lumped_magnetic)/stored,
                     lumped_inductive_fraction=lumped_magnetic/(magnetic+lumped_magnetic),
                     gap_voltage_spread_percent=(max(map(abs,voltages))-min(map(abs,voltages)))/abs(voltage)*100))
    print(rows[-1],flush=True)
    pd.DataFrame(rows).to_csv(_P('Design2_FTF_Andersen/results/general/legacy/layout3_full/readout_eigenmode/modal_linewidth.csv'),index=False)
