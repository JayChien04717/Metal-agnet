"""Isolated readout eigenmode fixtures from the restored Design2 mask."""
from pathlib import Path
import os
import sys

import numpy as np
from shapely.geometry import box, Point
from shapely import get_parts, prepare
from shapely.affinity import scale
from qiskit_metal import designs
from qiskit_metal.qlibrary.core import QComponent
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
from pyEPR.ansys import CalcObject
import pandas as pd

design_dir = Path(__file__).resolve().parents[1]
root = design_dir.parent
out = design_dir / "simulation/results/layout3_full/readout_eigenmode"
out.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(design_dir / "pattern design"), str(root / "Design5/simulation")]
from layout3_full import build_linear_design, recover
import save_aedt

src = design_dir / "gds/layout3 ftf.gds"
original = build_linear_design(src)
ground, nets = recover(src)

def fixture(side):
    sign = -1 if side == "left" else 1
    x0, x1 = sorted((sign * 1740, sign * 140))
    y0, y1 = -1280, 2660
    chip = box(x0, y0, x1, y1)
    # A grounded reference replaces the feedline in the isolated eigenmode fixture.
    # The original resonator CPW and its end clearances remain in the model.
    fixed_ground = ground.union(box(x0, 2481, x1, y1)).intersection(chip)
    signal = nets[f"res_{side}"]
    assert fixed_ground.distance(signal) > 9.9
    names = [f"res_{side}", f"fx_{side}_top", f"fx_{side}_bottom"]
    ids = {original.components[n].id for n in names}
    d = designs.DesignPlanar()
    d.chips.main.material = "silicon"
    d.chips.main.size.update(center_x=f"{(x0+x1)/2}um", center_y=f"{(y0+y1)/2}um",
                            size_x=f"{x1-x0}um", size_y=f"{y1-y0}um", size_z="-525um",
                            sample_holder_top="550um", sample_holder_bottom="1000um")
    class Fixture(QComponent):
        def make(self):
            for i, part in enumerate(get_parts(chip.difference(fixed_ground))):
                self.add_qgeometry("poly", {f"cut_{i}": scale(part, .001, .001, origin=(0,0))}, subtract=True)
            for _, row in original.qgeometry.tables["poly"].iterrows():
                if row.component in ids:
                    self.add_qgeometry("poly", {f"{row['name']}_{row.component}": row.geometry})
            for _, row in original.qgeometry.tables["junction"].iterrows():
                if row.component in ids:
                    self.add_qgeometry("junction", {row["name"]: row.geometry}, width=row.width,
                                       hfss_inductance=row.hfss_inductance)
    Fixture(d, "readout")
    prepare(signal)
    prepare(fixed_ground)
    bridges = []
    # Vertical straps cross each horizontal section of the meander. The two feet
    # must lie on ground throughout their width, with no second signal crossing.
    for x in (sign*685, sign*990, sign*1270):
        for y in range(-480, 2471, 1):
            if (signal.covers(box(x-5,y-2,x+5,y+2))
                    and fixed_ground.covers(box(x-5,y-34,x+5,y-27))
                    and fixed_ground.covers(box(x-5,y+27,x+5,y+34))
                    and all((x-a)**2+(y-b)**2 >= 150**2 for a,b in bridges)):
                bridges.append((x,y))
    assert len(bridges) >= 9, bridges
    regions = {"tee": (x0+30,x1-30,2360,2480),
               "meander": (x0+30,x1-30,-470,2360),
               "readout_pad": (sign*485-290,sign*485+290,-715,-480),
               "qubit": (x0+30,x1-30,-1260,-730)}
    return d, bridges, regions

def add_mesh_and_bridges(renderer, bridges):
    original_mesh = renderer.add_mesh
    def apply():
        original_mesh()
        modeler = renderer.modeler
        sheets = list(modeler._modeler.GetObjectsInGroup("Sheets"))
        metal = [s for s in sheets if s.startswith("metal_") or s.startswith("array_end")]
        assert metal, sheets
        modeler.mesh_length("ReadoutTraceSeed", metal, MaxLength="15um")
        straps = []
        for i,(x,y) in enumerate(bridges):
            for j,(yy,z,dy,dz) in enumerate(((y-30,0,0,10),(y-30,10,60,0),(y+30,0,0,10))):
                name = f"bridge_{i}_{j}"
                modeler.draw_rect_corner([(x-5)*1e-6,yy*1e-6,z*1e-6],
                    x_size=10e-6,y_size=dy*1e-6,z_size=dz*1e-6,name=name)
                straps.append(name)
        modeler.assign_perfect_E(straps, name="ReadoutGroundStraps")
    renderer.add_mesh = apply

def probe(sim, side, regions):
    setup = sim.renderer.pinfo.setup
    modeler = sim.renderer.modeler
    for name,(x0,x1,y0,y1) in regions.items():
        modeler.draw_box_corner([x0*1e-6,y0*1e-6,-120e-6],
            [(x1-x0)*1e-6,(y1-y0)*1e-6,220e-6],
            name="E_"+name,nonmodel=True,wireframe=True)
    sol = setup.get_solutions()
    f, _ = sol.eigenmodes()
    data = []
    for mode, freq in enumerate(f,1):
        sol.set_mode(mode)
        energies = {}
        for name in regions:
            setup.parent._fields_calc.CalcStack("Clear")
            electric = CalcObject([],setup).getQty("E")
            energies[name] = max(0,float(electric.dot(electric.conj()).real().integrate_vol("E_"+name).evaluate()))
        total = sum(energies.values())
        data.append(dict(mode=mode,frequency_GHz=freq,**{n:v/total for n,v in energies.items()}))
        print(side, data[-1], flush=True)
    pd.DataFrame(data).to_csv(out/f'local_eigen_{side}_fields.csv',index=False)

if __name__ == "__main__":
    sides = [arg for arg in sys.argv[1:] if arg in ("left", "right")] or ["left","right"]
    check = "--check" in sys.argv
    if not check:
        save_aedt.use_project(Path(os.environ.get("AEDT_PROJECT", "D:/HFSS_temp/D2_Layout3.aedt")))
        os.chdir(out)
    for side in sides:
        d, bridges, regions = fixture(side)
        print(side, bridges, flush=True)
        pd.DataFrame(bridges,columns=["x_um","y_um"]).to_csv(out/f"local_eigen_{side}_bridges.csv",index=False)
        if check:
            continue
        sim = EigenmodeSim(d,"hfss")
        add_mesh_and_bridges(sim.renderer, bridges)
        sim.setup.update(dict(min_freq_ghz=5.7,n_modes=3,max_delta_f=.30,max_passes=6,
                             min_passes=3,min_converged=2,pct_refinement=20,basis_order=1))
        sim.run_sim(name=f"D2_Readout_Local_{side}",components=list(d.components),open_terminations=[],box_plus_buffer=False)
        sim.convergence_f.to_csv(out/f"local_eigen_{side}_freqs.csv")
        sim.convergence_t.to_csv(out/f"local_eigen_{side}_mesh.csv")
        save_aedt._desktop().get_active_project()._project.Save()
        try:
            probe(sim,side,regions)
        except Exception:
            import traceback
            traceback.print_exc()
