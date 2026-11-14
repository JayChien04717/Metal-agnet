"""Render the FX1 R5 mutual-inductance fixture for inspection without solving."""
from pathlib import Path
from shapely.geometry import box

from run_verification import rounded_design, save_aedt
from hfss_mutual_fx import (
    JJ_TO_PORT, OPEN_PINS, PORT_LIST, QUBIT, SETUP, ScatteringImpedanceSim,
)

PROJECT = Path(__file__).resolve().parent / "port_review" / "Design3_2_FX1_Port_Review.aedt"
DESIGN_NAME = "FX1_R5_Protected_Ports"

def main():
    if PROJECT.exists():
        raise FileExistsError(f"Open the existing review project: {PROJECT}")
    design = rounded_design(for_process=True, fx_options={"flux_as_ground": "False"})
    # The incoming U leg extends 1 um past its pin into the port endcap.
    # Trim that tip to the port plane while keeping the other leg's ground short.
    fx = design.components[QUBIT]
    pin = fx.pins["flux_line"]
    x, y = pin.middle
    half_width = pin.width / 2
    polys = design.qgeometry.tables["poly"]
    idx, = polys.index[(polys.component == fx.id) & (polys.name == "flux_line")]
    line = polys.at[idx, "geometry"]
    trimmed = line.difference(box(line.bounds[0] - .001, y - half_width - .001,
                                  x, y + half_width + .001))
    if trimmed.geom_type != "Polygon" or not trimmed.is_valid:
        raise ValueError("Flux line must remain a connected polygon")
    polys.at[idx, "geometry"] = trimmed
    save_aedt.use_project(PROJECT)
    sim = ScatteringImpedanceSim(design, "hfss")
    sim.setup.update(SETUP)

    # Build the existing mutual-inductance setup and sweep without Analyze.
    def setup_only():
        sim.sim_setup_name, sim.sweep_name = sim.renderer.initialize_drivenmodal(**sim.setup)

    sim._analyze = setup_only
    name, _ = sim.run_sim(
        name=DESIGN_NAME,
        components=[QUBIT],
        open_terminations=OPEN_PINS,
        port_list=PORT_LIST,
        jj_to_port=JJ_TO_PORT,
        box_plus_buffer=True,
    )
    save_aedt._desktop().get_active_project()._project.Save()
    inspect(name)


def inspect(name):
    desktop = save_aedt._desktop()._desktop
    project = desktop.SetActiveProject(PROJECT.stem)
    native = project.SetActiveDesign(name)
    editor = native.SetActiveEditor("3D Modeler")
    ports = list(native.GetModule("BoundarySetup").GetExcitations())
    print("EXCITATIONS", ports, flush=True)
    for expected in ("LumpPort_FX1_flux_line", "LumpPort_FX1_rect_jj"):
        if not any(port.split(":")[0] == expected for port in ports):
            raise RuntimeError(f"Missing expected excitation {expected}")
    for label in ("flux_line", "rect_jj"):
        sheet = f"Port_FX1_{label}"
        vertices = [editor.GetVertexPosition(v) for v in editor.GetVertexIDsFromObject(sheet)]
        print(sheet, vertices, flush=True)
    # Probe points just inside the intended conducting edges in native HFSS.
    for label, points in {
        "JJ": [(-.340, -.79895), (-.340, -.80105)],
        "flux line": [(-.46995, -.760), (-.47605, -.760)],
    }.items():
        for x, y in points:
            bodies = editor.GetBodyNamesByPosition(
                ["NAME:Parameters", "XPosition:=", f"{x}mm",
                 "YPosition:=", f"{y}mm", "ZPosition:=", "0mm"]
            )
            print("CONTACT", label, x, y, list(bodies), flush=True)
    print("VALIDATION", native.ValidateDesign(str(PROJECT.with_suffix(".validation.log"))),
          flush=True)
    save_aedt._desktop()._desktop.RestoreWindow()
    editor.FitAll()
    print(f"OPEN {PROJECT} / {name}", flush=True)

if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["--inspect"]:
        inspect(DESIGN_NAME + "_hfss")
    else:
        main()
