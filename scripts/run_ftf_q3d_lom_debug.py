"""Run and debug FTF Q3D capacitance extraction from FTF.ipynb.

This script executes the geometry cells from FTF.ipynb, then executes the
notebook's LOM/Q3D cells. It is intentionally notebook-driven so the debug
path matches the interactive workflow.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import traceback


WORKSPACE = Path(__file__).resolve().parents[1]
NOTEBOOK = WORKSPACE / "FTF.ipynb"


def _cell_source(nb: dict, idx: int) -> str:
    return "".join(nb["cells"][idx].get("source", []))


def _notebook_code_without_magics(source: str) -> str:
    lines: list[str] = []
    for line in source.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("%"):
            continue
        if "MetalGUI" in line:
            continue
        lines.append(line)
    return "\n".join(lines)


def run() -> None:
    os.chdir(WORKSPACE)
    if str(WORKSPACE) not in sys.path:
        sys.path.insert(0, str(WORKSPACE))

    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8-sig"))
    ns: dict[str, object] = {"__name__": "__ftf_q3d_lom_debug__"}

    geometry_chunks = [
        _notebook_code_without_magics(_cell_source(nb, idx)) for idx in range(1, 9)
    ]
    exec(compile("\n\n".join(geometry_chunks), str(NOTEBOOK), "exec"), ns)

    for idx in (14, 15, 16):
        print(f"\n===== executing FTF.ipynb cell {idx} =====")
        if idx == 15:
            install_renderer_path_trace(ns["design"])
        exec(compile(_cell_source(nb, idx), f"{NOTEBOOK} cell {idx}", "exec"), ns)


def install_renderer_path_trace(design) -> None:
    """Print the QGeometry path being rendered before AEDT sweep calls."""
    from qiskit_metal.renderers.renderer_ansys.ansys_renderer import QAnsysRenderer

    install_ground_subtract_trace()

    if getattr(QAnsysRenderer.render_element_path, "_ftf_debug_wrapped", False):
        return

    original = QAnsysRenderer.render_element_path
    component_by_id = {component.id: name for name, component in design.components.items()}

    def traced_render_element_path(self, qgeom):
        component_name = component_by_id.get(int(qgeom.component), qgeom.component)
        geometry = qgeom.geometry
        length = getattr(geometry, "length", None)
        print(
            "RENDER PATH",
            f"component={component_name}",
            f"name={qgeom['name']}",
            f"subtract={bool(qgeom['subtract'])}",
            f"length={length}",
            flush=True,
        )
        return original(self, qgeom)

    traced_render_element_path._ftf_debug_wrapped = True
    QAnsysRenderer.render_element_path = traced_render_element_path


def install_ground_subtract_trace() -> None:
    """Subtract ground cutouts one at a time and print the failing body."""
    from qiskit_metal.renderers.renderer_ansys.ansys_renderer import QAnsysRenderer

    if getattr(QAnsysRenderer.subtract_from_ground, "_ftf_debug_wrapped", False):
        return

    def traced_subtract_from_ground(self):
        import pythoncom

        for chip, shapes in self.chip_subtract_dict.items():
            for shape in sorted(shapes):
                print(
                    "SUBTRACT GROUND",
                    f"chip={chip}",
                    f"shape={shape}",
                    flush=True,
                )
                try:
                    self.modeler.subtract(f"ground_{chip}_plane", [shape])
                except pythoncom.com_error as error:
                    print("FAILED SUBTRACT", f"chip={chip}", f"shape={shape}")
                    print("com_error:", error)
                    raise

    traced_subtract_from_ground._ftf_debug_wrapped = True
    QAnsysRenderer.subtract_from_ground = traced_subtract_from_ground


if __name__ == "__main__":
    try:
        run()
    except Exception as exc:
        print("\nFTF Q3D/LOM debug run failed.")
        print(f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
        raise
