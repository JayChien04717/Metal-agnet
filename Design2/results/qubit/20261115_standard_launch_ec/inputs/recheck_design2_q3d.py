from datetime import datetime
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path.cwd()
design = root / "Design2"
runner = design / "simulation/run_q3d.py"
run = design / "simulation/results" / sys.argv[1]
inputs = run / "inputs"
raw = run / "raw"
analysis = run / "analysis"
raw.mkdir(parents=True, exist_ok=False)
analysis.mkdir()
sources = [
    runner, design / "pattern design/layout3_full.py", design / "gds/layout3 ftf.gds",
    design / "gds/ftf_layout.gds", root / "Design5/simulation/save_aedt.py",
    root / "pattern design/bondpad.py", root / "note/cap_matrix.py",
    Path(__file__),
]
hashes = {}
for source in sources:
    if source.exists():
        relative = source.relative_to(root) if source.is_relative_to(root) else Path(source.name)
        target = inputs / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = source.read_bytes()
        target.write_bytes(payload)
        hashes[relative.as_posix()] = sha256(payload).hexdigest()

project = design / "simulation/aedt/Design2.aedt"
setup_name = "D2_Launch_EC"
os.environ["AEDT_PROJECT"] = str(project)
manifest = {
    "run_id": run.name,
    "started": datetime.now().astimezone().isoformat(),
    "command": [sys.executable, *sys.argv],
    "working_directory": str(root),
    "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "dirty_worktree": bool(subprocess.check_output(["git", "status", "--short"], text=True).strip()),
    "python": sys.version,
    "solver": "Ansys Electronics Desktop 2021 R2 / Q3D Extractor",
    "input_sha256": hashes,
    "geometry_revision": "Design2 standard base metal with Design6 launch pads, 525 um Si",
    "aedt_project": str(project.relative_to(root)),
    "aedt_design": setup_name + "_q3d",
    "status": "running",
}

def save_manifest():
    (run / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

save_manifest()
try:
    source = runner.read_text(encoding="utf-8")
    source = source.replace("out = design_dir / 'simulation/results/layout3_full'", "out = run_output")
    source = source.replace("name='D2_Layout3'", f"name='{setup_name}',components=list(d.components)")
    exec(compile(source, str(runner), "exec"), {"__file__": str(runner), "run_output": raw})

    import numpy as np
    import pandas as pd
    sys.path.insert(0, str(root / "note"))
    import cap_matrix as cap

    cmat = pd.read_csv(raw / "capacitance_fF.csv", index_col=0)
    assert list(cmat.index) == list(cmat.columns)
    matrix = cmat.to_numpy()
    assert np.allclose(matrix, matrix.T, rtol=1e-3, atol=1e-4)
    labels = ("fx_left_bottom", "fx_left_top", "transmon", "fx_right_top",
              "fx_right_bottom", "res_left", "res_right")
    nodes = [cap.find_node(cmat, "metal_" + label) for label in labels]
    assert len(set(nodes)) == len(labels)
    assert min(np.linalg.eigvalsh(cap.reduced_matrix(cmat, nodes))) > 0
    result = cap.analyse(cmat, nodes)
    result["C_sigma_fF"] = cap.total_capacitances(cmat, nodes)
    result["node_mapping"] = dict(zip(labels, nodes))
    conv = pd.read_csv(raw / "q3d_convergence.csv")
    result["q3d_passes"] = len(conv)
    result["q3d_final_delta_pct"] = float(conv["Delta %"].iloc[-1])
    result["q3d_target_delta_pct"] = 0.5
    (analysis / "electrical.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    manifest["status"] = "completed"
    manifest["convergence"] = {key: result[key] for key in
                              ("q3d_passes", "q3d_final_delta_pct", "q3d_target_delta_pct")}
except Exception:
    manifest["status"] = "failed"
    raise
finally:
    manifest["finished"] = datetime.now().astimezone().isoformat()
    save_manifest()
