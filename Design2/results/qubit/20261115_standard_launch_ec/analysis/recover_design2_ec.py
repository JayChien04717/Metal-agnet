from datetime import datetime
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

root = Path.cwd()
run = root / "Design2/results/qubit/20261115_standard_launch_ec"
old = root / "Design2/simulation/results/20261115_standard_launch_ec"
(run / "analysis").mkdir(exist_ok=True)
(run / "raw").mkdir(exist_ok=True)
for filename in ("capacitance_fF.csv", "q3d_convergence.csv", "solver.log"):
    if not (run / "raw" / filename).exists():
        shutil.copyfile(old / "raw" / filename, run / "raw" / filename)
cap_source = run / "inputs/note/cap_matrix.py"
spec = importlib.util.spec_from_file_location("cap_matrix", cap_source)
cap = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cap)
cmat = pd.read_csv(run / "raw/capacitance_fF.csv", index_col=0)
assert list(cmat.index) == list(cmat.columns)
matrix = cmat.to_numpy()
assert np.allclose(matrix, matrix.T, rtol=1e-3, atol=1e-4)
assert np.all(matrix[~np.eye(len(matrix), dtype=bool)] <= 0)
labels = ("fx_left_bottom", "fx_left_top", "transmon", "fx_right_top",
          "fx_right_bottom", "res_left", "res_right")
nodes = [cap.find_node(cmat, "metal_" + label) for label in labels]
assert len(set(nodes)) == len(labels)
assert min(np.linalg.eigvalsh(cap.reduced_matrix(cmat, nodes))) > 0
result = cap.analyse(cmat, nodes)
result["C_sigma_fF"] = cap.total_capacitances(cmat, nodes)
result["node_mapping"] = dict(zip(labels, nodes))
conv = pd.read_csv(run / "raw/q3d_convergence.csv")
result["q3d_passes"] = len(conv)
result["q3d_final_delta_pct"] = float(conv["Delta %"].iloc[-1])
result["q3d_target_delta_pct"] = 0.5
result["q3d_converged"] = result["q3d_final_delta_pct"] <= result["q3d_target_delta_pct"]
(run / "analysis" / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
(run / "analysis/electrical.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
manifest_path = run / "manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
manifest["status"] = "completed"
manifest["analysis_recovered"] = {
    "at": datetime.now().astimezone().isoformat(),
    "reason": "The solve and raw export completed before an import error. Analysis used the retained matrix.",
    "entry_point": f"analysis/{Path(__file__).name}",
}
manifest["convergence"] = {key: value for key, value in result.items() if key.startswith("q3d_")}
manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
