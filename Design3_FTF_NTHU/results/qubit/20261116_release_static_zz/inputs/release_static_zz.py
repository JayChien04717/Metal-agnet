import json
from pathlib import Path
import sys
from functools import lru_cache
from hashlib import sha256

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq, linear_sum_assignment, minimize_scalar

root, output = map(Path, sys.argv[1:3])
repo_model = root / "Design3/simulation/qubit/zz_analysis.py"
sys.path.insert(0, str(repo_model.parent))
from zz_analysis import fluxonium, transmon

sources = {
    "Design2": "Design2/results/qubit/20261115_standard_launch_ec/analysis/electrical.json",
    "Design3_3": "Design3_3/results/general/legacy/readout_flux_parameters.json",
}

class Circuit:
    def __init__(self, electric, nf, nt):
        self.ec = electric["Ec_TC [GHz]"]
        e1, n1 = fluxonium(electric["Ec_FX1 [GHz]"], 4.0, .6, n_keep=nf)
        e2, n2 = fluxonium(electric["Ec_FX2 [GHz]"], 4.0, .6, n_keep=nf)
        self.j = [electric[key] for key in ("J_1c [MHz]", "J_2c [MHz]", "J_12 [MHz]")]
        self.nf, self.nt = nf, nt
        i = np.eye(nf)
        t = np.eye(nt)
        self.e0 = (e1[:, None, None] + e2[None, None, :]) + np.zeros((1, nt, 1))
        self.static = 1e-3 * self.j[2] * np.kron(np.kron(n1, t), n2)
        self.n1, self.n2, self.i = n1, n2, i
        self.bare = [0, 1, nt*nf, nt*nf+1]

    def at(self, ej):
        et, nc = transmon(self.ec, ej, self.nt)
        h = (self.static + np.diag((self.e0 + et[None, :, None]).ravel())
             + 1e-3*(self.j[0]*np.kron(np.kron(self.n1,nc),self.i)
                     + self.j[1]*np.kron(np.kron(self.i,nc),self.n2)))
        w, v = eigh(h, subset_by_index=(0, 11))
        overlaps = abs(v[self.bare])**2
        indices = linear_sum_assignment(-overlaps)[1]
        e00, e01, e10, e11 = w[indices]
        return dict(EJ_eff_GHz=ej, bare_TC_f01_GHz=float(et[1]),
                    ZZ_kHz=float(1e6*(e11-e10-e01+e00)),
                    min_bare_overlap=float(min(overlaps[np.arange(4), indices])))

def optimum(function):
    grid = np.linspace(0., .5, 101)
    rows = [function(float(x)) for x in grid]
    zz = [row["ZZ_kHz"] for row in rows]
    candidates = [rows[0], rows[-1]]
    for a, b, fa, fb in zip(grid, grid[1:], zz, zz[1:]):
        if fa * fb < 0:
            candidates.append(function(brentq(lambda x: function(float(x))["ZZ_kHz"], a, b, xtol=1e-9)))
    for i, val in enumerate(zz[1:-1], 1):
        if abs(val) < abs(zz[i-1]) and abs(val) < abs(zz[i+1]):
            fit = minimize_scalar(lambda x: abs(function(float(x))["ZZ_kHz"]),
                                  bounds=(grid[i-1], grid[i+1]), method="bounded",
                                  options={"xatol": 1e-9})
            candidates.append(function(fit.x))
    return min(candidates, key=lambda row: abs(row["ZZ_kHz"])), rows

for name, reference in sources.items():
    path = root / reference
    matrix = json.loads(path.read_text(encoding="utf-8"))
    model = matrix.get("EcJ", matrix)
    record = dict(design=name, fluxonium_flux_Phi0=.5, fluxonium_EJ_GHz=4.,
                  fluxonium_EL_GHz=.6, source=reference, source_sha256=sha256(path.read_bytes()).hexdigest(),
                  hamiltonian_source=repo_model.relative_to(root).as_posix(),
                  hamiltonian_source_sha256=sha256(repo_model.read_bytes()).hexdigest(),
                  python=sys.version, inputs=dict((k,v) for k,v in model.items() if k.startswith(("Ec_","J_"))),
                  basis_checks=[])
    for nf, nt in [(8, 6), (12, 8), (14, 10)]:
        circuit = Circuit(model, nf, nt)
        if name == "Design2":
            row = circuit.at(16.87)
            record["basis_checks"].append(dict(fx_levels=nf, tc_levels=nt, **row))
            record["reference_point"] = row
            continue
        @lru_cache(maxsize=None)
        def at(phi):
            ej = np.hypot(23.056*np.cos(np.pi*phi), 9.746*np.sin(np.pi*phi))
            return dict(TC_flux_Phi0=phi, **circuit.at(float(ej)))
        best, grid = optimum(at)
        record["basis_checks"].append(dict(fx_levels=nf, tc_levels=nt, **best))
        record["minimum_absolute_ZZ"] = best
        record["sweet_spot"] = grid[0]
        record["half_flux"] = grid[-1]
        record["grid"] = grid
    output.mkdir(parents=True, exist_ok=True)
    (output / f"{name}_static_zz.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
    print(name, record["basis_checks"], flush=True)
