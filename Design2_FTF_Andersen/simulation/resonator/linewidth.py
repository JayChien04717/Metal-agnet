"""Convert Design2_FTF_Andersen Maxwell coupling capacitances to external readout linewidths."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import argparse
import csv
import math
from scipy.special import ellipk

def cpw_impedance(width_um=15.0, gap_um=10.0, silicon_eps=11.9):
    """Thick-substrate, thin-metal quasi-static CPW impedance."""
    modulus = width_um / (width_um + 2 * gap_um)
    return 30 * math.pi * ellipk(1-modulus**2) / (
        math.sqrt((silicon_eps+1)/2) * ellipk(modulus**2))

def linewidth(f_ghz, coupling_ff, z_line, z_res):
    omega = 2 * math.pi * f_ghz * 1e9
    c = coupling_ff * 1e-15
    # One side-coupling capacitor sees two matched directions of the through line.
    r = z_line/2
    g = (omega*c)**2*r / (1+(omega*c*r)**2)
    c_mode = math.pi / (2*omega*z_res)
    q_ext = omega*c_mode/g
    return q_ext, f_ghz*1000/q_ext

def read_matrix(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(row for row in csv.DictReader(handle) if row[""] == "metal_feed")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left_GHz", type=float)
    parser.add_argument("right_GHz", type=float)
    parser.add_argument("--matrix", type=Path, default=_P('Design2_FTF_Andersen/results/general/legacy/layout3_full/capacitance_fF.csv'))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--basis", default="isolated readout eigenmode")
    args = parser.parse_args()
    matrix = read_matrix(args.matrix)
    z = cpw_impedance()
    rows = []
    for side,f in (("left",args.left_GHz),("right",args.right_GHz)):
        c = -float(matrix[f"metal_res_{side}"])
        q,kappa = linewidth(f,c,z,z)
        rows.append(dict(side=side,frequency_basis=args.basis,frequency_GHz=f,Cc_fF=c,
                         Z_feed_ohm=z,Z_res_ohm=z,Q_ext=q,kappa_ext_over_2pi_MHz=kappa,
                         kappa_previous_45p1_ohm_MHz=linewidth(f,c,45.1,45.1)[1]))
    if args.output:
        with args.output.open("w", newline="", encoding="utf8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    for row in rows:
        print(row)

if __name__ == "__main__":
    main()
