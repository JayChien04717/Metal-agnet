"""Calculate transformed LOM capacitance parameters for the FTF matrix."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.constants import e, h, hbar, pi
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
MATRIX_CSV = WORKSPACE / "FTF_capacitance_matrix.csv"
PARAMS_CSV = WORKSPACE / "FTF_qubit_parameters_from_LOM.csv"
CHARGE_COUPLINGS_CSV = WORKSPACE / "FTF_charge_line_couplings_from_LOM.csv"
CTILDE_CSV = WORKSPACE / "FTF_Ctilde_fF.csv"
CTILDE_INV_CSV = WORKSPACE / "FTF_Ctilde_inverse_1_per_F.csv"


ORDER = [
    "pad_left_FL",
    "pad_right_FL",
    "cross_TC",
    "pad_left_FR",
    "pad_right_FR",
    "readout_launch_FL",
    "prime_cpw_LR_tee",
    "readout_launch_FR",
]

NEW_LABELS = [
    "FL_sum",
    "FL_diff",
    "TC",
    "FR_sum",
    "FR_diff",
    "res_FL",
    "res_c",
    "res_FR",
]

CHARGE_LINES = [
    "charge_launch_FL",
    "charge_line_TC",
    "charge_launch_FR",
]

M = np.array(
    [
        [1, 1, 0, 0, 0, 0, 0, 0],
        [1, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
    ],
    dtype=float,
)


def main() -> None:
    if not MATRIX_CSV.exists():
        raise FileNotFoundError(
            f"Missing {MATRIX_CSV}. Run scripts/run_ftf_q3d_lom_debug.py "
            "or the FTF.ipynb Q3D cells first to generate the raw capacitance matrix."
        )

    cap_ff = pd.read_csv(MATRIX_CSV, index_col=0)
    cap_f = cap_ff.loc[ORDER, ORDER].to_numpy(dtype=float) * 1e-15

    minv = np.linalg.inv(M)
    ctilde = minv @ cap_f @ minv
    ctilde_inv = np.linalg.inv(ctilde)

    ctilde_ff = pd.DataFrame(ctilde / 1e-15, index=NEW_LABELS, columns=NEW_LABELS)
    ctilde_inv_df = pd.DataFrame(ctilde_inv, index=NEW_LABELS, columns=NEW_LABELS)

    def energy_over_h_ghz(coef: float, i: int, j: int) -> float:
        return coef * e**2 * ctilde_inv[i, j] / h / 1e9

    params = [
        ("EC1/h", "GHz", energy_over_h_ghz(0.5, 1, 1)),
        ("ECc/h", "GHz", energy_over_h_ghz(0.5, 2, 2)),
        ("EC2/h", "GHz", energy_over_h_ghz(0.5, 4, 4)),
        ("J1c/h", "GHz", energy_over_h_ghz(4.0, 1, 2)),
        ("J2c/h", "GHz", energy_over_h_ghz(4.0, 2, 4)),
        ("J12/h", "GHz", energy_over_h_ghz(4.0, 1, 4)),
    ]

    z0 = 50.0
    prefactor = e * np.sqrt(pi * hbar / (2 * z0))

    def fluxonium_g_over_2pi_mhz(i: int, j: int) -> float:
        hbar_g = prefactor * ctilde_inv[i, j]
        return hbar_g / h / 1e6

    for qi, qname in [(1, "FL_diff"), (4, "FR_diff"), (2, "TC")]:
        for ri, rname in [(5, "res_FL"), (6, "res_c"), (7, "res_FR")]:
            params.append(
                (
                    f"g({qname},{rname})/2pi",
                    "MHz",
                    fluxonium_g_over_2pi_mhz(qi, ri),
                )
            )

    cxy = cap_ff.loc[CHARGE_LINES, ORDER].to_numpy(dtype=float) * 1e-15
    cxy_tilde = cxy @ minv
    charge_beta = cxy_tilde @ ctilde_inv
    charge_beta_df = pd.DataFrame(
        charge_beta,
        index=CHARGE_LINES,
        columns=NEW_LABELS,
    )

    charge_rows = []
    for line in CHARGE_LINES:
        for qname in ["FL_diff", "TC", "FR_diff"]:
            beta = charge_beta_df.loc[line, qname]
            # Hamiltonian term is 2 e V beta n. This is the frequency scale for
            # V = 1 mV, reported as (2 e beta 1mV) / h.
            per_mv_ghz = 2 * e * beta * 1e-3 / h / 1e9
            charge_rows.append((line, qname, beta, per_mv_ghz))

    pd.DataFrame(
        charge_rows,
        columns=[
            "charge_line",
            "coordinate",
            "beta_Cxy_Cinv",
            "drive_per_mV_over_h_GHz",
        ],
    ).to_csv(CHARGE_COUPLINGS_CSV, index=False)

    pd.DataFrame(params, columns=["parameter", "unit", "value"]).to_csv(
        PARAMS_CSV, index=False
    )
    ctilde_ff.to_csv(CTILDE_CSV)
    ctilde_inv_df.to_csv(CTILDE_INV_CSV)

    print("old order:", ORDER)
    print("new order:", NEW_LABELS)
    print("\nCtilde (fF):")
    print(ctilde_ff.round(6).to_string())
    print("\nParameters:")
    for name, unit, value in params:
        print(f"{name:28s} {value: .8g} {unit}")
    print("\nCharge-line drive couplings:")
    for line, qname, beta, per_mv_ghz in charge_rows:
        print(
            f"{line:18s} -> {qname:8s} "
            f"beta={beta: .8g}, (2e beta * 1mV)/h={per_mv_ghz: .8g} GHz"
        )
    print(f"\nsaved {PARAMS_CSV}")
    print(f"saved {CHARGE_COUPLINGS_CSV}")
    print(f"saved {CTILDE_CSV}")
    print(f"saved {CTILDE_INV_CSV}")


if __name__ == "__main__":
    main()
