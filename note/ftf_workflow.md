# FTF architecture: fluxonium – tunable transmon coupler – fluxonium

General architectural notes for planar FTF chips (two fluxonium qubits coupled through a central tunable transmon coupler). General coding, routing, GDS export, and figure rules live in `../skill/` ([skill.md](../skill/skill.md), [gds_export.md](../skill/gds_export.md), [layout_figures.md](../skill/layout_figures.md)). Chip-specific dimensions, entry scripts, and port maps belong in each Design's `README.md` and `notes/`.

---

## 1. Layout parametrization

* **Derived coupler geometry:** derive coupler `pad_arm_length` and height from the fluxonium pad dimensions and target paddle-to-pad gap (`PADDLE_SEP`) so resizing the fluxonium pads moves the coupler automatically.
* **Per-qubit overrides:** separate shared qubit geometry (`FX_OPTIONS`) from per-qubit overrides (`FX_PER_QUBIT`, e.g., readout pad gap per resonator frequency).
* **Resonator sizing inputs:** compute tee coupling length and target filleted length from per-resonator `f`, `kappa`, Q3D `C_pad`, any residual `C_extra`, and `hfss_scale` (`f_HFSS / f_target`).
* **Standard 8-port assignment:** feedline on top/bottom ports (P7/P8), fluxonium flux lines on P1/P6, fluxonium charge lines on P2/P5, coupler flux on P3, coupler charge on P4.

---

## 2. Simulation model vs. junction-free mask

Use **one fluxonium simulation model** across Q3D, HFSS eigenmode, and mutual-inductance runs: `FluxoniumArrayJJDolan` (`pattern design/fluxonium/fluxonium_nthu_dolan.py`).
* Pads, stems, pocket, and control lines match the mask (`MyFluxoniumDolan`).
* The loop (`MyFluxoniumDolan._loop`) consists of a short centre JJ gap (`L_j`, or port 2 in mutual-inductance runs) and a U-shaped array of two junction-table sides (`array_top` / `array_bot`, `l_inductance / 2` each) joined by a metal closing end (`array_end`). An HFSS lumped sheet must end on metal; in Q3D `array_end` is a floating net eliminated via Kron reduction by `code/qubit/cap_matrix.py`. Do not substitute a single straight inductor.
* Shorted flux lines are rendered as ground (`flux_as_ground='True'`) in Q3D to avoid net-overlap errors.

### What `jj_metal='False'` removes on the base-metal mask

| Shared component (`pattern design/`) | Dropped when `jj_metal='False'` | Kept on base-metal mask |
|---|---|---|
| `MyFluxoniumDolan` (`fluxonium/fluxonium_nthu_dolan.py`) | 2 µm lead line, metal U closing end (`array_end`) | `pad_top`, `pad_bot` (with stems), pocket, flux U-loop, readout pad, charge line |
| `MyTmonCouplerDolan` (`transmon/my_tmon_dolan.py`) | SQUID block (island electrode, loop hole, `squid_ground`), unless landing pads are explicitly retained | bar, paddles, stem, etch window |
| `FluxoniumDelftDolan` (`fluxonium/fluxonium_delft_dolan.py`) | 2 µm lead line, metal hairpin end (`array_end`) | round pads with stems, pocket, flux U-loop, readout pad, charge line |
| `Tmon_c` (`transmon/tmon_c.py`) | `JJpad_bottom` / `JJpad_up` | stem, arms, palms, etch |
| `HookedTestJunction` (`junction_test_coupon.py`) | junction (never drawn) | etched pocket + two floating pads |

For the full fab export sequence (`jj_metal=False`, `fab_contacts=True`, optional drawn e-beam junctions, corner F/T coupons, cheesing, and verification), follow [skill/gds_export.md §2.1](../skill/gds_export.md#21-fab-gds-workflow).
