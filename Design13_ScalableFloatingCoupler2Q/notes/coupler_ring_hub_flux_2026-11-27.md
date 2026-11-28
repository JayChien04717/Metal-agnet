# Ec, readout coupling and flux-line crosstalk tuning → ring-hub coupler (r7)

Date: 2026-11-27. Stage 2 checks requested by the user (small models; not full-chip).

Scripts:
- Q3D (capacitance, `g`, `χ`): [`simulation/qubit/q3d_ro_coupler_crosstalk.py`](../simulation/qubit/q3d_ro_coupler_crosstalk.py). Q0 + 4 couplers + 4 neighbour qubits; 8 passes, 0.5 %; single runs (convergence not separately verified).
- HFSS (flux-line `M`, crosstalk): `simulation/flux/hfss_mutual_arc.py` (now [`hfss_mutual_flux.py`](../simulation/flux/hfss_mutual_flux.py)). Q0 + `TC_E` + `TC_N` (SQUIDs 1.41 mm apart). DrivenModal Z-parameters, 12 passes, ΔS target 2e-4, discrete 10–100 MHz. The JJ1 ports are port-only (open-loop EMF); JJ2 carries the lumped L.
- `simulation/general/rotated_jj_render.py` (removed in r9, snapshot kept in the run folders): qiskit-metal's HFSS renderer draws junction sheets axis-aligned only. Design13's −45° JJs need rotated sheets, otherwise ports and RLC sheets are wrong and Analyze fails. Install it in every HFSS run of this design.

## 1. Ec, g_qc and readout coupling (Q3D)

Targets: `Ec(Q) 190–200 MHz`, `Ec(TC) 180–200 MHz`, `g_qc 70–85 MHz`, `|2χ| ≈ 1.5 MHz`, `n_crit > 60` (`f_q 4.65 GHz`, `f_r 7.17 GHz`, `f_c 6.2 GHz`; `χ = g² α/(Δ(Δ+α))` with `α ≈ −Ec`). `g` values are pad-level estimates.

| Tag (`results/qubit/…`) | Qubit / readout | Coupler hub | Ec Q / TC MHz | g_qc | g_qr | 2χ | n_crit | g_cr/g_qr |
|---|---|---|---|---|---|---|---|---|
| r6 as drawn (`…c28_ro26`) | Ø256, arc 26° gap 34 | split disc r42 | 219 / 265 | 94 | 58 | 0.21 | — | 8.8 % |
| `ec_v1` | Ø280, arc gap 10 | split disc r70 | 195 / 236 | 88 | 120 | 0.81 | 111 | 4.5 % |
| `ec_r2a` | Ø280, arc gap 5 w28 | split disc r120 | 195 / 198 | 79 | 154 | 1.35 | 67 | 3.6 % |
| `key_a14g20` | Ø280, keyhole r14 g20 | split disc r120 | 194 / 198 | 80 | 163 | 1.50 | 60 | 1.7 % |
| `ring_d60g6` | Ø280, keyhole r14 g20 | ring disc 60 / gap 6 / w16 | 192 / 188 | 78 | 164 | 1.51 | 59 | 1.6 % |
| **`ring_v3` (adopted r7)** | **Ø280, keyhole r13 g20** | **ring disc 60 / gap 6 / w10, hub gap 10, SQUID 24×70** | **192 / 197** | **81** | **159** | **1.41** | **63** | **1.7 %** |

- **Keyhole readout** (user suggestion): the island wraps the readout pad, so C_qr is about 10 fF in a small footprint, and the island shields the pad from the coupler arcs. Parasitic drops from 3.6–4.3 % (arc) to 1.6–1.7 %. The trade-off is a denser field in the 20 µm keyhole gap than in the 46 µm pocket gap (TLS participation not evaluated).

## 2. Flux line: own M and coupler-to-coupler crosstalk (HFSS)

Targets: own `M ≈ 0.7–1.0 pH`, crosstalk `M(neighbour)/M(own) ≤ 1 %`. Every completed run passed reciprocity (`< 1e-9`), `Re Z ≈ 0` and plateau-flatness (`< 1 %`) checks.

| Tag (`results/flux/…`) | Hub / SQUID | Flux line | Own M pH | I_Φ0 mA | Crosstalk TC_N←fl_E / TC_E←fl_N | ΔS (12 passes) |
|---|---|---|---|---|---|---|
| `hub120_arcL` | split disc r120, SQUID at centre | arc-L | 0.72 / 0.69 | 2.9 | **21.8 % / 27.0 %** | 5.1e-5 |
| `hub120_T` | split disc r120, SQUID at centre | rectangular T | 0.011 / 0.054 | 185 / 39 | (ratio of ~0 values) | 2.3e-4 |
| `ring_d60g6_arcT` | ring, rim SQUID 20×40 | arc-T (25°/50°) | 0.25 / 0.25 | 8.2 | 0.68 % / 3.98 % | 5.0e-5 |
| `ring_v2_arcT` | ring w10, hub gap 10, SQUID 20×64 | arc-T (40°/70°) | 0.61 / 0.64 | 3.3 | 0.44 % / 1.58 % | **4.6e-4 (not converged)** |
| **`ring_v3_arcT` (adopted r7)** | **ring w10, hub gap 10, SQUID 24×70** | **arc-T (40°/70°)** | **0.59 / 0.62** | **3.5 / 3.3** | **0.57 % / 1.65 %** | 1.9e-4 |

Interpretation:
- **arc-L with the SQUID at the centre of a large hub:** the return current circulates around the hub pocket, which makes a large magnetic dipole. Own M is fine, but crosstalk to the SQUID 1.41 mm away is 22–27 %, so it is unusable.
- **Rectangular T with the SQUID at the centre:** the symmetric arms cancel at a loop that is 100+ µm from the pocket edge, so own M ≈ 0.
- **Ring hub (inner disc + outer crescent) with the SQUID on the hub rim and an arc-T flux line:** quadrupole-like far field. Crosstalk drops about 10–40× (to ≤ 1.7 %), and own M is ~0.6 pH (15 % below the 0.7 pH lower bound; I_Φ0/2 ≈ 1.7 mA).
- **Enlarging the loop further (v2 → v3) did not raise M.** With the crescent inner radius fixed, a wider JJ spacing moves the flat SQUID gap inward, away from the flux line.

## 3. Status and limitations

- **Adopted for layout r7:** qubit Ø280 µm with keyhole readout (pad r 13 µm, gap 20 µm); ring hub (disc 60 µm, ring gap 6 µm, crescent 10 µm, hub gap 10 µm, SQUID gap 24 µm × spacing 70 µm, JJs at −45°); arc-T flux line (offset 40°, arms 70°).
- **Open:**
  - own M is ~15 % below target;
  - crosstalk in one direction is 1.65 % (> 1 %);
  - 2χ is 1.41 MHz (6 % below target);
  - all values are small-model, single-run numbers; the 12-vs-15-pass convergence check and full-chip return paths were not run;
  - L_loop reads ~0.66 × the lumped L, probably because the junction sheet overlaps lead metal. This does not affect the open-loop M.
