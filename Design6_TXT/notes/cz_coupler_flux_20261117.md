# Design6_TXT coupler-flux CZ: ZZ off/on points and dynamics (2026-11-17, exploration)

The first sections record historical revision-A parameter exploration. The final section records the r9 Hamiltonian runs.

All parameter sets are built on the revision-A HFSS v8 core fit ([`fit.json`](../results/qubit/2026-11-14/hfss_core_v8/fit.json), `needs_revalidation`):
- **Ec:** from that fit. It is assumed unchanged when g_qc is scaled; a larger C_qc would lower Ec by a few %.
- **Junctions:** re-solved by scqubits for the stated sweet-spot targets. EJ is a user design input, so these are proposals, not adopted values.
- **Couplings:** g in charge units (H term g n_i n_j). This is not the exchange coupling.
- **Errors:** coherent only (no T1/T2), 4 × 5 × 4 levels, dt 0.05 ns.

Entry point: [`simulation/qubit/cz_coupler_flux.py`](../simulation/qubit/cz_coupler_flux.py) (scqenv). It reuses `cz_scqubits.Model / evolve / analyse_U`.

## Gate scheme

- The qubits stay at their sweet spots.
- The coupler is idled at the ZZ-off bias: the ZZ zero, or the minimum |ZZ| when there is no zero, with f_c ≥ f_Q1 + 0.8 GHz.
- It is pulsed towards the qubits (flat top, cosine ramps), bounded at f_c ≥ f_Q1 + 0.3 GHz.
- The conditional phase is ∫ZZ dt, so the gate time is limited by ZZ_on at that bound: T_flat ≳ 1/(2 |ZZ_on|).

## Static result: the off point needs the straddling regime

These points come from a static scratch scan; the dynamics runs below confirm them.
- **Δ12 > |α|** (current targets: Q1 4.80, Q2 4.40 GHz, Δ12 = 400 MHz):
  - ZZ is negative for every coupler bias and never crosses zero.
  - The best "off" point is the coupler sweet spot (−5 kHz), and ZZ_on at the bound is only −0.84 MHz.
- **Straddling Δ12 = 100–150 MHz < |α| ≈ 200 MHz:** ZZ is positive and grows as the coupler approaches. An off point near the coupler sweet spot needs two things:
  - g12 tuned to the coupler-mediated coupling (here ≈ 4 MHz, current layout 2.7 MHz);
  - a high coupler f_max (6.8 GHz gives +5.4 kHz; 6.5 GHz gives +28 kHz).
- The off-point position is very sensitive to g12. A static scan at Q1 4.55 / Q2 4.40 GHz, g_qc 61 MHz, f_c,max 6.3 GHz shows this:
  - g12 = 2.7 MHz: no zero, ZZ at the sweet spot +424 kHz.
  - g12 = 6 MHz: zeros at 5.89 and 5.53 GHz, ZZ at the sweet spot +117 kHz.

## Dynamics runs

| Run | Sets | Ramps / times |
|---|---|---|
| `20261117_151704_cz_coupler_flux_3e09974c` (`results/qubit/20261117_151704_cz_coupler_flux_3e09974c/run.json`) | baseline, straddling f_c,max 6.8 and 6.5 GHz | 10 / 20 ns; baseline 200–700 ns, candidates 60–160 ns |
| `20261117_155016_cz_coupler_flux_5ad31183` (`results/qubit/20261117_155016_cz_coupler_flux_5ad31183/run.json`) | straddling f_c,max 6.8 GHz | 30 / 40 ns; 90–130 ns |

| Set | Q1 / Q2 / TC f_max | g_qc / g12 (MHz) | Idle ZZ | ZZ_on at bound | Best coupler-only CZ |
|---|---|---|---|---|---|
| baseline | 4.80 / 4.40 / 6.30 GHz | 47 / 2.7 | −5.0 kHz | −0.84 MHz | 700 ns, F 0.99996, leakage 3e-5 (ramp 20 ns) |
| strad, f_c,max 6.8 GHz | 4.55 / 4.40 / 6.80 GHz | 61 / 4.0 | +5.4 kHz | +6.5 MHz | **160 ns, F 0.99949, leakage 1.5e-4** (ramp 20 ns) |
| strad, f_c,max 6.5 GHz | 4.55 / 4.40 / 6.50 GHz | 61 / 4.0 | +28 kHz | +6.1 MHz | 160 ns, F 0.99916, leakage 5.4e-4 (ramp 20 ns) |

Findings:
- **Gate time is set by ZZ_on, not by the ramps.** At 60–130 ns the coupler sits at its bound and the conditional phase stops at 200–300° instead of 180° (positive ZZ: the phase counts down from 360°). Longer ramps (30 / 40 ns) do not help. A direct amplitude scan confirmed that 180° is not reached inside the bound.
- **The gate passes through strong |11⟩–|20⟩ hybridization.** In the straddling regime the |20⟩ level is only ≈ 50 MHz from |11⟩. During the flat top, |11⟩ oscillates with |20⟩ (peak population 0.19–0.46), and coupler-excited states reach about 2–3 %. The 160 ns gate works because the ramp-down returns the population. Its leakage is therefore timing-sensitive (1.5e-4 at 160 ns vs 4e-2 at 130 ns). A smoother adiabatic (e.g. Slepian / Martinis–Geller) pulse should be tried before quoting a gate time.
- **Comparison with the existing scheme.** The Q1 + coupler pulse CZ from the 2026-11-14 report (Δ12 400 MHz, current g) reached 80 ns, F 0.99962, coherent only (revision-A, historical). The coupler-only CZ in the straddling regime is about 2× slower at the same F, unless ZZ_on is raised further: larger g_qc, or a coupler bound closer than 0.3 GHz (more coupler excitation).

## Design implications (proposals, not adopted)

A coupler-only CZ needs all of the following:
- **Qubit sweet spots straddling.** Δ12 ≈ 150 MHz requires a Q1 junction change (Q1 f_max 4.80 → 4.55 GHz; EJ is a user decision). It also makes the 7.5 % XY crosstalk an active problem.
- **g_qc ≈ 60 MHz or more** (charge units; current fit 47 MHz). This needs a larger C_qc.
- **g12 ≈ 4 MHz** (current fit 2.7 MHz; the Q3D v8 matrix gives 2.7 MHz from ≈ 0.09 fF). The off point depends on it steeply, so it needs a dedicated Q3D check.
- **Coupler f_max ≈ 6.8 GHz.** Collision with the readout (≈ 6.8 / 7.0 GHz targets) has not been checked.

## r9 Hamiltonian

The r9 calculation uses the adopted pass-20 core Q3D E_C and charge couplings. The layout case uses the drawn L_J of 14.56/29.11 nH (Q1), 17.21/34.41 nH (Q2), and 11.64 nH per TC junction. Both straddling cases use hypothetical junctions to give bare sweet spots of 4.55/4.40/6.80 GHz (Q1/Q2/TC). The runs use 4 × 5 × 4 local levels, charge cutoff 15, dt = 0.05 ns, and no dissipation. [Entry point](../simulation/qubit/cz_r9.py) · r9 CZ manifest (`results/qubit/20261117_181834_cz_r9_20f32abc/run.json`) · Full results (`results/qubit/20261117_181834_cz_r9_20f32abc/cz_summary.json`) · ZZ and fidelity figure (`results/qubit/20261117_181834_cz_r9_20f32abc/cz_r9.png`)

| r9 set | Idle ZZ | Q1 + TC, best sampled | TC only, best sampled |
|---|---:|---|---|
| Drawn L_J, native g_n | −16.80 kHz | 60 ns, F = 0.999254, leakage 1.16 × 10⁻⁴, phase 179.88° | 700 ns, F = 0.999986, leakage 1.18 × 10⁻⁵, phase 180.23° |
| Straddling, native g_n | +17.53 kHz | 80 ns, F = 0.482872, phase 316.52° | 200 ns, F = 0.973681, phase 203.87° at drive bound |
| Straddling, g_qc = 61 MHz and g12 = 4 MHz | ZZ zero at f_TC = 6.6324 GHz | 80 ns, F = 0.569176, phase 295.18° | 160 ns, F = 0.998295, leakage 1.42 × 10⁻³, phase 180.73° |

The dressed idle frequencies with drawn L_J are 4.8070/4.4061/6.4592 GHz (Q1/Q2/TC), and Q1/Q2 anharmonicities are −203.25/−205.86 MHz. The 60 ns point uses 8 ns ramps, f_Q1 = 4.5832 GHz and f_TC = 5.0479 GHz during the gate. The 40, 50 and 80 ns points missed the desired conditional phase in this optimisation. The 700 ns coupler-only point uses 20 ns ramps, with ZZ_on = −1.18 MHz at the simulated drive bound. The tested downward Q1 pulse starts below its target anticrossing in the straddling sets. With the stronger-coupling proposal, peak |20⟩ population in the 160 ns result is 20.6%. Pulse robustness, truncation sensitivity and dissipative fidelity remain unknown.

Grouping repeated flat-top steps in a second computation reproduced the original fidelities to within 5.6 × 10⁻¹². Independent propagation run (`results/qubit/20261117_183534_cz_r9_84b84895/cz_summary.json`)

## QuTiP check of the r9 CZ (2026-11-23)

Independent check of the r9 "drawn L_J, native g_n" Q1 + TC point (60 ns, 8 ns ramps) on the same Hamiltonian: adopted r9 core Q3D pass 20 and the drawn L_J in `layout/design6_layout.py`. No EM rerun. Later mask revisions (r10 tails, r12/r13 SQUID ground tail and 4 µm strips, 40 µm cheese) did not rerun Q3D, so r9 pass 20 is still the newest Hamiltonian source. [Entry point](../simulation/qubit/cz_qutip_r9.py) · run `results/qubit/20261123_165526_cz_qutip_r9_636ae45f/` (`cz_qutip_summary.json`, `cz_qutip_r9.png`).

| Check | Result |
|---|---|
| QuTiP `sesolve` (adaptive ODE) vs parent piecewise propagator, 4 × 5 × 4 | F = 0.999254 vs 0.999254, phase 179.88°, leakage 1.16 × 10⁻⁴, max \|ΔU\| 7.5 × 10⁻⁵ |
| Same pulse, 5 × 6 × 5 (ncut 15) / 6 × 7 × 6 (ncut 20) | **F = 0.642 / 0.641**, phase 278°, leakage 4.9 × 10⁻² |
| Pulse re-optimised at 6 × 7 × 6, ncut 20 (60 ns, 8 ns ramps) | **F = 0.999559**, phase 180.10°, leakage 2.5 × 10⁻⁴, a_q 0.13538 (parent 0.13973), a_c 0.28781 |
| Gate-time error ±0.5 / ±1 / ±2 ns (4 × 5 × 4) | F ≥ 0.99892 / 0.99822 / 0.99531 |
| Q1 flux-amplitude error ±0.1 / ±0.3 / ±1 % | F ≈ 0.99878 / 0.99493 / 0.953 |
| Coupler flux-amplitude error ±0.1 / ±0.3 / ±1 % | F ≥ 0.99914 / 0.99817 / 0.98895 |

Decoherence (QuTiP `mesolve`, constant Lindblad T1 and pure dephasing on Q1, TC and Q2; T1 = T2 assumed, not measured; flux-noise dephasing away from the sweet spot not modelled). The added error comes from a reduced 3 × 4 × 3 basis, where the parent pulse is not a clean CZ (F = 0.9246 without loss), so it is an estimate of the time-integrated decoherence of a 60 ns gate.

| T1 = T2 (µs) | 200 | 100 | 50 | 30 | 20 |
|---|---:|---:|---:|---:|---:|
| Decoherence-added error | 3.5 × 10⁻⁴ | 7.1 × 10⁻⁴ | 1.4 × 10⁻³ | 2.4 × 10⁻³ | 3.5 × 10⁻³ |
| F with the re-optimised coherent 0.99956 | 0.99921 | 0.99885 | 0.99815 | 0.99721 | 0.99604 |

Findings:
- The CZ exists on the r9 Hamiltonian and reaches F ≈ 0.9996 coherent at 60 ns, but the parent run's 4 × 5 × 4 pulse parameters are not converged in level count. Use ≥ 5 levels per qubit (≥ 6 × 7 × 6 checked) for pulse calibration in simulation; the parent's 0.999254 is a model number for that basis only.
- The gate is a sharp |11⟩–|20⟩ return: Q1 flux must be held to ~0.1 % (about 1.4 × 10⁻⁴ Φ0) and gate time to ~0.5 ns for < 5 × 10⁻⁴ extra error. In the lab this is calibrated, not designed.
- With T1 = T2 ≈ 50–100 µs, decoherence dominates and F ≈ 0.998–0.999 is the expected range. Dissipative numbers are on the reduced basis with assumed coherence times; validation status: unknown.
