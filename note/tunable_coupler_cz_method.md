# Tunable-coupler CZ: methods

These are general methods for any Design with two flux-tunable transmons coupled through a flux-tunable transmon coupler.
- **Operating procedure and checklist:** [`skill/two_qubit_gate_design.md`](../skill/two_qubit_gate_design.md).
- **Design-specific numbers:** `Design<N>/notes/`.

## 1. Hamiltonian

The full charge-basis model used for all gate numbers is:

    H = Σ_i [4 Ec_i n_i² − EJ_i(Φ_i) cos φ_i] + g_1c n_1 n_c + g_2c n_2 n_c + g_12 n_1 n_2,

where:
- g_ij = 4e² (C⁻¹)_ij / h are the charge-unit couplings;
- C is the Maxwell capacitance matrix of the islands with all other conductors grounded;
- Ec_i = e² (C⁻¹)_ii / 2h.

SQUID junctions:
- Asymmetric SQUID: EJ(Φ) = EJΣ √(cos²(πΦ/Φ0) + d² sin²(πΦ/Φ0)).
- Symmetric coupler SQUID: EJ(Φ) = EJmax |cos(πΦ/Φ0)|.
- EJ = (Φ0/2π)² / L_J.

The model is diagonalised exactly in a truncated product basis (e.g. 4 × 5 × 4 levels, ncut 15). States are labelled by their maximum overlap with the bare states.

Converting to the exchange picture (for intuition and literature comparison only):
- g_ex,ij = g_ij ⟨0|n_i|1⟩⟨0|n_j|1⟩, with ⟨0|n|1⟩ ≈ (EJ/8Ec)^(1/4)/√2.
- The exchange values are larger than the charge-unit values.

## 2. Coupler-mediated coupling and ZZ (leading order, exchange picture)

These formulas are for intuition only; quote numbers from the full model.
- Detunings: Δ_i = f_i − f_c (< 0 with the coupler above the qubits), Δ12 = f1 − f2.
- Effective qubit–qubit coupling (Yan et al., PRApplied 10, 054062, 2018):

      g_eff ≈ g_12 + (g_1c g_2c / 2) (1/Δ1 + 1/Δ2).

  The coupler term is negative when the coupler is above the qubits. It cancels the direct term at one bias, and that point is near the ZZ-off point.
- Static ZZ of two coupled transmons, from the |11⟩–|20⟩ and |11⟩–|02⟩ repulsion:

      ZZ ≈ 2 g_eff² (α1 + α2) / ((Δ12 + α1)(α2 − Δ12)).

  - This two-level-pair formula ignores the direct coupler contributions to ZZ, which matter at g_qc ≈ 50 MHz. Use it only to see where the denominator changes sign (|Δ12| = |α|). Do not use it for the sign or size of ZZ.
  - **Full-model result (Design6):**
    - With |Δ12| > |α|, ZZ kept one sign (negative) over the whole coupler range, so no off point existed.
    - Inside the straddling window |Δ12| < |α|, ZZ was positive and grew strongly as the coupler approached. An off point existed only when g12 cancelled the coupler-mediated term.
  - Always locate the zero, and read the sign, from the full model.
- Sensitivity: near the zero, ZZ is roughly ∝ g_eff², and g_eff is linear in g_12. A 1 MHz change in g12 can therefore move or remove the off point. This is why g12 needs its own EM check.

## 3. Gate models

All gates use a flat-top envelope w(t) with cosine ramps, total time T. Propagation:
- piecewise-constant steps (dt 0.05 ns), each an exact exponential of H(t);
- the gate matrix is evaluated in the dressed computational basis at idle, in the rotating frame, with virtual Z removed.

1. **Coupler-only.** The control is Φ_c(t) = Φ_idle + a_c w(t), with the qubits at their sweet spots.
   - The conditional phase is ≈ 2π ∫ZZ(t) dt, so T_flat ≈ 1/(2 |ZZ_on|).
   - The coupler bound f_c ≥ f_Q1 + margin sets ZZ_on.
2. **Q1 + TC.** Q1 is pulsed to f_Q1 = f_Q2 − α1, the |11⟩–|20⟩ condition, while the coupler pulse raises g_eff.
   - The phase comes from one coherent |11⟩ → |20⟩ → |11⟩ oscillation, which gives the shortest gates.
   - The |11⟩–|02⟩ condition (f_Q1 = f_Q2 + α2) is avoided, because Q1 would cross Q2 during the ramps.

Fidelity: F = (|Tr(U_CZ† U)|² + Tr(U†U)) / 20, coherent only. Leakage = 1 − mean column norm of U on the computational subspace.

## 4. Known behaviour (from the Design6 explorations, 2026-11-17)

- **Coupler-only, |Δ12| > |α|:** with the qubit order Q1 above Q2, ZZ was negative everywhere. The gate needed hundreds of ns: Design6 needed 700 ns at Δ12 = 400 MHz.
- **Straddling with g12 matched:** an off point exists, and 160 ns gates reach F ≈ 0.998–0.9995. However, |11⟩ hybridises strongly with |20⟩ during the flat top, with peak population 0.2–0.5.
  - Leakage is very timing-sensitive, because the ramp-down must return that population.
  - Smoother adiabatic pulses are the next step.
- **Q1 + TC at Δ12 ≈ 400 MHz:** 60–80 ns gates reach F ≈ 0.9992–0.9996, coherent only.

Design-specific values and runs: `Design6/notes/cz_coupler_flux_20261117.md`.
