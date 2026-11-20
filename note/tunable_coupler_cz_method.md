# Tunable-coupler CZ: physics and methods

General physical model and gate-design methods for two flux-tunable transmons (Q1, Q2) coupled through a flux-tunable transmon coupler (TC).

* **Operating procedure & checklist:** [`skill/two_qubit_gate_design.md`](../skill/two_qubit_gate_design.md).
* **Design-specific numbers and runs:** each Design's `notes/` and `design.json`.

---

## 1. Hamiltonian

All gate numbers come from exact diagonalisation of the full charge-basis Hamiltonian:

$$
H = \sum_i \left[ 4 E_{c,i} n_i^2 - E_{J,i}(\Phi_i)\cos\varphi_i \right] + g_{1c} n_1 n_c + g_{2c} n_2 n_c + g_{12} n_1 n_2,
$$

where:
- `g_ij = 4 e² (C⁻¹)_ij / h` (charge-unit couplings) and `E_{c,i} = e² (C⁻¹)_ii / (2 h)` from the Maxwell capacitance matrix `C` of the islands with all other conductors grounded;
- Asymmetric SQUID: `E_J(Φ) = E_{JΣ} sqrt(cos²(π Φ/Φ0) + d² sin²(π Φ/Φ0))`, with `d = (E_{Ja} − E_{Jb}) / (E_{Ja} + E_{Jb})`;
- Symmetric coupler SQUID: `E_J(Φ) = E_{J,max} |cos(π Φ/Φ0)|`, with `E_J = (Φ0 / 2π)² / L_J`.

Diagonalise in a truncated product basis (e.g. `4 × 5 × 4` levels, `ncut = 15`), labelling dressed states by maximum overlap with bare states.

**Conversion to exchange picture** (for intuition and literature comparison only):
`g_{ex,ij} = g_{ij} ⟨0|n_i|1⟩⟨0|n_j|1⟩`, with `⟨0|n|1⟩ ≈ (E_J / 8 E_c)^{1/4} / √2` (~1.2–1.5 for `E_J/E_c ≈ 50–150`, making `g_ex` ~1.5–2× larger than charge-unit `g`).

---

## 2. Coupler-mediated coupling and ZZ (leading-order intuition)

Use these perturbative expressions for intuition only; always extract ZZ sign and zeros from the full model:
- Detunings: `Δ_i = f_i − f_c` (`< 0` with coupler above qubits), `Δ12 = f1 − f2`.
- Effective qubit–qubit coupling (Yan et al., *Phys. Rev. Applied* 10, 054062, 2018):

  $$
  g_{\mathrm{eff}} \approx g_{12} + \frac{g_{1c} g_{2c}}{2} \left(\frac{1}{\Delta_1} + \frac{1}{\Delta_2}\right).
  $$

  With the coupler above the qubits, the coupler term is negative and cancels direct `g_12` near the ZZ-off point.
- Two-level-pair static ZZ estimate:

  $$
  ZZ \approx \frac{2 g_{\mathrm{eff}}^2 (\alpha_1 + \alpha_2)}{(\Delta_{12} + \alpha_1)(\alpha_2 - \Delta_{12})}.
  $$

  This omits direct coupler contributions to ZZ (significant at `g_qc ~ 50 MHz`). Use it only to identify the straddling boundary `|Δ12| = |α|`.
- **Sensitivity:** near the ZZ zero, `ZZ ∝ g_eff²` and `g_eff` is linear in `g_12`, so a ~1 MHz shift in `g_12` can move or eliminate the off point—always verify `g_12` with Q3D.

---

## 3. Gate models and regime behaviour

Time-domain evolution uses a flat-top envelope `w(t)` with cosine ramps (total duration `T`), propagated via piecewise-constant exact matrix exponentials (`dt ≈ 0.05 ns`), evaluated in the idle dressed computational basis in the rotating frame with virtual Z removed:
- **Coherent average gate fidelity:** `F = (|Tr(U_CZ† U)|² + Tr(U† U)) / 20`.
- **Leakage:** `1 − mean column norm` of `U` on the computational subspace.

### General regime behaviour

1. **Coupler-only CZ (`Φ_c(t) = Φ_idle + a_c w(t)`, qubits at sweet spots):**
   - Outside straddling (`|Δ12| > |α|`): ZZ keeps one sign over the full coupler range, so no exact ZZ-off point exists and gates are slow (hundreds of ns).
   - Inside straddling (`|Δ12| < |α|`) with `g_12` matched: a ZZ-off point exists and gates reach ~100–200 ns, constrained by `T_flat ≳ 1 / (2 |ZZ_on|)` at the coupler proximity bound (`f_c ≥ f_Q + margin`). Strong `|11⟩–|20⟩` hybridisation during the flat top makes leakage timing-sensitive unless adiabatic pulse shaping is used.
2. **Q1 + TC CZ (Q1 pulsed to `f_Q1 = f_Q2 − α1` for `|11⟩–|20⟩` plus TC pulse):**
   - Works for `|Δ12| > |α|` and achieves the shortest gates (~60–80 ns coherent) via a single `|11⟩ → |20⟩ → |11⟩` oscillation.
   - Avoid `|11⟩–|02⟩` (`f_Q1 = f_Q2 + α2`) when Q1 would have to cross Q2 during the ramp (`|01⟩–|10⟩` swap error).
