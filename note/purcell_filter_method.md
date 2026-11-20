# Readout channel with an individual Purcell filter: methods and layout

General circuit and electromagnetic methods for any Heinsoo / Wallraff / Cai 2024 style readout channel (`feed line – C_p – λ/4 Purcell filter – C_j hammer – λ/4 readout resonator`).

* **Operating procedure & checklist:** [`skill/purcell_resonator_design.md`](../skill/purcell_resonator_design.md).
* **Standard component & circuit model:** `pattern design/resonator/purcell_readout.py` (`PurcellReadout`) and [`code/resonator/purcell_readout_model.py`](../code/resonator/purcell_readout_model.py).
* **Design-specific numbers and validation runs:** live in each Design's `README.md` and `notes/`.
* **References:** Wallraff group PhD thesis (ETH Zurich), sec. 2.5 "Qubit Readout", pp. 119–127 (Fig. 2.29, eqs. 2.119–2.137, Table 2.4); Cai et al., *Chin. Phys. B* 33, 068501 (2024).

---

## 1. First estimates of couplers and lengths

1. **Coupling capacitor for a target Q_c:** `Q_c = π / (Z0 Zr (Cc ω)²)` ⇒ `Cc = sqrt(π / (Q_c Z0 Zr)) / ω`. A two-port transmission filter loaded equally by both ports has `Q_c = 2 Q_F`.
2. **Interdigital capacitance (Alley / Bahl):** for `N` fingers, overlap `l` (µm), finger width `W`, substrate height `h`:
   `C [pF] = (ε_r + 1) l [(N − 3) A1 + A2]`, with `A1 = 4.409e-6 tanh(0.55 (h/W)^0.45)`, `A2 = 9.92e-6 tanh(0.52 (h/W)^0.5)`. Replace with Q3D once drawn.
3. **Loaded λ/2 line:** `f = v / (2 L) (1 − 2 Cc / (C' L))` when both ends are loaded; iterate meanders on filleted length.
4. **Q3D capacitor scaling:** two geometry points fix the local trend (`hammer C_j ≈ a / gap + b`; `enclosing cap C` ≈ linear in plate width).

### Quarter-wave frequency and short-tail trim

1. **Unloaded estimate:** for a uniform λ/4 line (one open, one short), `f_{λ/4} = v / (4 L)`, where `L` is the filleted centre-line length. For a small open-end shunt capacitance `C_end`, `L_eff ≈ L + C_end / C'` with `C' = 1 / (Z0 v)`.
2. **Branched circuit frequency:** split each line at its tee into open-side arm, shorted arm, and hammer arm at their drawn reference planes. With `e^{+jωt}`, `Y0 = 1 / Z0`, `θ = 2π f ℓ / v`:

   $$
   Y_{\mathrm{open}} = j Y_0 \tan\theta, \qquad
   Y_{\mathrm{short}} = -j Y_0 \cot\theta, \qquad
   Y_{\mathrm{in}} = Y_0 \frac{Y_L + j Y_0 \tan\theta}{Y_0 + j Y_L \tan\theta}.
   $$

   Transform each termination admittance along its arm and sum at the tee, inserting Q3D-derived `C_p` and `C_j` π networks (de-embedding CPW capacitance already represented by line sections). Bare frequencies are the zeros of the branch driving-point admittance (real part of complex root when feed-loaded) with off-diagonal transfer suppressed; restore coupling to solve for hybrid poles.
3. **Reference for trimming:** anchor trims to an HFSS simulation with matching topology, termination loads, and fixture. Recover bare frequencies from resolved complex S21 poles (or loaded eigenmode fields); check adaptive mesh stability before converting a model–HFSS offset into a physical trim.
4. **Local trim slope:** hold the tee, hammer, and termination fixed; perturb the shorted arm by `±h` in the circuit to compute `s = ∂f / ∂L_short ≈ [f(L_short + h) − f(L_short − h)] / (2h)`:

   $$
   \Delta L_{\mathrm{short}} \approx \frac{f_{\mathrm{target}} - f_{\mathrm{ref}}}{s}.
   $$

   Record `s` in MHz/µm (`s ≈ −f_ref / L_eff` as a first estimate; negative trim shortens the arm and raises frequency).

---

## 2. Wallraff-group circuit model (no curve fit)

* **Topology (thesis Fig. 2.29, eqs. 2.122–2.136):**
  * Resonator: qubit-side load (`C_qr` in series with `C_qg`, plus `C_qrg` and air-bridge `C_AB`) → line `L_rr` → hammer T → shorted arm (`L_rrg`, `L_rre`).
  * Hammer: resonator node (`C_prg`) – series `C_pr` – filter node (`C_prg`).
  * Filter: feed T → `C_pff` – series `C_pf` (hugger) – `C_pfg` → line `L_pf` → shunt hammer branch → shorted arm (`L_pfg`, `L_pfe`), giving branch admittance `Y_rc`. Put the hugger frame on the **filter** side so its larger ground capacitance loads the filter rather than the feedline.
  * Single-channel matched-feed transmission: `S21 = 2 / (2 + Z0 Y_rc)`.
* **Exact 2-mode extraction (eq. 2.137, zero internal loss):**
  `ω_r` is the real zero of `Y_rc` (peak between the two dips); `l1, l2` are the complex poles of `S21` (zeros of `2 + Z0 Y_rc`), equal to eigenvalues of `[[ω_p + i κ_p/2, J], [J, ω_r]]`:
  - `ω_p = Re(l1 + l2) − ω_r`, `κ_p = 2 |Im(l1 + l2)|`, `J² = (ω_p + i Im(l1 + l2)) ω_r − l1 l2` (`Im J² ≈ 0` is the consistency check).
  - Sign convention: `e^{+jωt}` gives poles with `Im > 0`.
  - Unknown phase velocity `v`, air-bridge capacitance `C_AB`, and lumped hugger parasitics leave bare `f_p`/`f_r` uncertain by tens of MHz, and single-channel models omit full-feed `C_in` standing waves; take final `f_p`, `f_r`, and `κ_p` from EM / full-feed models.
* **Closed-form J sanity check** ([`code/general/purcell_j.py`](../code/general/purcell_j.py), `purcell_readout_model.j_hammer`):
  `J = (f / 2) · C_j · cos θ_p · cos θ_r / sqrt(C_eff,p · C_eff,r)`, with open-end loads mapped to `x = C_end / C'` and `C_eff = C' (L + x) / 2`.
  Use this as the first check of any circuit/EM model: if EM `J` exceeds this closed form by much more than ~30 %, a parasitic coupling path outside `C_j` exists (Section 6).

---

## 3. Fitting HFSS S21 with a circuit model

1. **ABCD network:** `P1 – TL – C_in(π) – TL – n0 – TL – P2` with shunt branch at `n0` (`stub – C_p(π) – filter – tap – C_j(π) – resonator`), where each π network uses Q3D series C and electrode-to-ground capacitances.
2. **Residual & optimiser:** fit complex `S21` and weighted complex `S11` in linear scale (not dB, so lossless −60 dB nulls do not dominate), with soft Gaussian priors on Q3D capacitance and velocity scales, using global search (`differential_evolution`) followed by `least_squares`.
3. **Practical points:**
   - When fitting eq. 2.137 (`e^{-iωt}`) directly to `e^{+jωt}` complex S21 (HFSS/ABCD), fit `conj(S21)` so `κ` stays positive.
   - Parameterize by the two S21 zeros `z1, z2` (`ω_p = z1 + z2 − ω_r`, `J² = ω_p ω_r − z1 z2`) to lock sharp nulls.
   - `J` is well-constrained only when both dips are resolved near `f_p ≈ f_r` (far detuned `Δ ≫ J` gives only `J²/Δ` shift).
   - In whole-channel Q3D, render shorted ends open so nets do not merge with ground, and use only mutual capacitances.

---

## 4. Two-mode inversion of HFSS eigenmodes

From two complex eigenfrequencies `z1, z2` near resonance, solve `[[f_p + i κ_p/2, J], [J, f_r]]` per point (least squares on trace and determinant). Identify modes by perturbing one line length by ~20 µm.

---

## 5. General EM lessons for Purcell channels

* **Ground bonds everywhere:** unbonded etched stubs + hammer pockets sever the ground plane; return-current detours create parasitic coupling invisible to Q3D (inflating `J` up to ~3×). Auto-bonds skip segments shorter than `wb_threshold`, so add explicit PEC straps on short stubs and tees.
* **Side-by-side CPW runs couple directly:** long parallel filter–resonator runs add capacitive and inductive coupling not represented by `C_j` (Section 6).
* **Eigenmode vs. driven `κ_p`:** eigenmode `κ_p` typically runs ~1.5–1.8× higher than driven S21 / Q3D; trust driven S21 and Q3D for `C_p` and `κ_p`.
* **Far-detuned resonators:** invisible in driven S21 due to narrow linewidth; locate via eigenmode or tune near `f_p = f_r`.

---

## 6. PurcellReadout layout: hammer-only coupling

[`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py) enforces **frequency from length** and **coupling from the hammer `C_j` only**, eliminating the three parasitic coupling mechanisms seen in legacy side-by-side layouts:
1. Long parallel CPW runs between filter and resonator (both capacitive `k_C` and inductive `k_L` coupling).
2. Ground-plane cuts across hammer stubs and pockets (return-current coupling).
3. Neighbouring meander coupling.

### 6.1 Layout principles in the module

1. **Opposite-side routing:** each main line runs straight through its T; the filter's open arm goes to `C_p` while its shorted arm jogs and meanders away from the resonator (and vice versa). The only parallel run is the short T region `2 × (tee_half + lead)` (300 µm default).
2. **Spacing & meanders:** `line_spacing`, `jog`, and `meander_pitch` default to 400 µm (not hard limits). Verify any remaining parallel run with `purcell_readout_model.parallel_j`.
3. **Six ground straps per module (`strap_points()`):** one across each hammer stub and one across each main line on both sides of each T (PEC straps in EM, air bridges on mask).
4. **Perpendicular transitions near the hammer:** stubs, jogs, and facing meander steps stay perpendicular so they add no parallel run.

### 6.2 2-D distributed CPW pair coupling

Extracted with [`code/resonator/cpw_pair_cap2d.py`](../code/resonator/cpw_pair_cap2d.py) (finite-volume Laplace solver, zero-thickness metal, 500 µm Si substrate; tabulated as `CPW_PAIR_2D` in `purcell_readout_model.py`, grid-converged within 1.3 % between 1 µm and 0.25 µm; take single-line self `C'` from the conformal value):

| Centre spacing (µm) | 150 | 200 | 300 | 400 | 600 | 800 |
|---|---:|---:|---:|---:|---:|---:|
| `c'_m` 10/6 µm, ε_r 11.9 (fF/mm) | 0.400 | 0.229 | 0.106 | 0.062 | 0.028 | 0.015 |
| `c'_m` 16/9 µm, ε_r 11.9 (fF/mm) | 0.990 | 0.564 | 0.260 | 0.151 | 0.068 | 0.036 |

Because `k_L ≈ k_C`, the inductive contribution is comparable to the capacitive one and cannot be neglected when bounding parallel-run `J`.

### 6.3 Modelling caveats

* Qubit-loaded channel lengths require a qubit-loaded Q3D end capacitance `C_end` for the specific CPW cross-section.
* HFSS uses thin PEC straps rather than physical 3-D air bridges (real air-bridge capacitance `C_AB` adds slight loading).
* Legacy side-by-side pairs (`TapJunction` + `HammerCap`) required coupled-line sections, full-Q3D mutual corrections, and claw-length offsets (`CLAW_DL`); `PurcellReadout` avoids those corrections by taking `C_end` directly in `length_for_f`.
