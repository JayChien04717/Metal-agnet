# Readout channel with an individual Purcell filter: methods and layout

General methods for any Design (Heinsoo / Wallraff / Cai 2024 style channel: feed line - C_p - lambda/4 filter -
C_j hammer - lambda/4 resonator). Design-specific numbers live in `Design<N>/notes/` and `Design<N>/README.md`.

References: PhD thesis (ETH Zurich, Wallraff group), sec. 2.5 "Qubit Readout", pp. 119–127 (Fig. 2.29, eqs.
2.119–2.137, Table 2.4); Cai et al., Chin. Phys. B 33, 068501 (2024). The local page screenshots and their summary
(`note/wallraff reference/`) were removed on 2026-11-19; the parts this project uses are in section 2 and in git
history.
Operating procedure (order of work, checklist): [`skill/purcell_resonator_design.md`](../skill/purcell_resonator_design.md).
Since 2026-11-19 new channels use the shared QComponent `pattern design/resonator/purcell_readout.py` and its circuit
model [`code/resonator/purcell_readout_model.py`](../code/resonator/purcell_readout_model.py). The layout keeps the
hammer as the only significant filter-resonator coupling (section 6). The TL frequency (section 1), the closed-form J
(section 2) and the inversions (sections 2-4) below are used as they are. The coupled-line, full-Q3D-mutual and
`CLAW_DL` corrections are needed only for the archived side-by-side pairs.

## 1. First estimates of couplers and lengths

1. Coupling capacitor for a bandwidth (lab formula, FXTFX24 cell 42): Q_c = pi / (Z0 Zr (Cc w)^2)
   -> Cc = sqrt(pi / (Q_c Z0 Zr)) / w. A two-port (transmission) filter loaded equally by both ports has Q_c = 2 Q_F.
2. Interdigital capacitance (Alley / Bahl), N fingers, overlap l (um), finger width W, substrate h:
   C [pF] = (eps_r + 1) l [(N - 3) A1 + A2], A1 = 4.409e-6 tanh(0.55 (h/W)^0.45), A2 = 9.92e-6 tanh(0.52 (h/W)^0.5).
   Replace by Q3D once the layout exists.
3. Loaded line length: f = v / (2 L) (1 - 2 Cc / (C' L)) for lambda/2 with both ends loaded; iterate the meander on
   the filleted length; check f with HFSS (the cap metal is not counted).
4. Capacitor tuning in Q3D: one run per geometry change, no sweeps. Two points fix the trend:
   hammer C_j ~ a / gap + b; enclosing cap C ~ linear in plate width.

### Quarter-wave frequency and short-tail trim

Use this method for the quarter-wave Purcell filter and readout resonator layouts in Design4, Design5, Design6, Design7, and Design8. Use the actual boundary conditions for any other topology, including the half-wave case above. Take the CPW cross section, phase velocity, filleted lengths, and loads from the selected geometry and solver configuration of each Design.

1. **Estimate the unloaded line frequency.** For a uniform line with one ideal open and one ideal short, $f_{\lambda/4}=v/(4L)$, where $L$ is the filleted CPW centre-line length between the termination reference planes and $v$ is its phase velocity. For a small excess shunt capacitance $C_{\mathrm{end}}$ at the open end, use $L_{\mathrm{eff}}\approx L+C_{\mathrm{end}}/C'$, where $C'=1/(Z_0 v)$. Use this approximation for a local open-end load. Include longer branches as transmission-line sections.
2. **Solve the circuit frequency with the actual branches.** Split each line at its tee and retain the open-side arm, shorted arm, and hammer arm at their drawn electrical reference planes. With the $e^{+j\omega t}$ convention, use $Y_0=1/Z_0$ and $\theta=2\pi f\ell/v$.

   $$
   Y_{\mathrm{open}}=jY_0\tan\theta,\qquad
   Y_{\mathrm{short}}=-jY_0\cot\theta,\qquad
   Y_{\mathrm{in}}=Y_0
      \frac{Y_L+jY_0\tan\theta}{Y_0+jY_L\tan\theta}.
   $$

   Transform each termination admittance along its arm and sum the arm admittances at the tee. Use the Q3D-derived $C_p$ and $C_j$ π networks at their reference planes and deembed any CPW capacitance already represented by the line sections. The filter sees the feed-side $C_p$ loading, and the resonator sees the open-end loading of the selected fixture, including the claw and qubit when present. Find a bare frequency from the zero of its driving-point admittance near the expected fundamental. For a lossless branch, this is a susceptance zero. With feed loading, use the real part of the complex root. In the two-mode circuit, suppress the off-diagonal transfer when calculating the bare modes and retain the chosen diagonal loads. For an isolated eigenmode comparison, use the loads corresponding to the actual isolated geometry. Restore the coupling and solve for the hybrid poles when predicting the coupled spectrum.
3. **Choose a reference for a length trim.** Use the HFSS frequency for the same topology, termination loads, and fixture. A resonator-only eigenmode can anchor an isolated resonator calculation. For the assembled channel, identify the modes with loaded eigenmode fields or recover the bare frequencies through the inversion of resolved complex S21 poles. A hybrid dip frequency and a bare frequency serve different design targets. Check frequency stability with the adaptive history before translating a model–HFSS difference into a physical trim. A persistent difference between isolated and assembled fixtures needs a configuration-matched comparison before it can set a trim.
4. **Calculate the local trim for the segment being changed.** Hold the tee, hammer, other arms, and termination geometry fixed. In the circuit, perturb the selected shorted arm by a small length $h$ and estimate $s=\partial f/\partial L_{\mathrm{short}}$ from $[f(L_{\mathrm{short}}+h)-f(L_{\mathrm{short}}-h)]/(2h)$. For a small target shift use

   $$
   \Delta L_{\mathrm{short}}\approx \frac{f_{\mathrm{target}}-f_{\mathrm{ref}}}{s}.
   $$

   Record $s$ in MHz/µm and the signed trim in µm. A negative trim shortens the arm and typically raises its frequency. For a simple quarter-wave path with fixed end loading, $s\approx-f_{\mathrm{ref}}/L_{\mathrm{eff}}$ is a first estimate. Recalculate the slope after changing a tee, hammer, termination, or CPW cross section. Confirm the resulting frequency in the matching HFSS fixture.

## 2. Wallraff-group circuit model (no curve fit)

* **Circuit (thesis Fig. 2.29, eqs. 2.122–2.136), written as admittances.** Every "C to ground" is a shunt; C_AB is
  the air-bridge capacitance to ground (value not given).
  * Resonator, qubit side: C_qr in series with C_qg, plus C_qrg + 2 C_AB at the open end, line L_rr to the hammer T
    (+ C_AB). Ground side: line L_rrg, 2 C_AB, line L_rre, short.
  * Hammer: resonator node (+ C_AB + C_prg), series C_pr, filter side (C_prg + C_AB).
  * Filter from the feed line: C_AB + C_pff, series C_pf (hugger), C_pfg + C_AB, line L_pf, shunt hammer branch,
    lines L_pfg, L_pfe, short -> branch admittance Y_rc at the feed-line T.
  * Feed line: L_in, C_in (~50 fF), feed sections, shunt Y_rc, L_out; single channel S21 = 2 / (2 + Z0 Y_rc).
* **Design scale from Table 2.4** (five channels, 6.7–7.4 GHz): C_pr ≈ 2.4–2.8 fF gives J_sim ≈ 15–16 MHz (measured
  1.05–1.30× higher); κ_p ≈ 35–63 MHz, so κ_p / 4 < J and the two hybrid dips with the ω_r peak between them are
  visible; ω_p − ω_r ≈ +11 … +23 MHz. Most of the hugger's ground capacitance is on the **filter** side
  (C_pfg ≈ 80–100 fF, C_pff ≈ −6 fF): put the frame on the filter side. The Q35 column used by Design5 W1 is in
  `Design5/layout/design5_layout.py`.
* **ω_p, ω_r, J, κ_p without a curve fit** (eq. 2.137 = 2-mode model,
  |S21| = |B (1 − A e^{iφ} 2iκ_p(ω_r − ω) / [4J² + (2i(ω_p − ω) + κ_p) 2i(ω_r − ω)])|, internal loss zero):
  * ω_r = real zero of Y_rc (S21 = B, the peak between the two dips);
  * l1, l2 = complex poles of S21 (zeros of 2 + Z0 Y_rc, complex secant) = eigenvalues of [[ω_p + iκ_p/2, J], [J, ω_r]];
  * ω_p = Re(l1 + l2) − ω_r, κ_p = 2|Im(l1 + l2)|, J² = (ω_p + i·Im(l1+l2)) ω_r − l1 l2. Im J² ≈ 0 is the consistency check.
  * Sign: the Y model uses e^{+jωt}, so the poles have Im > 0. Using ω_p − iκ/2 gives Im J² ~ 10⁵: wrong.
* **Reproducing Table 2.4** (re-examined 2026-11-17, run [`wallraff_check_adfc47c2`](../Design6/results/resonator/20261117_115343_wallraff_check_adfc47c2/run.json)). The old full-circuit reproduction (`Design6/simulation/wallraff_model.py`, git `2eb2bef`, removed in `fafa3ac`) used one global v = 1.147e8 m/s (their CPW is not given) and C_AB = 0.
  * **J:** 4–10 % below J_sim per channel (14.0–14.9 vs 15.0–15.8 MHz). This is the only quantity it reproduces.
  * **ω_p and ω_r: not reproduced.** ω_r is 20–29 MHz low and ω_p is 5–64 MHz high in every channel. ω_p − ω_r comes out at +45 to +108 MHz against +11 to +23 MHz in the table (4.1–6.2×).
  * **κ_p: not reproduced.** Per-channel errors run from −55 % to +45 % with no common sign. The earlier summary "27–51 vs 35–63 MHz" compared sorted ranges and hid this. The table values come from the full feed line (C_in, standing waves), which the single-channel model leaves out.
  * **Implication:** an unknown v, the air bridges and the lumped hugger make f_p uncertain by tens of MHz. Take f_p and f_r from EM, not from this model.
* **Closed-form J check** ([`code/general/purcell_j.py`](../code/general/purcell_j.py), also `purcell_readout_model.j_hammer`):
  * Formula: J = (f/2)·C_j·cos θ_p·cos θ_r / sqrt(C_eff,p·C_eff,r), with each open-end load converted to a line extension x = C_end/C′ and C_eff = C′(L + x)/2.
  * **Table 2.4:** it reproduces J_sim within ±5 % (15.1–16.2 vs 15.0–15.8 MHz). The measured J is 5–30 % above J_sim, so even the Wallraff model underestimates the real coupling slightly.
  * **Use it as the first sanity check of any circuit model.** If an EM J exceeds the closed form by much more than ~30 %, a coupling path other than C_j exists. Find that path before tuning C_j (section 6).

## 3. Fitting HFSS S21 with the circuit model

1. Model: lossless ABCD of P1 - TL - C_in(pi) - TL - n0 - TL - P2 with the shunt branch at n0:
   stub - C_p(pi) - filter (lambda/4, open at C_p, shorted) - tap at D_cj - C_j(pi) - resonator (lambda/4).
   Pi network = Q3D series C + the two electrode-to-ground C (C_ii + C_ij). Z0 = 50.9 Ohm, v = 1.2016e8 m/s.
2. Free parameters: dLf, dLr (effective length corrections, mm), s_in, s_p, s_j (scales on the Q3D series C;
   s_in / s_p also scale the ground C), kv (scale on v). Soft Gaussian priors: s_in, s_p sigma 0.15, s_j 0.4, kv 0.03.
3. Residual: complex S21 and 0.3 x complex S11 in linear scale (not dB, so the lossless -60 dB nulls do not dominate).
   Optional eigenmode constraint: the model eigenfrequencies with the feed cut at n0 (both n0 half lines open,
   0.15 mm each) are the zeros of Im Y(n0) (Foster: B crosses - -> +); residual (f_model - f_HFSS) / 10 MHz.
4. Optimiser: differential_evolution (global) then least_squares.
5. Derived numbers: kappa_p = FWHM of the |S21|^2 notch with C_j -> 0; J = half the minimum splitting of the
   two eigen-cut modes when dLr is scanned through the anticrossing.
   Mode identification: shift dLr / dLf by 20 um and see which mode moves.
Caveat: when the resonator is far detuned (Delta >> J) the data constrain J only through Delta shift J^2/Delta
(< 1 MHz), so J then comes from the prior (Q3D C_j). J is measured only once the two dips are resolved near resonance.

Practical points:
* For the direct eq. 2.137 fit to an e^{+jwt} complex S21 (HFSS / ABCD), fit conj(S21); otherwise kappa runs negative.
* Parametrise by the two S21 zeros z1, z2 (w_p = z1 + z2 - w_r, J^2 = w_p w_r - z1 z2) so the sharp nulls are fixed.
* kappa_p on a chip depends on the real environment at the T (C_in on one side, standing waves): about 1.8x the
  matched-line value. Compare like with like.
* Whole-channel Q3D: render shorted ends open (no net merges with ground); use only the mutual C.

## 4. Two-mode inversion of eigenmodes (J, kappa_p from HFSS eigenmode)

From the two complex eigenfrequencies z1, z2 near the crossing, solve [[f_p + i kappa_p/2, J], [J, f_r]] per point
(least squares on the trace and the determinant; `point_J` in the eigenmode scripts). A 3-point scan fit with fixed
f_p fails for weak coupling; the per-point inversion does not. Identify modes by shifting one length by 20 um.

## 5. Lessons that hold for every channel

* **Ground bonds on every CPW**, including hammer stubs, tee arms and short stubs. Unbonded, the etched stub + hammer
  pocket cuts the ground; the return-current detour gives a coupling Q3D cannot see, and an S21 / eigenmode fit puts it
  into C_j (J can come out 3x too large). qiskit auto bonds only cover
  straight segments longer than wb_threshold; add straps on the short pieces.
* **Side-by-side runs couple directly.** In the old separately routed pairs, long parallel filter–resonator runs
  at 200–300 um gave ~1-2 MHz (Design4: 6 of 12 MHz) that no capacitor in the model represents. The fix is the
  layout, not a spacing limit: avoid long parallel runs (section 6) and estimate any remaining run with
  `purcell_readout_model.parallel_j`.
* **Driven S21 does not see a far-detuned resonator** (kHz linewidth, coarse mesh along the undriven line); a fitted C_j
  then only describes the missing line. Get J from eigenmode, or from S21 with f_p = f_r.
* Eigenmode kappa_p came out ~1.5x the driven / Q3D value at 10/6 um CPW; driven S21 matched Q3D C_p within
  1-3 %. Trust driven / Q3D for C_p.
* AEDT can refuse a Mesh Link silently (ImportMesh stays false): check it, or seed the driven mesh directly.
* A full chip with several channels at a 10 um seed mesh does not fit in 24 GB; tune per channel or per pair.

## 6. PurcellReadout layout: hammer-only coupling

The goal of [`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py) is that **frequency comes from
length** and **filter–resonator coupling comes from the hammer C_j only**. Then J follows from one small Q3D coupon,
and f_p, f_r follow from the filleted length plus Q3D end loads. Design-specific numbers and runs are in
[Design8](../Design8/README.md).

### 6.1 Where unwanted coupling comes from

Earlier Designs showed three paths besides the hammer (Design4 notes, section 5):

| Path | Seen in | Size |
|---|---|---|
| Side-by-side CPW runs (capacitive + inductive) | Design4 at 200 µm, ~1.5 mm of parallel run | 6.0 MHz of 12.0 MHz HFSS J; the hammer gave 2.9 MHz |
| Ground cut by stub + hammer pocket (return current) | Design4 at 250 µm without straps | J 68 MHz → 23 MHz with 2 straps → 20 MHz with straps on every CPW |
| Meander self/neighbour coupling, parallel runs at 230–300 µm | earlier meander layouts | 1–2 MHz |

### 6.2 Layout rules in the module

1. **Opposite-side routing.** Each main line passes straight through its T. The filter's open arm goes to C_p on
   one side; its shorted arm turns away from the resonator (jog, then a meander on the far side). The resonator is
   the point mirror. The only side-by-side run is the T region: 2 × (`tee_half` + `lead`) = 300 µm at the default
   settings.
2. **Line spacing at the T** (`line_spacing`, default 400 µm; no fixed minimum). With opposite-side routing the
   T region is the only side-by-side run, so the spacing is set by the stub length the hammer needs; check the
   run with `parallel_j` (0.30 MHz, 2.2 % of J at 400 µm). At 400 µm the ground strip between the lines is 378 µm
   wide (10/6 µm CPW) and is cut only by the stubs and the hammer pocket.
3. **Meanders** (`jog` and `meander_pitch`, default 400 µm each): `jog` puts the first meander run 800 µm from the
   other line. These are defaults, not limits; the meander runs belong to one line, so they do not add
   filter–resonator coupling.
4. **Ground straps / air bridges at the hammer** (`strap_points()`): one across each stub, halfway between the T
   cut and the pocket, and one across each main line on both sides of each T (Wallraff C_AB). That is six bridges
   per module. In HFSS they are PEC straps; on the mask they are air bridges. Every other CPW gets qiskit auto bonds
   every 350 µm (skill.md §13).
5. **Perpendicular pieces only near the hammer.** The stubs run perpendicular to the main lines. Jogs and meander
   steps that face the other line are perpendicular or end-on, so they add no parallel run.

Alternatives that were considered:
- **Perpendicular main lines** (filter T rotated 90°): no parallel run at all, but the hammer bar must bend, so
  the coupon is no longer a straight two-pin piece.
- **Wider spacing:** at 600 µm, c′_m halves again (table below), but the longer stubs move the hammer away from
  the lumped-element regime.
- **Ground fingers between the T regions:** not needed. The 378 µm strip is already solid ground. What matters is
  that it stays one piece, which the straps ensure.

### 6.3 Numbers

**Distributed coupling (2-D extraction).** Cross-section values from
[`code/resonator/cpw_pair_cap2d.py`](../code/resonator/cpw_pair_cap2d.py); the full rows (c′_m, k_C, k_L) are
`CPW_PAIR_2D` in [`purcell_readout_model.py`](../code/resonator/purcell_readout_model.py), which records the run
provenance. The solver is a finite-volume Laplace solve with zero-thickness metal, a 500 µm substrate and natural
outer boundaries. As a check, it gives 0.546 fF/mm for a 16/9 µm pair at 200 µm (ε_r 11.5), against 0.563 fF/mm from
the Design4 Q3D pair coupons. Its c′_m changes by less than 1.3 % between grid steps of 1 µm and 0.25 µm. The
single-line self C′ converges slowly, about 5 % high at the finest grid, so take C′ from the conformal value.

| Centre spacing (µm) | 150 | 200 | 300 | 400 | 600 | 800 |
|---|---|---|---|---|---|---|
| c′_m 10/6 µm, ε_r 11.9 (fF/mm) | 0.400 | 0.229 | 0.106 | 0.062 | 0.028 | 0.015 |
| c′_m 16/9 µm, ε_r 11.9 (fF/mm) | 0.990 | 0.564 | 0.260 | 0.151 | 0.068 | 0.036 |

k_L ≈ k_C (for example 3.3e-4 against 3.6e-4 at 400 µm, 10/6 µm). The inductive part is therefore as large as the
capacitive part and cannot be neglected.

**Design8 module (10/6 µm, 400 µm, f 6.5 GHz).** Q3D run
[`20261118_232405_q3d_module_b3d273d0`](../Design8/results/resonator/20261118_232405_q3d_module_b3d273d0/summary.json)
used 8 passes. No case reported solver convergence; last-pass changes were 0.8–4.9 %.

| Quantity | Value | Source |
|---|---|---|
| Hammer C_j (coupon mutual − straight-line reference mutual) | 1.971 fF | Q3D (fit to Design4 data: 1.908 fF) |
| Full-module static mutual (ends rendered open) | 2.076 fF | Q3D |
| Static non-hammer mutual | 0.105 fF (5.3 % of C_j) | Q3D difference |
| Same, side-by-side runs only (2-D c′_m × run lengths) | 0.037 fF | 2-D model |
| J from the hammer (closed form, Q3D C_j) | 13.19 MHz | model |
| Parallel-run J bound (\|J_C\| + \|J_L\|) | ≤ 0.30 MHz (2.2 %) | model |

The static non-hammer mutual (0.105 fF) is about 3× the 2-D side-by-side estimate. The rest is near-field coupling
around the T regions and pocket, plus Q3D mesh noise: the straight reference alone has 0.032 fF of mutual, and its
two lines differ by 2.6 fF to ground although they are symmetric. Static mutual capacitance also overstates J,
because most of the extra part sits where the mode voltages are small.

**Design8 HFSS check** (open fixture, lengths from the Q3D-loaded model, nothing fitted; [Design8](../Design8/README.md#2-validation-results-open-fixture-revision-2026-11-18_r1)):

| | Model | HFSS driven, 6 passes | HFSS eigenmode, 6 passes |
|---|---|---|---|
| J | 13.19 MHz (hammer only) | 13.39 MHz (×1.02) | 16.71 MHz (×1.27) |
| f_p / f_r | 6.500 / 6.500 GHz | +8.2 / −9.4 MHz | +10.0 / +15.8 MHz |
| κ_p | 26.8 MHz | 25.6 MHz | 47.2 MHz |

The driven J matching the hammer-only value within 2 % is the direct evidence that the layout removes the other
coupling paths. In Design4 (old pair) the same comparison gave HFSS 12.0 MHz against 2.9 MHz from the hammer.
The bound of 0.30 MHz from the side-by-side runs (2.2 %) is consistent with the residual.

### 6.4 What is not simulated

- The qubit-loaded channel. The lengths for the fluxonium and transmon chips use the Design3 Q3D pad capacitances,
  which were extracted for the old 15/5.9 µm CPW. Run a qubit-loaded Q3D before using those lengths.
- Air-bridge loading on the mask. HFSS uses thin PEC wires, not real bridges.
- The perpendicular-routing variant, and line spacings other than 400 µm, in Q3D or HFSS.
- Q3D convergence beyond 8 passes.
