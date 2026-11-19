# Readout channel with an individual Purcell filter: methods

General methods for any Design (Heinsoo / Wallraff / Cai 2024 style channel: feed line - C_p - lambda/4 filter -
C_j hammer - lambda/4 resonator). Design-specific numbers live in `Design<N>/notes/`.
References: `note/wallraff reference/` (Wallraff-group thesis sec. 2.5, Table 2.4; Cai et al., Chin. Phys. B 33, 068501 (2024),
password-protected PDF: read with `pdftotext -layout`).
Operating procedure (order of work, checklist): [`skill/purcell_resonator_design.md`](../skill/purcell_resonator_design.md).
Since 2026-11-19 new channels use the shared module `pattern design/resonator/purcell_readout.py`. Its layout keeps
the hammer as the only significant filter-resonator coupling ([layout note](purcell_readout_parasitics.md)). The TL
frequency (section 1), the closed-form J (section 2) and the inversions (sections 2-4) below are used as they are.
The coupled-line, full-Q3D-mutual and `CLAW_DL` corrections are needed only for the archived side-by-side pairs.

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

Use this method for the quarter-wave Purcell filter and readout resonator layouts in Design4, Design5, Design6, and Design7. Use the actual boundary conditions for any other topology, including the half-wave case above. Take the CPW cross section, phase velocity, filleted lengths, and loads from the selected geometry and solver configuration of each Design.

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

The reference pages are in `note/wallraff reference/` (thesis sec. 2.5, Table 2.4).

* **Circuit:** eqs. 2.122–2.134, written as admittances. Branch admittance Y_rc at the feed-line T; single channel S21 = 2 / (2 + Z0 Y_rc).
* **ω_p, ω_r, J, κ_p without a curve fit** (eq. 2.137 = 2-mode model):
  * ω_r = real zero of Y_rc (S21 = B, the peak between the two dips);
  * l1, l2 = complex poles of S21 (zeros of 2 + Z0 Y_rc, complex secant) = eigenvalues of [[ω_p + iκ_p/2, J], [J, ω_r]];
  * ω_p = Re(l1 + l2) − ω_r, κ_p = 2|Im(l1 + l2)|, J² = (ω_p + i·Im(l1+l2)) ω_r − l1 l2. Im J² ≈ 0 is the consistency check.
  * Sign: the Y model uses e^{+jωt}, so the poles have Im > 0. Using ω_p − iκ/2 gives Im J² ~ 10⁵: wrong.
* **Reproducing Table 2.4** (re-examined 2026-11-17, run [`wallraff_check_adfc47c2`](../Design6/results/resonator/20261117_115343_wallraff_check_adfc47c2/run.json)). The old full-circuit reproduction (`Design6/simulation/wallraff_model.py`, git `2eb2bef`, removed in `fafa3ac`) used one global v = 1.147e8 m/s (their CPW is not given) and C_AB = 0.
  * **J:** 4–10 % below J_sim per channel (14.0–14.9 vs 15.0–15.8 MHz). This is the only quantity it reproduces.
  * **ω_p and ω_r: not reproduced.** ω_r is 20–29 MHz low and ω_p is 5–64 MHz high in every channel. ω_p − ω_r comes out at +45 to +108 MHz against +11 to +23 MHz in the table (4.1–6.2×).
  * **κ_p: not reproduced.** Per-channel errors run from −55 % to +45 % with no common sign. The earlier summary "27–51 vs 35–63 MHz" compared sorted ranges and hid this. The table values come from the full feed line (C_in, standing waves), which the single-channel model leaves out.
  * **Implication:** an unknown v, the air bridges and the lumped hugger make f_p uncertain by tens of MHz. Take f_p and f_r from EM, not from this model.
* **Closed-form J check** ([`code/general/purcell_j.py`](../code/general/purcell_j.py)):
  * Formula: J = (f/2)·C_j·cos θ_p·cos θ_r / sqrt(C_eff,p·C_eff,r), with each open-end load converted to a line extension x = C_end/C′ and C_eff = C′(L + x)/2.
  * **Table 2.4:** it reproduces J_sim within ±5 % (15.1–16.2 vs 15.0–15.8 MHz). The measured J is 5–30 % above J_sim, so even the Wallraff model underestimates the real coupling slightly.
  * **Use it as the first sanity check of any circuit model.** If an EM J exceeds the closed form by much more than ~30 %, a coupling path other than C_j exists. Find that path before tuning C_j.

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
* **Every CPW pair >= 400 um apart** (skill.md section 11): closer runs give a direct filter-resonator coupling of
  ~1-2 MHz that no capacitor in the model represents.
* **Driven S21 does not see a far-detuned resonator** (kHz linewidth, coarse mesh along the undriven line); a fitted C_j
  then only describes the missing line. Get J from eigenmode, or from S21 with f_p = f_r.
* Eigenmode kappa_p came out ~1.5x the driven / Q3D value at 10/6 um CPW; driven S21 matched Q3D C_p within
  1-3 %. Trust driven / Q3D for C_p.
* AEDT can refuse a Mesh Link silently (ImportMesh stays false): check it, or seed the driven mesh directly.
* A full chip with several channels at a 10 um seed mesh does not fit in 24 GB; tune per channel or per pair.
