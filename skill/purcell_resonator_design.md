# Designing a Purcell filter + readout resonator channel efficiently

Operating procedure for a channel of the form feed line – C_p – λ/4 Purcell filter – C_j hammer – λ/4 readout resonator – claw – qubit. It covers the order of work, the cheapest tool for each quantity, and the end-load bookkeeping (`CLAW_DL`).

- The [shared frequency and short-tail trim method](../note/purcell_filter_method.md#quarter-wave-frequency-and-short-tail-trim) applies to the quarter-wave layouts in Design4, Design5, Design6, and Design7. The same note also covers the Wallraff circuit model.
- Design-specific numbers belong in `Design<N>/notes/`. The worked example is [Design6 readout_s21_20261117.md](../Design6/notes/readout_s21_20261117.md).

## 1. What sets what

| Quantity | Set by | Fast tool | Notes |
|---|---|---|---|
| κ_p (filter linewidth) | C_p to the feed line | Q3D coupon, then the circuit model | HFSS κ_p ran 12 % above the Q3D model in Design6 |
| J (filter–resonator) | hammer C_j **plus** parallel-line coupling and ground return | HFSS S21 (hybridized pair) | Hammer-only models give J 3–4× too low; with coupled lines and the full-Q3D mutual, ≈ +10 % (section 2 step 3, section 4) |
| f_p | L_filter and the filter's end loads (C_p electrode ground C, tee) | HFSS S21 pole fit | see section 3 |
| f_r | L_res and the resonator's end load (claw, qubit) | HFSS S21 pole fit plus the `CLAW_DL` rule | see section 3 |
| κ_r, χ | J, κ_p, Δ = f_r − f_p, the qubit–claw g | derived | design for f_p = f_r **at the loaded operating point** |

Design intent is **f_p = f_r in the final, loaded geometry** (claw and qubit present). A detuning that shows up in a simulation means either the lengths are wrong, or the simulation's end loads differ from those the lengths were designed for. Section 3 explains the second case.

## 2. Efficient order of work

1. **Q3D coupons (minutes each, no full chip).** Solve:
   - C_p: mutual plus both ground Cs.
   - Hammer C_j: mutual, plus excess ground per side relative to a straight-CPW reference of the same arm length.
   - Claw end load: the claw to ground and to the qubit pad, with the qubit present.

   Render shorted ends open and use only the mutual C (method note §3).
2. **End-load lengths (closed form, no solver).** An open-end excess capacitance C_end shortens a λ/4 line by

   ΔL ≈ C_end / C′, with C′ = 1 / (Z0 · v)

   For 16/9 µm CPW on Si (Z0 49.9 Ω, v 1.199e8 m/s), C′ ≈ 167 fF/mm, so 17 fF ≈ 0.10 mm. Bookkeep the filter's and the resonator's end loads separately. Two lines of equal length are **not** at the same frequency even when the hammer is symmetric, because the filter end sees the C_p electrode and lead and the resonator end sees the claw.
3. **Circuit model for the first lengths only.** Use a nodal TL model of the real feed line (both channels, 50 Ω ports, Q3D pi networks); Design6 entry point: `Design6/simulation/resonator/readout_model.py`. Trust only its trends.
   - In Design6 the model's f_r − f_p was off by +24 to +41 MHz.
   - Model f_p was 38–44 MHz below HFSS, and f_r was 3–14 MHz below.
   - Model J was 3.4× low (hammer only).

   Without the full-Q3D corrections below, do not tune final lengths with the model, and do not reuse a correction factor such as `HFSS_CORR` from another Design without re-deriving it.
   - **Making the model predictive (Design6, `readout_model_q3d.py`).** Run one full readout Q3D of the same geometry (`readout_q3d_full.py`: shorted and feed ends rendered open, one net per line) and compare the model's static matrix entry by entry (`readout_q3d_check.py`). Then add, from Q3D only:
     - coupled-line sections along the parallel run (c′m from a two-line Q3D);
     - the rest of the full-Q3D filter–resonator mutual;
     - any direct feed–resonator mutual;
     - the filter-ground excess δ (filter minus per-length-scaled resonator difference; robust to Q3D mesh drift);
     - the Q3D line C′ with the conformal L′.

     For Design6 this gave f_p −4 / +1 MHz, f_r +15 / +3 MHz, J +10 % and κ_p −5 % vs HFSS, with no HFSS fitting.
4. **One readout-only HFSS driven S21, close to resonance.** Include the feed line, filters, hammers and resonators, but no qubits.
   - Give each resonator the end load the lengths assume (section 3).
   - Keep |f_r − f_p| ≲ 2J, so the S21 shows two resolved dips per channel. Then a 1 MHz wide sweep is enough and no fine sweep is needed.
   - Use 10 µm seeds and ground straps on every CPW, including the stubs (skill.md §12–13).
   - Use the [shared HFSS driven ΔS criterion](skill.md#hfss-driven-acceptance). Check pole resolution and fit quality before quoting bare frequencies or J.
5. **Four-pole fit and inversion** (`Design6/simulation/resonator/readout_poles.py --auto-guesses`).
   - Fit a cubic background plus one complex pole per mode.
   - For each channel, f_r = Im(ab)/Im(a+b), f_p = Re(a+b) − f_r and J² = f_p f_r − Re(ab), with a, b the complex poles (Im = half-width).
   - Assumptions: one lossy mode, real J, no intrinsic loss, no cross-channel coupling.
   - Pair the poles by line length (the longer pair is the lower one), not by the order of the dips in the list.
6. **Trim the selected shorted arm.** Use a frequency reference from the matching HFSS geometry and loading configuration. In the circuit, perturb that arm alone to find its local slope s = ∂f/∂L_short and use ΔL_short ≈ (f_target − f_ref) / s. The simple estimate s ≈ −f_ref/L_eff applies to a quarter-wave line with fixed end loading. Use the [shared method](../note/purcell_filter_method.md#quarter-wave-frequency-and-short-tail-trim) when tee and hammer branches affect the slope. Resolve an isolated-versus-assembled frequency discrepancy before using it to set a trim. Re-run the matching HFSS check to confirm the updated frequency.
7. **Loaded check.** Run one eigenmode or driven run with the claw and qubit (junction as a lumped L) to confirm f_p = f_r and to get χ and κ_r at the operating point.

## 3. `CLAW_DL`: matching the end load between simulation and fabrication

- `CLAW_DL` is the resonator shortening that compensates the claw's end capacitance. Design6 uses `CLAW_DL` = 0.100 mm for ≈17 fF (`Design6/layout/design6_layout.py`). `tune_lengths()` designs f_p = f_r with the claw, then sets the fabricated L_res = L_res(no claw) − `CLAW_DL`.
- **Rule: a simulation without the claw must use L_res(fab) + `CLAW_DL`.** A simulation with the claw and qubit uses L_res(fab) as drawn.
- **Record the choice:** write `res_extra_mm` / `claw_dl_mm` and the reason into run.json (`readout_s21.py --res-extra-mm`).
- **What it costs to forget:** Design6 r9 had the claw removed and the shortening kept. The resonators sat 165–197 MHz above the filters, the two modes did not hybridize, the readout dips were kHz–sub-MHz wide on a 1 MHz grid, and J could only be guessed dispersively (21.7 MHz, an overestimate). With 0.1 mm restored, the same model gave resolved pairs, J = 16.1 / 16.6 MHz and f_r − f_p = +24 / +7 MHz.
- **Back-conversion:** a length correction found in a no-claw simulation (step 6) carries over 1:1 to the fabricated L_res, as long as `CLAW_DL` stays the same.
- **Validate `CLAW_DL` itself:**
  - Q3D the claw with the qubit pad present and convert with C_end / C′.
  - Compare against a loaded HFSS run.
  - If the claw or qubit geometry changes, re-derive `CLAW_DL`; do not carry it over.
- **The filter has an analogous end load:** the C_p electrode ground C and the tee stub. If the C_p geometry changes, re-tune L_filter. Do not assume the old f_p.

## 4. Known model gaps

- **First check: closed-form J** (`code/general/purcell_j.py`). It reproduces the Wallraff Table 2.4 J_sim within ±5 %, and the Wallraff measured J was 5–30 % above it.
  - In Design6 r9 the closed form (5.0 / 4.8 MHz) equals the full circuit model (4.9 / 4.7 MHz), but HFSS gives 16.6 / 16.1 MHz.
  - The equivalent lumped C_j is therefore ≈ 2.2 fF, 3.3× the Q3D hammer (0.66 fF).
  - A ratio far above 1.3 means a coupling path outside C_j. Find it before resizing the hammer.

- **Parallel-line coupling:** at 200 µm centre spacing it gives a substantial part of J. Design6 r9 has ≈ 1.64 mm of filter–resonator side-by-side run per channel: 0.815 mm above the tee, 0.30 mm through it and 0.52 mm below it. Above the tee the filter's open end runs beside the resonator's shorted tail; below it, the reverse. Capacitive and inductive coupling therefore do not cancel as they would in an aligned (combline) pair. Q3D C′m is about 0.56 fF/mm, and including it was worth ≈6 MHz of J in Design4. Deembedding only the hammer misses it.
- **Ground return around the hammer:** without ground straps on the stubs, the etched pocket cuts the ground, and an S21 fit puts that extra coupling into C_j (see the [Purcell C_j mismatch](../note/purcell_filter_method.md) lessons, §5).
- **f_p error:** the filter frequency is the least reliable model output, because the C_p / tee junction is a 3-D region that a pi network does not capture. Take f_p from HFSS.
- **Far-detuned modes:** with Δ ≫ J, the S21 constrains J only through κ_r ≈ J²κ_p/Δ², which is sensitive to any extra decay path. Re-tune toward resonance before quoting J.

## 5. Checklist before quoting channel parameters

- [ ] The simulation's resonator end load matches the length convention (`CLAW_DL` added if the claw is absent) and is recorded in run.json.
- [ ] Both dips of each channel are resolved: FWHM > 2 grid steps (`widths_resolved` in the pole-fit output).
- [ ] Poles are paired by line length, and the labels are checked against the targets.
- [ ] Bare f_p is unchanged between runs where only L_res changed (consistency check).
- [ ] Adaptive Delta S, the pass count and the fast-sweep caveat are stated. Validation status is "unknown" unless convergence evidence exists.
- [ ] Circuit-model numbers are labelled as model, and HFSS numbers as HFSS fit, never mixed.
