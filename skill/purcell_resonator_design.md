# Designing a Purcell filter + readout resonator channel

Operating procedure for a channel of the form feed line – C_p – λ/4 Purcell filter – C_j hammer – λ/4 readout
resonator – qubit. **Since 2026-11-19 every new channel uses the shared module
[`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py)** with its model
[`purcell_readout_model.py`](../pattern%20design/resonator/purcell_readout_model.py). The layout keeps the
hammer as the only significant filter–resonator coupling ([layout note](../note/purcell_readout_parasitics.md)).
With that layout, one Q3D coupon gives J, and length plus Q3D end loads give f_p and f_r.

- The general methods (TL frequency, short-tail trim, Wallraff model, two-mode inversion) are in
  [note/purcell_filter_method.md](../note/purcell_filter_method.md).
- The reference implementation and validation is [Design8](../Design8/README.md). Design5 and Design6 use the same
  module (`Design5/layout/chip6_prr.py`, `Design6/layout/design6_layout.py` `readout='prr'`).
- Old channels built from TapJunction + HammerCap / InterdigitalCap with side-by-side lines are archived and not
  recommended ([archive](../pattern%20design/archive/README.md), section 7 below).

## 1. What sets what

| Quantity | Set by | How to get it |
|---|---|---|
| f_p | filter filleted length + end loads (C_p π network, hammer arm) | `length_for_f(f, tap, C_p mutual + C_p plate ground, s/2, C_j + head excess)` |
| f_r | resonator filleted length + end loads (qubit pad or claw, hammer arm) | `length_for_f(f, tap (+ open extension), C_end, s/2, C_j + head excess)` |
| J | hammer C_j only | `j_hammer(...)` with the Q3D coupon C_j |
| κ_p | C_p to the feed line and the feed environment (C_in, standing waves) | `kappa_hanger` as a first estimate; HFSS driven for the value |
| κ_r, χ | J, κ_p, Δ, qubit coupling | `kappa_readout(J, κ_p, Δ)` |

Design intent: f_p = f_r at the loaded operating point. Taps (1.2 / 1.3 mm from the open ends), line spacing
(400 µm) and the hammer set J. The lengths set the frequencies and nothing else.

## 2. Layout rules (built into the module)

1. The two main lines pass straight through their own T. The filter and resonator then leave in opposite
   directions. The only side-by-side run is 2 × (`tee_half` + `lead`) = 300 µm at `line_spacing` ≥ 400 µm.
2. Meanders keep ≥ 400 µm to everything (`jog`, `meander_pitch`). Use `*_first_run='back'` for tall, narrow
   channels and `mirror` to put the filter on the other side.
3. Air bridges / ground straps go at `strap_points()`: across both hammer stubs and across both main lines on each
   side of each T. In HFSS and Q3D they are PEC straps. Every other CPW gets bonds every 350 µm (skill.md §13).
4. The qubit's readout pin abuts the `qubit` pin. If it cannot (a side claw), route a short extension and declare
   its filleted length in `resonator_open_extension`. The extension counts in the tap distance and in
   resonator_length.
5. C_p (EnclosingCap pin b) abuts the `cp` pin. Put C_in after the input launch when the chip has one
   (Design5, Design6).

## 3. Efficient order of work

1. **Q3D coupons (minutes, 8 passes)**, with the shared drivers (`code/resonator/prr_sim.py`;
   `Design<N>/simulation/resonator/prr_validate.py q3d`):
   - `cj`: module `draw_mode='coupon'`. `ref`: two straight lines with the same outer ends.
     C_j = mutual(cj) − mutual(ref); head excess = ground(cj) − ground(ref) − C′·s/2.
   - `cp`: C_p alone. The filter end load is the mutual plus the plate-side ground.
   - Resonator end load, with the qubit present and its island floating (Kron): the claw case for an Xmon claw;
     the readout-pad Q3D for a fluxonium or a transmon pad.
   - `module`: the whole module with `short_ends='open'`. mutual − C_j is the non-hammer coupling. Keep it below
     ~10 % of C_j (Design8: 5.3 %).
2. **Lengths from the model** (no solver): `length_for_f` per line with the Q3D loads above. Check the drawn
   `length_report` against the targets (exact by construction).
3. **HFSS check of the same geometry** (`prr_validate.py eigen` / `driven`, or Design8 `hfss_module.py`):
   - Eigenmode with 6 passes: two complex modes per channel → `invert` → f_p, f_r, J, κ_p.
   - Driven S21: the zeros + S21 = 1 point method, and a pole fit with the zeros fixed (`s21_channel`), per
     channel window.
   - Driven pass budget: 6 passes, as an explicit deviation from skill.md §12 that the run records. It matches the
     eigenmode mesh (~525k elements for one channel). In Design8, 3 passes left the dips 20 MHz off, and 10 passes
     ran out of memory.
   - Report the convergence as recorded. Compare with the prediction without fitting.
4. **Trim if needed:** ΔL = Δf / s with `trim_slope` (≈ −1.5 MHz/µm at 6.5–7 GHz). Change only the shorted arm
   length (the meander absorbs it). The hammer and taps stay fixed, so J does not change.
5. **Loaded check** with the qubit (junction as lumped L) for χ and κ_r at the operating point. In Design6 the
   validation fixture is already loaded.

## 4. Expected agreement (Design8, 10/6 µm, ε_r 11.9, f 6.5 GHz)

| | Model (Q3D in) | HFSS driven, 6 passes | HFSS driven, 3 passes | HFSS eigenmode, 6 passes |
|---|---|---|---|---|
| J | 13.19 MHz | **13.39 MHz (×1.02)** | 13.11 MHz (×0.99) | 16.7 MHz (×1.27) |
| f_p / f_r | 6.500 / 6.500 GHz | +8 / −9 MHz | −12 / −18 MHz | +10 / +16 MHz |
| κ_p | 26.8 MHz | 25.6 MHz | 25.3 MHz | 47.2 MHz (eigenmode κ_p runs ~1.5–1.8× high, method note §5) |

- **ΔS:** driven ΔS did not settle (0.2–0.9 from pass to pass, skill.md §10).
- **Frequencies:** they move ~20 MHz between 3 and 6 passes, so frequency agreement is at the mesh-uncertainty
  level.
- **J and κ_p:** use driven S21; they are stable between meshes.
- **Pass budget:** 10 passes ran out of 32 GB.

The Design5 / Design6 checks are in the Design READMEs and the
[layout note](../note/purcell_readout_parasitics.md). Compare: in Design4 (old pair, 200 µm, ~1.5 mm side by
side) the hammer gave 2.9 MHz of the 12.0 MHz HFSS J.

## 5. Feed with a series C_in: tee placement and C_p compensation

A series input capacitor C_in reflects the input side of the feed. With Q3D C_in = 40 fF (Cai T-type cap), the
off-resonance |S21| is about −16 dB. Each filter's decay then interferes with its own reflection off C_in. κ_p
follows the real part of the feed impedance at the channel's C_p node. It is no longer the matched-feed (hanger)
value, and it changes along the feed with period λ/2 = v / 2f. On chip v ≈ 1.2e8 m/s, so λ/2 ≈ 8.6–10 mm at
6–7 GHz: about one period across a 10 mm chip.

**Evidence (Design5 prr_r1, 1-pass HFSS S21):**
- **κ_p from HFSS:** κ_p was 0.07–1.47× the hanger value. ch6 at 7.0 GHz sat 3.3 mm from C_in, near a node, with
  κ_p 2.6 MHz.
- **Circuit model:** the full-feed circuit model with Q3D C_in / C_p π networks and no HFSS input
  (`Design5/simulation/resonator/prr_circuit_model.py`) gave 0.7 MHz for ch6. It reproduced the pattern on all six
  channels.
- **Without C_in:** the same model gives ≈ 1.0× everywhere.

Rules:
1. **Model the whole feed, not one channel.** Compute κ_p of every channel with all channels present, using:
   - ports, C_in (series + both shunt grounds), each tee with its 150 µm lead and C_p (frame shunt on the feed
     side), and the filter–hammer–resonator network;
   - `prr_circuit_model.py` (S21 + the same pole fit as HFSS), or its fast form κ_p = κ_hanger · Re Z_node / (Z0/2)
     at the C_p frame node.

   The other channels' C_p frame shunts (119 fF each) shift the standing wave. A single-channel map is only a
   guide.
2. **Avoid the minima, not the exact nodes.** Choose tee positions that maximise the worst κ_p / κ_hanger:
   - **Robustness corners:** the worst case is taken over f × (0.99, 1, 1.01), v × (0.98, 1.02),
     C_frame × (0.8, 1.2) and C_in × (0.9, 1.1). The design must survive the 1 % frequency offsets and the model
     uncertainty.
   - **Constraints:** channel footprints (same-side pitch ≥ 1.7 mm needs taller meanders; opposite-side tees
     ≥ 0.6 mm) and the 400 µm spacing rule.
   - **Freedoms:** the frequency order along the feed and the side of each channel are free variables.
3. **Compensate with C_p.**
   - **Where:** a channel that remains below its target (ratio < 1) but is clear of a minimum.
   - **How much:** κ_p ∝ C_p², so C_p,new = C_p · √(κ_target / κ_model). Iterate with the full model, because C_p
     also loads the feed.
   - **Geometry:** get the cap geometry and its loads (series, plate and frame ground) from a Q3D coupon
     (`prr_validate.py q3d --cases cp` with `ENC_PLATE_UM`).
   - **Lengths:** recompute that channel's filter length from the new C_p loads (`length_for_f` end load
     = mutual + plate ground). The resonator length, the hammer and J do not change.
4. **Record** the placement, the per-channel C_p geometry and loads, the objective and the robustness corners in
   the revision config. Design5: `parameters/readout_prr_r2.json` `placement` and `placement_method`.
   `chip6_prr.py` reads them.
5. **Validate** with a full-feed HFSS S21 and compare κ_p per channel with the circuit model.

**Design5 prr_r2 (worked example, `parameters/readout_prr_r2.json`):**

| Channel | x (mm) | Side | C_p (fF) |
|---|---|---|---|
| ch6 (7.0 GHz) | −3.39 | below | 44.03 (207 µm plate) |
| ch3 (6.4 GHz) | −2.77 | above | 40.10 |
| ch5 (6.8 GHz) | +0.76 | above, tall meanders | 44.75 (211.5 µm plate) |
| ch4 (6.6 GHz) | +1.36 | below, tall meanders | 40.10 |
| ch2 (6.2 GHz) | +2.50 | above, tall meanders | 40.10 |
| ch1 (6.0 GHz) | +3.17 | below, tall meanders | 40.10 |

- **Lengths:** the resonator lengths are those of r1. ch5 / ch6 filters are 32 / 27 µm shorter for the larger C_p.
- **Full-feed circuit model, κ_p / κ_hanger(40.10 fF):** r1 0.02–1.47 → r2 0.92–1.33.
- **Robust worst case:** 0.64 at the robustness corners, before C_p compensation.
- **HFSS S21, 1 pass, six channels:** κ_p / target was 0.07–1.47 (r1) and is 0.73–1.19 (r2). HFSS sits 0–30 % below the model with C_in ([Design5 README](../Design5/README.md)).

The same applies to Design6. There C_in is 40 fF at 0.10 mm after P8 on a two-channel feed (not computed yet).

## 6. Checklist before quoting channel parameters

- [ ] Module from `PurcellReadout`; `length_report` equals the target lengths; `resonator_open_extension` declared if used.
- [ ] Q3D run named for C_j, head excess, C_p loads and the resonator end load (with the qubit present); module parasitic < 10 % of C_j.
- [ ] Straps at `strap_points()` in every HFSS / Q3D run, and air bridges at the same positions on the mask.
- [ ] Both dips of each channel resolved; poles grouped by channel window, not by list order.
- [ ] With a series C_in on the feed: κ_p from the full-feed model (section 5), tee positions clear of the minima, C_p compensation recorded.
- [ ] Pass budget, final ΔS / ΔF and any deviation from skill.md §12 recorded; validation status "unknown" unless convergence evidence exists.
- [ ] Model numbers labelled as model, HFSS numbers as HFSS fit; no fitted lengths or capacitances.

## 7. Legacy channels (archived pattern, not recommended)

Design4, the Design5 r3–r7 wafer and Design6 r9–r13 used TapJunction + HammerCap with the two lines side by side
(190–250 µm, 1–1.6 mm). Their results stay valid as history. The method they needed was:
- **Coupled lines:** add coupled-line sections for every side-by-side run. Take c′_m from a two-line Q3D or
  `code/resonator/cpw_pair_cap2d.py`, including the inductive part.
- **Full-Q3D mutual:** add the rest of the full-channel Q3D mutual, and any feed–resonator mutual.
- **Ground straps** on the stubs and every T.
- **J and f_p from HFSS**, using the S21 pole fit.
- **`CLAW_DL`:** a simulation without the claw had to add the claw shortening back (0.100 mm in Design6).

Without these corrections the hammer explained only ~¼ of J, and model f_p was 20–45 MHz off
([Design4 note](../Design4/notes/hammer_fixed_length_validation_20261116.md)). The old models are
`Design6/simulation/resonator/readout_model_q3d.py` and `Design4/simulation/resonator/hammer_coupled_model.py`.
With `PurcellReadout`, the end load is a model input (`length_for_f` takes C_end directly), so no `CLAW_DL`
bookkeeping is needed. A simulation without the qubit simply uses the matching C_end.
