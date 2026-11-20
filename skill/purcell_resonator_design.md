# Skill: Designing a Purcell filter + readout resonator channel

Operating procedure for a channel of the form `feed line – C_p – λ/4 Purcell filter – C_j hammer – λ/4 readout resonator – qubit`.

* **Standard component & model:** [`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py) with [`code/resonator/purcell_readout_model.py`](../code/resonator/purcell_readout_model.py). The opposite-direction layout keeps the hammer as the only significant filter–resonator coupling so Q3D coupon C_j sets J and filleted lengths + Q3D end loads set f_p and f_r.
* **Physics & formulas:** [note/purcell_filter_method.md](../note/purcell_filter_method.md).
* **Solver defaults:** [em_simulation.md](em_simulation.md) (Q3D 8 passes, HFSS eigenmode 6 passes / 10 µm, HFSS driven 4 passes with explicit length mesh). Design-specific pass overrides (e.g. for large multi-channel chips) belong in that Design's `README.md` / `notes/`.
* **Legacy parts:** old side-by-side pairs (`TapJunction` + `HammerCap` / `InterdigitalCap`) are archived and not recommended ([archive README](../pattern%20design/archive/README.md)).
* **Run naming:** use `purcell_resonator_pair_*` tags for new runs (earlier immutable run folders named `*_prr_*` remain unchanged).

---

## 1. What sets what

| Quantity | Set by | Model helper (`purcell_readout_model.py`) |
|---|---|---|
| f_p | Filter filleted length + end loads (C_p π network, hammer arm) | `length_for_f(f, tap, C_p mutual + C_p plate ground, s/2, C_j + head excess)` |
| f_r | Resonator filleted length + end loads (qubit pad/claw, hammer arm) | `length_for_f(f, tap (+ open extension), C_end, s/2, C_j + head excess)` |
| J | Hammer C_j only | `j_hammer(...)` with Q3D coupon C_j |
| κ_p | C_p to feed line + feed environment (C_in, standing waves) | `kappa_hanger` (first estimate); full-feed model / HFSS driven for actual value |
| κ_r, χ | J, κ_p, Δ, qubit coupling | `kappa_readout(J, κ_p, Δ)` |

Design intent: `f_p = f_r` at the loaded operating point. Taps, `line_spacing` (default 400 µm), and the hammer set J; the shorted-arm lengths set frequencies without altering J.

---

## 2. Layout rules (built into `PurcellReadout`)

1. **Opposite-direction routing:** both main lines pass straight through their T and leave in opposite directions. The only side-by-side run is `2 × (tee_half + lead)` at `line_spacing` (default 400 µm, not a hard limit; check with `parallel_j`).
2. **Meanders:** `jog` and `meander_pitch` default to 400 µm. Use `*_first_run='back'` for tall, narrow footprints and `mirror` to flip the filter side.
3. **Ground straps / air bridges:** place at `strap_points()` (both hammer stubs + both sides of each T; PEC straps in HFSS/Q3D, air bridges on mask). Every other CPW gets wire bonds every 350 µm ([em_simulation.md §3](em_simulation.md#3-ground-wire-bonds-in-every-hfss--q3d-run)).
4. **Pin connections:** abut the qubit's readout pin to `qubit` (if using a side-claw extension, declare its filleted length in `resonator_open_extension` so it counts in tap distance and resonator length). Abut `EnclosingCap` pin `b` to `cp`.

---

## 3. Order of work

1. **Q3D coupons (8 passes, [em_simulation.md §2](em_simulation.md#2-default-pass-counts-and-initial-mesh))**, via `code/resonator/purcell_resonator_pair_sim.py`:
   - `cj` (`draw_mode='coupon'`) vs. `ref` (two straight lines with same outer ends): `C_j = mutual(cj) − mutual(ref)`; `head_excess = ground(cj) − ground(ref) − C'·s/2`.
   - `cp` (`EnclosingCap` alone): filter end load = `mutual + plate_ground`.
   - Resonator end load with the qubit present and its island floating (Kron-reduced): claw or readout-pad Q3D.
   - `module` (`short_ends='open'`): verify non-hammer mutual (`mutual − C_j`) stays **< 10 % of C_j**.
2. **Compute lengths from model (no solver):** run `length_for_f` per line with Q3D loads; verify `length_report` equals targets.
3. **HFSS verification of the same geometry:**
   - **Eigenmode (6 passes, 10 µm mesh):** extract two complex modes per channel → `invert` → f_p, f_r, J, κ_p. (Note: eigenmode κ_p typically runs ~1.5–1.8× higher than driven/Q3D; use driven S21 for J and κ_p.)
   - **Driven S21 (4 passes, explicit length-based mesh, 10 µm default):** fit per channel window using zeros + S21 = 1 point / fixed-zero pole fit (`s21_channel`). Record the mesh length, final ΔS, and any pass-budget deviation in `run.json`.
   - Compare against model predictions without fitting.
4. **Trim if needed:** `ΔL = Δf / s` via `trim_slope` (≈ −1.5 MHz/µm around 6.5–7 GHz), adjusting only the shorted arm (meander absorbs it; hammer and taps stay fixed so J is unchanged).
5. **Loaded check:** include the qubit (junction as lumped L) to verify χ and κ_r at the operating point.

---

## 4. Feed with a series C_in: tee placement and C_p compensation

A series input capacitor C_in reflects the input side of the feed, so each filter's decay interferes with its reflection off C_in. κ_p tracks `Re Z_node` along the feed with period `λ/2 = v / (2f)` (~8.6–10 mm at 6–7 GHz on Si) and can drop near zero at a standing-wave node.

1. **Model the whole feed, not one channel:** include ports, C_in (series + both shunt grounds), every tee with its lead and C_p frame shunt, and all filter–hammer–resonator branches (other channels' C_p shunts shift the standing wave). Use full-feed S21 + pole fit or `κ_p ≈ κ_hanger · Re Z_node / (Z0/2)`.
2. **Avoid minima across robustness corners:** choose tee positions (frequency ordering and feed side are free variables, respecting channel footprint clearances) that maximise worst-case `κ_p / κ_hanger` over perturbations such as `f × (0.99, 1, 1.01)`, `v × (0.98, 1.02)`, `C_frame × (0.8, 1.2)`, and `C_in × (0.9, 1.1)`.
3. **Compensate remaining low-κ_p channels with C_p:** for channels clear of a minimum with ratio < 1, scale `C_p,new = C_p · √(κ_target / κ_model)`, iterate in the full-feed model, extract the new C_p π-network from a Q3D coupon, and recompute only that channel's filter length from the updated C_p loads (resonator length and J stay unchanged).
4. **Record and validate:** record placement, per-channel C_p geometry/loads, and robustness corners in the Design's parameter JSON; validate with full-feed HFSS S21.

---

## 5. Checklist before quoting channel parameters

- [ ] Built with `PurcellReadout`; `length_report` matches target lengths; `resonator_open_extension` declared if used.
- [ ] Q3D runs cited for C_j, head excess, C_p loads, and qubit-loaded resonator end load; non-hammer mutual < 10 % of C_j.
- [ ] Straps at `strap_points()` in every HFSS / Q3D run and air bridges at the same sites on mask.
- [ ] Both dips per channel resolved; poles grouped by channel window, not list order.
- [ ] With series C_in: κ_p from full-feed model, tees clear of minima, C_p compensation recorded.
- [ ] Pass budget, explicit mesh length, final ΔS / ΔF, and any deviation from [em_simulation.md](em_simulation.md) recorded in `run.json`; validation status `unknown` without convergence evidence.
- [ ] Model predictions and HFSS fits labelled separately; no fitted lengths or capacitances.
