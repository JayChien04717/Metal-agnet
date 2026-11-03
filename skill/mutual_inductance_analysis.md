# Skill: Analysing the flux-line → qubit-loop mutual inductance

How to get the mutual inductance M between a flux (bias) line and a SQUID / fluxonium loop from an
HFSS driven-modal simulation, check it, and turn it into the numbers you need (bias current per
flux quantum, crosstalk). This file is the procedure; the results for the current chip are in
`note/mutual_inductance_method.md`, the script is `simulation/hfss_mutual_fx.py`.

---

## 1. What M means and what you derive from it

* Flux through the loop for a bias current: `Φ = M · I`.
* Current for one flux quantum: `I_Φ0 = Φ0 / M` (Φ0 = 2.068e-15 Wb; M = 0.5 pH → 4.1 mA).
* Flux-line power / heating: a fluxonium needs ~ ½ Φ0 sweep range, so budget `I_Φ0/2` through the
  50 Ω line, the attenuators and the on-chip short.
* Crosstalk: `M_ij` of line j into loop i; ratio `M_ij / M_ii` is the flux crosstalk.
* Typical target: M ≈ 0.1 – 1 pH (too large → flux noise from the line, too small → mA currents,
  heating).

---

## 2. Method: Z-parameters of a two-port

Two ports:

* **Port 1:** the flux-line pin (lumped port at the pocket wall / launch end, 50 Ω reference).
* **Port 2:** the loop, across the junction gap (the point where the EMF of the loop appears).

With port 2 open (I2 = 0): `V2 = jω M I1` ⇒

```
M = Im(Z21) / (2π f)          # HFSS report: im(Z(2,1))/(2*pi*Freq)
```

### The port-2 trap

* Use the JJ gap as a **port only** (`jj_to_port=[(comp, 'rect_jj', 50, False)]`): the port voltage
  is then the open-loop EMF of the whole loop.
* If the junction inductance L_j is placed in parallel with the port (`True`), the loop closes
  through L_j + L_array and port 2 sees only `EMF · L_j / (L_j + L_array)` ⇒ M is underestimated by
  that ratio (≈ 0.15 for 34 nH / 200 nH).
* With port 2 as port only, no current flows in the loop, so the lumped L values do not change M.
  They matter for eigenmode runs.

### Frequency

Use the quasi-static range, far below any structure resonance (loop L with pad C, ≈ GHz):
discrete sweep 10 – 100 MHz. Too low an f makes the pH signal numerically small, too high picks
up capacitive coupling and resonances.

---

## 3. Setup that converged

* Geometry: simulation version of the qubit (`FluxoniumArrayJJ`), real array size and wire widths,
  pocket, flux line, readout and charge line. Readout / charge pins open-terminated.
* Ground plane: bounding box + 0.2 mm buffer, 525 µm Si.
* Adaptive solution at 100 MHz, max ΔS 0.0002, **≥ 12 passes**, 3 converged passes.
* ΔS alone is a weak criterion at these frequencies: ΔS = 0.002 stopped after 4 passes (7 k
  elements) with M 7 % low. Always compare M between two pass counts (12 vs 15 passes).
* Discrete sweep 10 – 100 MHz, 10 points (~3 min).
* Output: `<tag>_Z.csv` (M, L_fl = Im Z11/ω, L_loop = Im Z22/ω, Re/Im Z21, Z21−Z12), `_summary.json`,
  `_convergence.csv`.
* The JJ port is the one with the larger Im Zii (≈ 200 nH vs ≈ 80 pH). HFSS `GetExcitations`
  order does not match the exported Z-matrix order, so identify ports by Im Zii, not by index.

---

## 4. Sanity checks

| check | expected |
|---|---|
| reciprocity | `Z21 = Z12` to < 1e-6 relative |
| lossless | `Re Z21 ≈ 0` |
| flat plateau | `Im Z21 / ω` flat over 10 – 100 MHz (a few % rise ∝ f² is capacitive flux-line→pad coupling; quote the lowest-f value) |
| L_loop | ≈ `l_inductance + L_geometric` (e.g. 200 nH + ~2 nH) |
| L_fl | tens of pH for a short-circuit U-loop; hundreds of pH for a long one |
| convergence | M changes < 1 % between 12 and 15 passes |
| broadband | `Im Z21` rises with slope 1 (log-log) at low f; pole where L_loop resonates with the pad C; sign flips above it |

Interpolating sweeps can show artifacts (kinks); the discrete sweep is the reference.

### Independent analytical cross-check

Biot–Savart of thin filaments (flux-line centre lines vs the loop area, no ground plane):

```
Φ = ∮∮  μ0 / 4π · (dl1 · dl2) / |r1 − r2|     # Neumann formula, M = Φ / I
```

Agreement within ~10 % confirms the geometry and the port setup (0.44 pH vs 0.44 pH for the thin
filament case here). HFSS can be up to ~40 % above free space when wide PEC pads sit next to a
narrow slot (flux focusing); a mismatch in that direction is plausible, the opposite direction or a
factor ≥ 2 is not.

---

## 5. What drives M (tuning knobs)

* **Loop area** enclosed by the array / junction branch: M ∝ area facing the flux-line current.
* **Flux-line tip distance** to the loop: 20 µm → 15 µm raised M from 0.44 pH to 0.575 pH in this
  design (roughly ∝ 1/d at close range).
* **U-loop geometry:** the two legs (supply + return) cancel at large distance; a shorted U-loop
  with the tip close to the loop is much more effective than a straight line.
* **Return path:** a real chip has the return current in the chip ground and the wirebond; the
  isolated-qubit simulation uses only the local ground. Expect M to change by a few – tens of %
  with full routing; verify with the full-chip run if M is critical.
* **Material model:** PEC sheets ignore kinetic inductance / film penetration of the narrow
  (2 µm) wires. For thin Al lines this slightly reduces the effective loop area; mention it when
  quoting M.

---

## 6. Using the result

```
I_Φ0 = Φ0 / M
```

* Check that the shorted end of the flux line (its narrowest section) can carry `I_Φ0/2` without
  exceeding the film's critical current; compute it from your film thickness / width and the
  measured Jc of your process.
* Dissipation at the fridge stage: P = I² · 50 Ω for AC bias.
* Flux crosstalk between neighbouring qubits: repeat with port 1 on the neighbour's flux line.
* Coupler / SQUID loops: the same method, port 2 across the SQUID junction gap.

---

## 7. Reporting checklist

- [ ] Geometry version (tag), which pins are open, ground-plane size
- [ ] Ports: port 2 as "port only", port identification by Im Zii
- [ ] Passes / elements / ΔS and the M of two different pass counts
- [ ] M at 10 MHz and its spread over the sweep
- [ ] Reciprocity / Re Z21 / L_loop checks passed
- [ ] Biot–Savart cross-check
- [ ] `I_Φ0 = Φ0 / M` and model limits (PEC, local return path)
