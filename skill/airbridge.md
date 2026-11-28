# Skill: Airbridges & Octilinear CPW Routing

Design rules for planar reflowed-resist airbridges (signal crossovers and ground-stitching bridges) and octilinear CPW routing on superconducting qubit chips. All numerical examples assume a standard `10 / 6 µm` `50 Ω` CPW (`W = 10 µm`, `S = 6 µm`, total trench width `W + 2S = 22 µm`) on high-resistivity Si; scale proportionally for other CPW dimensions and apply stricter fab limits where required.

---

## 1. CPW Bend Angles & Octilinear Routing

- **No acute (`< 90°`) bends anywhere on the chip.** Every Layer-1 CPW bend must have an interior angle `≥ 90°` (enforced by hard DRC: `interior_deg ≥ 88°` on the final QGeometry `path` table).
- **Prefer `135°` obtuse (octilinear `45°` diagonal) bends** over sharp `90°` bends for fanout buses, diagonal jogs, and corner transitions. Use fillets with `R ≥ 3W` (`50–70 µm` for `10 / 6 µm` CPW).
- **Smooth circular-arc bends (zero sharp unfilleted corners):** Every CPW bend must have a smooth circular arc fillet (never `0 µm` sharp corners), preferably `R ≥ 28–30 µm` for `10 / 6 µm` CPW (`W/2 + S = 11 µm`) so both inner and outer trench edges are visibly rounded. When splitting a route at airbridge pins, clamp the sub-route fillet (`fillet ≤ 0.28 × min_seg`) while keeping `min_seg` long enough (`≥ 70 µm`) to preserve smooth rounded bends without Qiskit Metal short-segment warnings.
- **CPW-to-CPW clearance (hard crosstalk & ground-stitching rule):** Any two distinct CPW centerlines (control–control **and** control–readout/Purcell) **must** stay `≥ 180 µm` apart everywhere (`≥ 158 µm` solid ground between `10/6 µm` trenches; `≥ 250 µm` preferred), so **both** lines have room for `35 µm` ground-stitching airbridges and the ground strip between them is a real, stitched ground. Readout resonator meanders keep `≥ 300 µm` trench-to-trench (see `skill/chip_drc.md` Rule B7).
  - **Definition (no loopholes):** clearance is the **true minimum distance between centerline polylines** (e.g. shapely `distance`), not only the perpendicular gap between exactly parallel segments. Converging diagonals, corners, fanout diagonals, `90°`-crossing jogs and segments next to bridges all count.
  - **Only exemption:** within a disc of radius `d_min / sin θ` around a legitimate bridged crossing (`θ ≥ 85°`), where the two lines necessarily approach each other.
  - **Bridges of different lines:** two crossover bridges (`L20`) of different control lines over the same underpass must be `≥ 180 µm` apart.
  - **Parallel bundles cross together:** when a group of parallel lines crosses a diagonal line, the whole bundle turns together (keeping its pitch); never jog each line independently.
  - **DRC gate equals the rule:** the hard full-chip DRC threshold must equal this rule (`180 µm`); never ship a weaker assert (e.g. `70 µm`) than the written rule.
- **XY charge-line isolation from BOTH couplers AND readout resonators:** Place every XY (charge) drive open stub on an **empty diagonal corner (`≥ 90°` away from the qubit's readout resonator pad)**, pointing radially at `45°` directly toward the circular qubit center (`> 200 µm` away from all 4 adjacent floating coupler pads). Never place an XY open stub on the same corner as the readout resonator, or the microwave XY drive will directly couple into the readout resonator.
- **Monotonic perimeter pad sorting (minimize series crossovers):** Collect **all** control lines exiting each die edge (XY and flux **together**) and sort them monotonically by exit coordinate (`x` for N/S edges, `y` for E/W edges) before assigning perimeter launch pads. This eliminates all same-side criss-crossing, minimizes series crossover airbridges per line (improving fab yield), and prevents multi-crossing pile-ups.
- **Ground-enclosed launch pads:** Inset perimeter launch pads far enough from the die edge (`≥ taper_height + pad_height + back_gap + 140 µm`, e.g., `650 µm` for a `200 µm` taper + `260 µm` pad + `50 µm` back gap) so solid continuous ground plane completely surrounds the outer back edge of every launch port.

---

## 2. Where Airbridges Are Allowed vs. Forbidden

| Structure | Signal Crossover (`L20/L21`) | Ground-Stitching (`L22/L23`) | Reason |
|---|---|---|---|
| Readout resonators & Purcell filter bodies | **Forbidden** (must stay continuous in Layer 1) | **Forbidden** (unless explicitly requested) | Preserves internal quality factor $Q_i$ and avoids frequency/loss shifts |
| Wide qubit pads, coupler pads, `EnclosingCap` pockets | **Forbidden** | **Forbidden** | Long spans (`> 60–100 µm`) sag/collapse and add parasitic capacitance/TLS loss |
| Narrow `10/6 µm` CPW necks & `50 Ω` control lines (XY, flux, external feed) | **Allowed** (overpass only) | **Required** | Clean short spans, low loss, suppresses slotline modes and crosstalk |
| Flip-chip top qubit tier (Q-tier) | **Forbidden** | **Forbidden** | Q-tier holds only qubits, couplers, and ground-stitching micro-bumps; put all airbridges on the C-tier interposer |

- **Priority rule at crossings:** Assign `priority = 0` to Purcell filters and readout resonators, and `priority ≥ 1` to `50 Ω` control lines (XY, flux). When a control line crosses a readout/Purcell line, the `priority = 0` line **always** remains unbroken in Layer 1 and the `priority ≥ 1` line jumps over it.

---

## 3. Signal Crossover Airbridges (`L20` arch / `L21` piers)

### 3.1 Physical Span & Strict 100% `90°` Perpendicular Crossings
- **Hard reflowed-resist arch span limit:** **`25–35 µm` suspended pier-to-pier arch** for **100% of airbridges** on the chip (zero exceptions):
  - Ground-stitching (`L22/L23`): **`25 µm`** suspended arch (`35 µm` total `L22` with `2 × 5 µm` piers).
  - Signal crossover (`L20/L21` over `10 / 6 µm` CPW at `90°`): **`35 µm`** suspended arch (`45 µm` total `L20` with `2 × 5 µm` piers: `22 µm` underpass trench + `2 × 4 µm` intact ground strip + `2 × 2.5 µm` pad end-gap).
- **Mandatory 100% `90°` perpendicular crossings (zero `45°` bridges):** Every single crossover airbridge **must** cross the underlying CPW **perpendicularly at `90°`** (`≥ 85°` hard DRC) so mutual inductance vanishes ($M \propto \cos 90^\circ = 0$), underpass ground strips remain unbroken, and every arch stays strictly `≤ 35 µm`.
  - Whenever a control line approaches a `45°` diagonal line (e.g. a diagonal readout ray or Purcell bus), it must already be travelling perpendicular to it before entering that line's clearance buffer (the §4a router enforces this). Local after-the-fact jogs are a last resort and must be re-checked against the §1 clearance rule.
  - **No triple-junction fallback to `45°`:** If two crossings sit too close on a segment to insert a `90°` jog (e.g., a horizontal and vertical control line meeting on top of a `45°` readout line), **reroute one of the control lines** (`≥ 500 µm` away) rather than falling back to an oblique `45°` bridge.
- **No bends or meanders underneath airbridges:** Readout resonators must exit **straight** radially past all control lines before starting their meanders in the open outer bay (never wiggling or bending underneath control-line airbridges).
- **Never cross over a Purcell–resonator coupling leg:** Control lines crossing a Purcell filter bus must cross in the open single-CPW bus segments between coupling legs — never spanning over both the Purcell bus and a parallel `32 µm`-spaced resonator coupling leg in one oversized bridge.
- **Minimum spacing between consecutive crossover airbridges (`≥ 500 µm` per line):** Along any single CPW line, do not place more than one crossover airbridge (`L20`) within a `500 µm` path-length window. Enforcing `≥ 500 µm` between consecutive `L20` bridges along every route prevents cascaded impedance discontinuities, dense bridge clusters, and accidental multi-crossings of a readout resonator meander U-turn.

### 3.2 Preserving Underpass CPW Impedance (Zero CPW Damage)
- **Flush landing pads:** Keep the overpass landing pad width equal to the CPW trace width (`pad_width = W = 10 µm`, `pier_size = 8 µm`, `bridge_width = 8 µm`) to eliminate capacitive width-step impedance discontinuities.
- **Unbroken underpass ground strip:** Never allow the overpass `top_etch` / `bot_etch` pockets to intrude into the underpass CPW trench or destroy its adjacent ground plane. Clip overpass etch polygons against an explicit underpass protection corridor of half-width:
  $$w_{\text{keepout}} = \frac{W}{2} + S + g_{\text{strip}}$$
  For `10 / 6 µm` CPW with a minimum `g_strip = 6 µm` intact ground strip on both sides, `w_keepout ≥ 17 µm` (`clear_span ≥ 34 µm` at `90°`); widen `g_strip` (e.g. up to `15–20 µm`) when stricter alignment tolerance is specified and arch length permits.
- **Zero silently skipped crossings & max-length cap (Hard DRC):**
  - Never drop real crossings with overly aggressive endpoint or global-proximity skip filters.
  - **Never chain-merge multiple crossings into an oversized mega-bridge.** Only cluster two crossings if their total cluster span stays within the fab arch-length limit; otherwise reroute/spread the lines so each crossing gets its own compact bridge.
  - Always run hard full-chip DRC assertions on the final QGeometry tables verifying **`0` unbridged Layer-1 CPW-to-CPW intersections**, **`0` oblique (`< 85°`) crossover bridges**, AND **`0` Layer-20 airbridges with pier-to-pier suspended arch `> 35 µm`** (`L20 > 45 µm` for `5 µm` piers).

---

## 4. Ground-Stitching Airbridges (`L22` arch / `L23` piers)

- **Purpose:** Tie the two ground planes across every `50 Ω` control CPW (XY drive, flux bias, external feedlines) to suppress parasitic slotline modes and reduce crosstalk.
- **Dimensions:** `span_length = 35 µm` (`25 µm` suspended arch over the `22 µm` CPW trench + `2 × 5 µm` piers on solid ground, `bridge_width = 10 µm`).
- **Spacing & parallel runs:**
  - Place every **`~400–550 µm`** along straight control CPW segments and near bends (`≥ 75 µm` from corners, `≥ 95 µm` from crossover bridges).
  - **Parallel CPW pairs:** Both lines of every parallel control CPW bus **must** receive ground-stitching airbridges to prevent crosstalk. Enforce same-route minimum spacing (`~250 µm`) separately from cross-route non-overlap clearance (`~65 µm`) — never use a large global dedup radius (`> 150 µm`) that suppresses ground bridges on the second parallel line.
  - Verify that both piers land on solid ground plane and never overlap any neighboring CPW trench, qubit/coupler pocket, or crossover etch box.
- **Coverage is a hard DRC, not best effort:** a placement routine must never silently skip a bridge. After placement, every control-line run longer than `550 µm` (measured along the route, excluding `≤ 120 µm` around crossover bridges and component pins) must contain a ground-stitching bridge; otherwise the build fails and the lines are re-spaced.

---

## 4a. Control-Line Routing Algorithm (clearance-aware octilinear maze routing)

Hand-tuned fanout leads, per-line offsets and after-the-fact `90°` jogs do not converge on dense chips (fixing one crossing creates the next). Route every multi-crossing control line with a **clearance-aware octilinear grid router** (shared tool: `code/general/octilinear_router.py`; checker: `code/general/cpw_clearance.py`):

1. **Fixed geometry first:** place components, Purcell filters, readout resonators and fixed short stubs. Rasterize components (pockets + pads + caps) and resonator meander / coupling-leg hulls as **hard obstacles**.
2. **Oriented clearance buffers:** every already-routed CPW contributes a `≥ 180 µm` (+ grid margin) buffer that remembers the line orientation. A new route may enter such a buffer **only straight and perpendicular** to the occupying line (`→` a legitimate `90°` bridged crossing); turning inside a buffer, running parallel, or crossing near a corner (two orientations) is impossible by construction.
3. **Octilinear A\*:** moves along `0° / 45° / 90°` multiples only; `45°` turns (`135°` interior) preferred, `90°` turns penalized, `> 90°` turns forbidden; a minimum straight run between bends leaves room for `R ≥ 50 µm` fillets. Costs: length + bends + per-crossing penalty + soft penalty inside the `250 µm` preferred zone.
4. **Pins and pads:** each net starts with a fixed straight stub along the pin normal and ends with a straight approach along the launch-pad axis; pending nets' stubs are reserved before routing so earlier nets cannot block them. Grid-snapped end segments are re-intersected exactly onto the pin / pad axes.
5. **Order and pad assignment:** keep monotonic perimeter pad sorting; route innermost nets first; rip up and reroute a net that fails.
6. **Bridges after routing:** insert one `≤ 35 µm`-arch `L20` bridge at each (now perpendicular) crossing, then ground-stitching bridges, then run the hard DRC: `0` oblique crossings, `0` unbridged crossings, `0` clearance violations (`< 180 µm`), `0` oversized `L20`, ground-bridge coverage, and `≥ 500 µm` between consecutive `L20` bridges along a line (report any exceptions explicitly).

---

## 5. GDS Layering & Visual Verification

- **Separate GDS layers** so crossover and ground airbridges can be inspected, modified, or selectively omitted independently:
  - `L1`: Base metal (Nb / Al ground plane, CPWs, pads)
  - `L20` / `L21`: Signal crossover airbridge arch / piers
  - `L22` / `L23`: Ground-stitching airbridge arch / piers
  - `L30`: Flip-chip Indium micro-bump UBM pads (when applicable)
- **Color-coded figures:** Always render signal crossover airbridges (`L20`) and ground-stitching airbridges (`L22`) in distinct high-contrast colors (e.g. orange vs. purple, plus green for Indium bumps) with a legend on chip and unit-cell figures.
- **KLayout / close-up check:** Before approving Stage 1, zoom in on representative `90°` and `45°` crossings in the exported GDS to verify intact underpass ground strips, flush `50 Ω` pads, `≥ 90°` bends, and `0` unbridged crossings.
