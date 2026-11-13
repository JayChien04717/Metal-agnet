# Skill: Layout DRC for superconducting qubit chips

Design-rule check for planar qubit chips (transmon, fluxonium, couplers, readout resonators,
feedline) drawn with Qiskit Metal. Part A is the physical reason for each rule, part B is the
rule table to tick off before exporting a GDS, part C is what is specific to the qubit
(junction, pads, flux / charge lines), part D is a pre-tape-out procedure.
Values are for high-resistivity Si (εr ≈ 11.9, 500 µm) and a 10 / 6 µm CPW (W = 10 µm, S = 6 µm,
50.9 Ω); take fab-specific minimum feature sizes from your fab's design rules and use those
where they are stricter.

---

## A. Physics behind the rules

### A.1 Transmission lines
* **No right-angle bends.** Corners give impedance mismatch, local field concentration and loss.
  Use fillets with inner radius `R ≥ 3W` (≥ 5W for high coherence).
* **50 Ω geometry.** W = 10 µm, S = 6 µm on Si (50.9 Ω, 163.4 pF/m, v = 1.2016e8 m/s). Wider
  lines (W = 20, S = 11 µm) reduce edge dielectric loss for high-power feedlines.
* **Tapers.** Between different widths (launch pad → CPW) use a linear taper, angle < 15°.
* **Ground symmetry.** An asymmetric bend lets the two ground planes differ in potential and
  excites parasitic slotline modes. Tie the grounds with airbridges / wire bonds at bends and every
  ~λ/10 (~1 mm at 6 GHz).
* **No crossings** in a single layer: signal-to-signal overlap is a short. Crossings need
  airbridges; otherwise reroute.

### A.2 Qubit loss
* Sharp metal corners concentrate E-field at the metal–substrate / substrate–air interfaces,
  where TLS defects live → lower T1. Round every pad corner.
* Josephson junction (AlOx) is lossy: place it where the qubit's E-field is lowest (in the narrow
  bridge between the pads), not at a pad edge.
* Dicing edges (lattice damage, debris) and the package walls are lossy / couple capacitively.

### A.3 Crosstalk and parasitic modes
* Continuous ground between qubits, resonators and lines; spurious ground islands resonate.
* Meander self-coupling shifts f_r and degrades Q.
* A 10 mm Si chip in a metal box has chip/box modes at several GHz (first mode ≈
  `c / (2 L √εeff)` ≈ 6 GHz for L = 10 mm): check that no qubit or resonator frequency sits on one.

---

## B. Rule table (generic, whole chip)

| # | Category | Check | Spec | Consequence |
|---|---|---|---|---|
| B1 | CPW | Fillet radius at every bend | `R ≥ 3W` (here 99 µm) | reflections, loss |
| B2 | CPW | Shortest straight segment between two bends | `≥ 2 × R` (assert in the build script) | fillet cannot be drawn, kinks |
| B3 | CPW | Cross-section | W / S as design (10 / 6 µm), not changed at tees/launch without taper | impedance step |
| B4 | CPW | Taper angle | `< 15°` | reflection |
| B5 | Routing | Overlaps of different nets | 0 | short |
| B6 | Routing | Route crosses another route / line | none without airbridge | short / box |
| B7 | Routing | Distance between parallel CPW sections | `≥ 5 × (W + 2S)` = 110 µm for 10/6 (meander spacing rule: 150 µm) | self-coupling |
| B8 | Resonator | Total length on the **filleted** centre line | equals target length (iterate `total_length`) | wrong f_r |
| B9 | Ground | Isolated ground islands | none; every island tied to main ground (airbridge / bond / second connection) | slotline modes, spurious resonators |
| B10 | Ground | Ground continuity between lines | wide ground between control lines and qubits | crosstalk |
| B11 | Chip | Keep-out from dicing edge | no active structure within 500 µm | TLS, box coupling |
| B12 | Launch pads | Size / transition | ~200 × 120 µm pads, taper < 15° to CPW | bonding yield, reflection |
| B13 | Ports | Pin normals | point outward; route leaves the component away from the body | route through qubit |
| B14 | Geometry | Self-intersections, zero-area slivers, duplicate points | none (`_clean`) | Ansys/GDS errors |
| B15 | Units | GDS units | µm, chip spans the expected size (here ±5000 µm) | scale error |

---

## C. Qubit-specific DRC

### C.1 Capacitor pads / island
| # | Check | Spec | Why |
|---|---|---|---|
| C1 | Pad corner radius | ≥ 10 µm (10 – 50 µm) on **every** corner incl. inner corners of notches | E-field hot spots → TLS |
| C2 | Pad–ground gap (pocket) | ≥ a few × the narrowest gap of the lines; large enough that Ec/C_Σ matches the Q3D target; not so small that gap-edge field dominates loss | Ec, T1 trade-off |
| C3 | Pad–pad gap at the junction | as designed for the fab (junction site is **empty** in this GDS; the fab adds it) | junction fits |
| C4 | Pad geometry symmetric where assumed | left/right pads match the Ec, charge-matrix assumptions | wrong Ec/asymmetric coupling |
| C5 | Pocket fully contains pad + gap | pocket ⊃ pad by `connection_gap`; subtract geometry present | short to ground |
| C6 | Pad area vs target | Ec target met in the Q3D matrix | frequency |

### C.2 Junction and array
| # | Check | Spec | Why |
|---|---|---|---|
| C7 | Junction placement | in the lowest-field region (narrow bridge between pads) | TLS in barrier |
| C8 | Overlap with pads | leads/pads overlap the area the fab uses for bandage or lead overlap; both ends of the junction line touch metal (simulation) | open circuit / no contact |
| C9 | Junction direction | all junctions of a device drawn parallel (Dolan bridges need the same angle) | shadow evaporation |
| C10 | Test junctions | present in the chip corners (pockets 800 × 300 µm, pads 225/295 × 129 µm, gap 211 µm); not within the cheese region (50 µm keep-out) | process monitoring |
| C11 | Array (fluxonium) | symmetric, straight, minimal loop area; wire width ≥ fab minimum (here 2 µm) | flux noise, yield |
| C12 | Min. feature | smallest wire/gap in the layout ≥ fab limit (list the minimum in the build log) | lithography |
| C13 | Cheese keep-out | no holes within 50 µm of any etched feature (set in the GDS settings) | holes near junctions |

### C.3 Flux and charge lines
| # | Check | Spec | Why |
|---|---|---|---|
| C14 | Flux line vs loop | 100 % galvanic isolation; clearance 10 – 20 µm (≥ 2 × fab tolerance) | defect-induced short |
| C15 | Flux line termination | shorted to ground next to the loop; return path continuous | bias current path |
| C16 | M from simulation | `skill/mutual_inductance_analysis.md`; M in the 0.1 – 1 pH range, `I_Φ0 = Φ0/M` reachable | bias current, heating |
| C17 | Flux-line crosstalk | line of qubit j couples to loop j ≫ loop i | flux crosstalk |
| C18 | Charge line | open end near one pad only; Cc ≈ 0.1 – 1 fF (β small), T1 limit ≥ 10⁷ Q_d | drive-port decay |
| C19 | Charge-line distance to resonator/feedline | not closer than the meander rule (B7) | crosstalk |

### C.4 Resonators and couplers
| # | Check | Spec | Why |
|---|---|---|---|
| C20 | Readout pad to qubit pad | gap sets g; same as in the Q3D model | g |
| C21 | Coupling tee length | from the Cc vs length fit; verify with full-circuit κ run | κ |
| C22 | Resonator spacing between neighbours (frequency collision) | ≥ 100 MHz difference between resonators, qubits far from resonators (dispersive), no qubit on a chip/box mode | collisions |
| C23 | Coupler island | pad sizes consistent with target Ec,TC and J from the matrix | J |
| C24 | Qubit-qubit direct coupling | stays below target J (sum/difference matrix) | unintended coupling |

---

## D. Procedure before exporting

1. **Build:** run the Design's layout script → no warnings about short segments; B2 asserted.
2. **Pins:** print every pin's middle and normal (`design.components[name].pins`), B13.
3. **Plot from the GDS** (not from Metal) so fillets and junction gaps look like the mask; look at
   `fig/layout/*.png` for crossings (B6), islands (B9), corners (C1), the junction sites (C3, C8).
4. **Geometry checks** (script on the GDS / Shapely):
   * min. feature width and gap (B3, C11, C12);
   * min. distance between different nets and from the dicing edge (B7, B11);
   * isolated ground polygons (B9);
   * self-intersections, zero-area polygons (B14);
   * resonator filleted length (B8).
5. **Electrical checks:** Q3D matrix (Ec, J, g, β) – `skill/cap_matrix_analysis.md`; κ; M –
   `skill/mutual_inductance_analysis.md`; HFSS eigenmode for f_r.
6. **GDS:** export with `fab=True` for cheesing (junction-free, see `gds_export.md`), check units
   (B15), no junction cell, test junctions present, edge keep-out of the cheese 0.5 mm.
7. **Log** the result of each check in the commit / `Design<N>/design_parameters.txt`.

---

## E. Common DRC failures seen in this workflow

* A fillet shortens a route, so the resonator ends up too short → iterate on the filleted length.
* Two short segments closer than `2R` → rounding fails or creates a kink.
* Zero-area slivers from `draw.subtract` → Ansys error; drop with `_clean`.
* Grounded flux line drawn as a separate net → Q3D net overlap.
* GDS in mm labelled as µm → chip looks 1000× small in KLayout.
* Junction lines drawn in different directions in one device.
* Hand-drawn polygons that self-intersect → `TopologyException`.

---
*Based on lab design rules and the layout workflow in this repository.*
