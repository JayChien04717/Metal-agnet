# Fluxonium Redesign Agenda

Date: 2026-06-01

## Goal

Redesign the local `fluxonium.py` and `fluxonium_module.py` components into a cleaner Qiskit Metal planar fluxonium module for FTF layouts. The module must expose:

- flux bias line
- charge drive line
- readout port
- FTF-facing coupled pad for a transmon coupler

---

## Reference Guidelines & Specs

For detailed design specifications and programming best practices, please refer to:
- 📐 **[Qubit Chip Design Guidelines (Physics)](file:///c:/Users/QEL/Desktop/test%20fluxonium/qubit_chip_design_guidelines.md)**: Planar superconducting qubit design, CPW bends, impedance, rounded corners (fillets) for TLS mitigation, and package edge clearance.
- 🛠️ **[Qiskit Metal Coding Skills (Software)](file:///c:/Users/QEL/Desktop/test%20fluxonium/skill.md)**: QComponent class templates, parameter parsing with units (`self.p`), QGeometry registration (poly/path/junction), and Pin declarations.

---

## Automated Conversation Design Workflow

This agenda is the entry point for automated design conversations. The assistant should read `Agend.md` first to understand the current objective, accepted design rules, known files, and unresolved tasks. When a stage requires specialized rules, the assistant should read the corresponding file in `skill/`.

### Current Skill Files

- `skill/skill.md`: Qiskit Metal component, QGeometry, routing, and notebook coding rules.
- `skill/skill_lom.md`: Ansys Q3D, pyEPR compatibility, Q3D-safe geometry selection, capacitance matrix transforms, EC/coupling calculations, charge-line coupling, and resonator/charge-line loss comparison.
- `qubit_chip_design_guidelines.md`: physical layout, CPW, resonator, package, coherence, and design-rule guidance.

### Entry Protocol

1. Read `Agend.md` at the start of a design or redesign task.
2. Read `qubit_chip_design_guidelines.md` when physical layout, resonator, CPW, package, coupling, or loss questions are involved.
3. Read `skill/skill.md` before editing Qiskit Metal components, routes, notebooks, or geometry generators.
4. Read `skill/skill_lom.md` before running Ansys Q3D, LOM, capacitance matrix transforms, EC/coupling calculations, or charge-line/resonator loss comparisons.
5. Ask for missing target specifications before changing geometry unless the current notebook already records them.

### Required User Inputs

The assistant should collect these through conversation before building or modifying a layout:

- `qubit_number`: number of qubits.
- `qubit_type`: fluxonium, transmon, coupler, or mixed.
- `target_EC_over_h`: target charging energy for each qubit/coupler.
- `target_EJ_over_h` and `target_EL_over_h`: needed for fluxonium parameter/loss estimates when available.
- `target_qubit_frequency`: approximate operating frequency or sweet-spot frequency for each qubit.
- `resonator_frequency`: readout resonator frequency for each qubit.
- `resonator_linewidth`: `kappa/2pi` or `Qext` for each resonator.
- `charge_line_assignment`: which charge/control line drives which qubit/coupler.
- `readout_assignment`: which resonator reads which qubit.
- `loss_priority`: resonator Purcell loss, charge-line loss, flux-line loss, layout area, or a ranked combination.

If values are missing, ask a compact follow-up such as:

```text
Before I build the next layout, I need:
1. qubit_number
2. target EC/h for each qubit
3. resonator frequency and kappa/2pi or Qext
4. approximate qubit frequency for loss comparison
```

### Stage Flow

1. **Specification**
   - Confirm qubit count, target EC, resonator frequencies, resonator linewidths, and control/readout assignments.
   - Deliverable: target-spec table.

2. **Layout Construction**
   - Use or update `FTF.ipynb`.
   - Follow `skill/skill.md` for QComponent and QGeometry rules.
   - Deliverables: layout notebook, preview, geometry sanity checks.

3. **Q3D / LOM Selection**
   - Select Q3D-safe components and open pins.
   - Follow `skill/skill_lom.md`.
   - Exclude geometry that should not enter capacitance extraction: flux-bias pickup paths, external PCB feed routes, and fine superinductor/array path cuts unless explicitly needed.
   - Deliverables: selected component list and open-termination list.

4. **Raw Capacitance Extraction**
   - Run Ansys Q3D through `LOManalysis(design, "q3d")`.
   - Use elevated execution if COM automation is required.
   - Deliverables: raw capacitance matrix and publication-style capacitance table in fF.

5. **Matrix Transform**
   - Define node order and transformation matrix `M`.
   - Compute `C_tilde = M^-1 C M^-1` and `C_tilde_inv = inv(C_tilde)`.
   - Deliverables: transformed matrix and inverse transformed matrix.

6. **Parameter Table**
   - Compute EC, qubit-qubit/coupler coupling, qubit-resonator coupling, charge-line coupling, and MHz/mV drive scales.
   - Deliverables: GHz/MHz parameter table and charge-line coupling table.

7. **Loss Comparison**
   - Compare resonator Purcell loss and charge-line loss.
   - For charge-line loss, use `beta` and the user's effective charge-line impedance/filtering when available.
   - If impedance/filtering is unknown, present a conservative `50 ohm` estimate and label it clearly.
   - Deliverable: `T1_resonator` versus `T1_charge_line` comparison table.

8. **Design Verdict**
   - State whether the layout meets target EC, coupling, drive, and loss priorities.
   - Suggest geometry changes only after comparing current numbers to the target spec.

### Skill Routing Rules

| User intent | Read first | Then use |
|---|---|---|
| design layout, modify component, fix geometry | `Agend.md` | `skill/skill.md` |
| run Ansys, Q3D, capacitance matrix, LOM | `Agend.md` | `skill/skill_lom.md` |
| calculate EC, coupling, transform matrix | `Agend.md` | `skill/skill_lom.md` |
| compare charge-line loss vs resonator loss | `Agend.md` | `skill/skill_lom.md` and `qubit_chip_design_guidelines.md` |
| make a new qubit/coupler class | `Agend.md` | `skill/skill.md` |

### Standard Deliverables

For every complete automated design pass, produce or update:

- Target specification table.
- Layout notebook.
- Geometry preview / sanity check.
- Raw capacitance matrix table in fF.
- Transformed capacitance matrix.
- Inverse transformed matrix.
- EC/coupling table in GHz/MHz.
- Charge-line coupling table, including dimensionless `beta` and MHz/mV.
- Loss comparison table when qubit frequency, resonator linewidth, and charge-line impedance/filtering are known.

### Stop Conditions

Do not claim the design is complete until these are true:

- Required user targets are recorded.
- Notebook builds the layout without geometry errors.
- Q3D-safe selection is documented.
- Ansys Q3D either completes successfully or the exact AEDT/geometry blocker is identified.
- Raw capacitance matrix is produced.
- Matrix transform and MHz/GHz parameter tables are produced.
- Loss comparison is either calculated or clearly marked as waiting for missing physical inputs.

---

## Design Review Notes

- A fluxonium should be represented as a small Josephson junction shunted by a superinductor, normally realized as a Josephson-junction array. The revised core follows the provided device image: two top capacitor pads, a small top junction bridge, and a U-shaped two-leg junction-array superinductor.
- High-coherence and heavy-fluxonium designs usually operate near the half-flux sweet spot, so the layout needs a real loop area and a local inductive flux-bias line close to the loop.
- Readout/control should be separated from the flux-bias function: readout and charge drive are capacitive ports; flux bias is an inductive line placed near the loop.
- FTF layouts couple two fluxonium qubits through a central transmon coupler. The module therefore has a wide side coupler pad with a `transmon_coupler` pin, intended to face the central transmon island.
- Readout topology is waveguide-coupled: a common feedline/waveguide couples to readout resonators, and each resonator branch terminates at the corresponding qubit readout port.
- FTFTFTF rule: `flux_line_options.side='outer'` by default, so the flux loop is automatically placed on the side opposite `transmon_coupler_options.loc_W` to avoid overlap with the transmon-coupler pad/bridge corridor.

---

## Implemented Files

- **`custom/fluxonium.py`**
  - Replaced the old abstract loop with a planar U-array fluxonium core.
  - Added two rounded left/right shunt capacitor pads.
  - Increased the main fluxonium capacitor pad fillet to `12um`, satisfying the qubit-pad rounding guideline (`>=10um`) for TLS reduction.
  - Added a main Josephson junction between the two top nodes.
  - Added two vertical junction-array legs plus a bottom U-turn return.
  - Default `array_mode` is `lumped` for HFSS: only `main_junction` and `array_lumped_inductor` are emitted as junction qgeometry.
  - Optional `array_mode='explicit'` draws 30 junctions per leg (60 total) for visual/layout review.
  - Added a rounded ground pocket cutout.

- **`custom/fluxonium_module.py`**
  - Extends the core component.
  - Adds routable pins: `readout`, `charge`, `flux`, and `transmon_coupler`.
  - Uses capacitive pads for readout/charge.
  - Increased auxiliary pad fillets: readout pad `8um`, charge pad `6um`, and FTF transmon-coupler pad `12um`.
  - Adds `readout_line_options.launch_axis`; default remains `vertical`, while `FTF.ipynb` uses `horizontal` so the readout pin exits outward instead of appearing as a vertical line through the qubit region.
  - Uses an Andersen `fluxonium_pocket.py` style U-shaped flux-bias line near the inductor side.  It enters from the pocket edge, runs beside the array/loop, and exits to an outer-side pin without touching the fluxonium metal.
  - Geometry validation: the flux-bias line is galvanically isolated from the fluxonium pads/array and is routed to the outside of the FTF coupler corridor.
  - Uses a wide FTF-facing transmon coupler pad.

- **`fluxonium.ipynb`**
  - Creates a Qiskit Metal `DesignPlanar`.
  - Instantiates `FluxoniumModule`.
  - Opens `MetalGUI` when available and includes a lightweight matplotlib preview.

- **`resonator_6GHz.ipynb`**
  - Designs a capacitively coupled CPW quarter-wave resonator near 6 GHz on a 10 mm x 10 mm, 500 um thick silicon chip.
  - Target: `kappa / 2pi ~= 1 MHz` (loaded $Q \approx 6000$) with calculated electrical length of ~4.957 mm.
  - Sets initial coupling geometry to `coupling_length=180um` and `coupling_space=6um` to be swept in HFSS.
  - Revised the resonator from a narrow meander to a straight quarter-wave CPW from `CoupledLineTee` to `ShortToGround`; this removes artificial short meander segments that were incompatible with the `99um` fillet rule.  The built length is about `4.9568mm`.

- **`custom/grounded_transmon_coupler.py`**
  - Adds a simplified grounded transmon coupler for FTF cells (left/right pads, central island, grounded stem, one junction, and `drive` pin).
  - Adds a centered coupler flux-bias loop and `flux` pin for tunable-coupler operation.  The loop sits below the grounded junction/stem, matching the intended inductive-bias placement.
  - Moves the coupler drive to an offset capacitive drive pad, since drive coupling does not need to sit directly under the junction.
  - Increased grounded-transmon coupler pad fillet to `12um`.

- **`FTF.ipynb`**
  - Builds a Fluxonium-Transmon-Fluxonium layout in a 10.4 mm x 10.4 mm PCB window with 8 launch pads.
  - Assigns `P1` to the left fluxonium flux-bias line and `P2` to the left fluxonium charge line to prevent flux/charge overlap on the left edge.
  - Reads resonators (quarter-wave, 4.957 mm length, `350um` meander spacing) coupled via `CoupledLineTee`.
  - Lowers the common readout waveguide from `y=3.75mm` to `y=3.35mm` to reduce awkward resonator entry geometry while keeping both resonators at about 4.957 mm.
  - Uses horizontal outward readout launches on both fluxonium modules for the FTF layout.
  - Uses outer-side closed flux pickup loops and outer-side charge pads to route control lines around the perimeter, keeping the central coupler corridor clean.
  - Connects F-T-F coupler bridges (`FL_to_TC_coupler` and `FR_to_TC_coupler`).
  - Employs `ManualCPWRoute` waypoint routes for control lines to avoid diagonal traces, separating P3 (transmon drive) and P4 (transmon-coupler flux) bottom routes with explicit jogs.
  - P3 routes to the offset transmon drive pad; P4 routes to the centered transmon-coupler flux loop below the junction.
  - Current geometry validation reports zero non-endpoint transmission-line crossings.
  - Current previews: `FTF_preview.png` for the full chip and `FTF_center_preview.png` for the central detail.

---

## DRC / Review Notes

- `FTF.ipynb`: zero non-endpoint CPW crossings after the P1/P2 swap and horizontal readout launch update.
- `FTF.ipynb`: active high-coherence geometry (fluxonium modules, coupler, resonators, tee couplers) satisfies the `500um` chip-edge keep-out.  Launch pads and package-bound routes intentionally sit near the chip edge.
- `FTF.ipynb`: pin normals are non-zero and point outward enough for Qiskit Metal routing.
- `FTF.ipynb`: remaining Qiskit warning is local to `LR_resonator` and `RR_resonator`: some automatically generated `RouteMeander` segments are short for `fillet='99um'`.  This is not a crossing, but before HFSS export the readout resonators should be re-laid out with larger horizontal span or a custom waypoint resonator so every bend has enough run-up for the 99 um fillet.
- `resonator_6GHz.ipynb`: revised to a straight resonator to avoid unnecessary meander bends; built length is `4.9568mm`.

---

## References Checked

- Fluxonium overview: https://qubitzoo.org/Zoo/fluxonium
- FTF transmon-coupler architecture: https://arxiv.org/abs/2304.06087
- PRX FTF device overview notes: https://doi.org/10.1103/PhysRevX.13.031035
- Recent fluxonium package language mentioning readout resonators, charge lines, and local flux lines: https://journals.aps.org/prx/pdf/10.1103/PhysRevX.13.031035

---

## Next Calibration Tasks

- Tune shunt pad dimensions to target the desired $E_C$.
- Tune small junction inductance/capacitance to target $E_J$.
- Tune array junction count and per-junction inductance to target $E_L$.
- EM-simulate the readout, charge, flux, and transmon-coupler capacitance/mutual-inductance values.
- Keep the flux line galvanically isolated from the qubit pads/loop.
- Build a two-fluxonium plus central-transmon notebook after the single-module geometry is accepted.
