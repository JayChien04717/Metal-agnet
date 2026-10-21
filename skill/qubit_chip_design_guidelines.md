# Superconducting Qubit Chip Design Guidelines

This document outlines the core physical and geometric guidelines for designing planar superconducting qubit chips (e.g., Transmon, Fluxonium, resonators, and couplers) using CAD tools such as **Qiskit Metal**.

---

## 1. Coplanar Waveguides (CPW) & Transmission Lines

### 1.1 Fillets & Bends
* **🚫 No Right-Angle Bends**: CPW lines must never feature 90° sharp corners. Sharp bends cause significant impedance mismatches (microwave reflections), local electric field concentrations (leading to dielectric breakdown or loss), and increased insertion loss.
* **📐 The 3W Fillet Rule**: All bends must use smooth arcs (fillets). The inner fillet radius $R$ should satisfy:
  $$R \ge 3W$$
  where $W$ is the width of the central signal trace.
  * *Best Practice*: For high-coherence devices, $R \ge 5W$ is recommended. In Qiskit Metal layouts, standard values like $50\,\mu\text{m}$ or $99\,\mu\text{m}$ are commonly used for a $W = 10\,\mu\text{m}$ line.
* **📈 Tapered Transitions**: When transitioning between different trace widths (e.g., from a large Launch Pad to a narrow CPW line), use a linear taper with an angle of less than $15^\circ$ to ensure a smooth impedance gradient.

### 1.2 Impedance Matching & Dimensions
* **⚡ 50 Ω Characteristic Impedance**: All microwave feedlines and control lines must maintain a strict geometry ratio to achieve a $50\,\Omega$ characteristic impedance. This prevents microwave reflections and standing wave patterns.
* **💎 Silicon Substrate Standard Dimensions**: For high-resistivity silicon substrates ($\epsilon_r \approx 11.9$, thickness $500\,\mu\text{m}$), typical $50\,\Omega$ CPW geometries are:
  * **Signal trace width $W = 10\,\mu\text{m}$**
  * **Gap (spacing) to ground $S = 6\,\mu\text{m}$**
  * *Note*: For high-power feedlines, a scaled geometry like $W = 20\,\mu\text{m}$ and $S = 11\,\mu\text{m}$ is often used to minimize edge dielectric losses.
  * *This FTF design*: W = 15 µm, S = 5.909 µm (45 Ω, 184.5 pF/m, v = 1.2016e8 m/s), 99 µm fillets on all routes.

### 1.3 Overlaps & Crossings
* **🚫 Zero-Overlap Principle**: In single-layer planar designs, signal traces and qubit pads must **never overlap or cross** directly. This leads to an immediate electrical short.
* **🌉 3D Crossovers (Airbridges)**: If a crossing is topologically unavoidable (e.g., a control line crossing a readout feedline), micro-fabricated **airbridges** must be implemented. For strictly planar layouts, crossings must be entirely avoided by rerouting lines around the chip perimeter.

### 1.4 Ground Symmetry & Slotline Mode Suppression
* **⚠️ Ground Plane Symmetry**: Asymmetric paths or bends in a CPW cause a phase difference between the two ground planes, exciting parasitic **slotline modes** that scatter and dissipate energy.
* **🛠️ Strapping/Airbridges**: Place **airbridges** or **wire bonds** at the entrance/exit of every bend and at regular intervals along the transmission line (approximately every $\lambda/10$, which is $\sim 1\,\text{mm}$ at $6\,\text{GHz}$) to equalize the ground potential.

---

## 2. Superconducting Qubits (Transmon / Fluxonium)

### 2.1 Capacitor Pad Geometry
* **⭕ Rounded Corners (Fillets)**: All edges, corners, and internal corners of the qubit capacitor pads (e.g., Transmon shunts, Fluxonium capacitor pads) **must be rounded** with a radius of at least $10\,\mu\text{m} - 50\,\mu\text{m}$.
  * *Physics Mechanism*: Sharp corners concentrate high electric fields, which strongly couple to **Two-Level System (TLS) defects** at the metal-substrate-air (MA/MS/SA) interfaces. Rounding corners is crucial to maximizing the energy relaxation time $T_1$.
* **🕳️ Ground Cutout (Pocket)**: The ground plane cutout around the qubit pads must be carefully simulated to balance the target charging energy $E_C$ against the exposure of the lossy substrate surface.

### 2.2 Josephson Junctions & Superinductor Arrays
* **🧼 Low-Field Placement**: Josephson junctions (the highly lossy $\text{AlO}_x$ barriers) must be placed in regions of minimum electric field strength (e.g., within the narrow bridge corridor connecting the capacitor pads) to mitigate TLS relaxation.
* **⛓️ Superinductor Array (for Fluxonium)**:
  * In Fluxonium qubits, the superinductor is formed by a series Josephson junction array (JJ Array).
  * **Layout Rule**: Keep array paths highly symmetric and straight (e.g., a two-leg U-shape array). Minimize the physical loop area of the array itself to prevent spurious coupling of ambient magnetic flux noise.

### 2.3 Flux Bias Loops
* **🔄 Optimized Loop Area**: The qubit loop area must be carefully balanced. A large area makes the qubit susceptible to environmental dephasing noise, while an excessively small area requires large currents to bias, causing local heating.
* **🔌 Galvanic Isolation**: The local **Flux Bias Line and the qubit loop must be 100% galvanically isolated** (no metal-to-metal connection). They should couple only via spatial mutual inductance ($M \approx 0.1 - 1\,\text{pH}$).

---

## 3. Readout Resonators

### 3.1 Resonator Meanders & Self-Coupling
* **🎨 Large Fillet Radii**: Meandering paths ($\lambda/4$ or $\lambda/2$ CPWs) used to conserve chip area must feature large fillets (typically $99\,\mu\text{m}$) at every turn.
* **↔️ Wide Spacing**: The spacing between adjacent parallel segments of a meander line ($S_{\text{meander}}$) must satisfy:
  $$S_{\text{meander}} \ge 5 \times (W + 2S)$$
  (Typically $\ge 150\,\mu\text{m}$ for $W=10, S=6$). Insufficient spacing causes parasitic capacitive and inductive self-coupling, resulting in unpredictable resonance frequency shifts and degraded quality factors.

### 3.2 Coupling Optimization
* **🔗 Feedline Coupling (Coupled Line Tee)**: The capacitive/inductive coupler between the resonator and the main feedline must be simulated to set the desired external quality factor $Q_e \approx 10^4 - 10^5$, giving a readout bandwidth $\kappa/2\pi \approx 1 - 5\,\text{MHz}$.

---

## 4. Control Lines & Crosstalk Prevention

### 4.1 Control Port Couplings
* **🔋 Charge Drive Lines (Capacitive)**:
  * Terminate with an open-ended CPW trace or a small coupling pad close to one of the qubit pads.
  * The coupling capacitance must be extremely small ($\approx 0.1 - 1\,\text{fF}$), yielding a drive-line quality factor $Q_d \ge 10^7$ to prevent qubit energy relaxation into the drive port (Purcell decay).
* **🧲 Flux Bias Lines (Inductive)**:
  * Terminate with a short-to-ground CPW loop placed near the qubit loop.
  * Keep a safe clearance (e.g., $10 - 20\,\mu\text{m}$) to prevent manufacturing defects from causing a galvanic short.

### 4.2 Shielding & Isolation
* **🛡️ Continuous Ground Planes**: Keep wide, continuous ground planes between control lines, qubits, and resonators to block microwave crosstalk.
* **🗺️ Spatial Channel Routing**: Route drive and readout lines away from each other. For instance, in a Fluxonium-Transmon-Fluxonium (FTF) cell, route qubit controls around the perimeter, and run the coupler flux line perpendicular to the cell to minimize crosstalk.

---

## 5. Chip Boundaries & Packaging

### 5.1 Edge Clearance
* **🚷 Dicing Keep-Out Zone**: Do not place any active qubit pads, resonators, or high-coherence components within $500\,\mu\text{m}$ of the dicing edge.
  * *Reason*: Dicing introduces lattice damage, strain, and micro-debris (TLS defect hotspots). Furthermore, components close to the edge will couple capacitively to the metallic packaging box, causing spurious modes.

### 5.2 Launch Pads (Bond Pads)
* **🛹 Dimensions**: Interface pads near the chip edge should be large (e.g., $200\,\mu\text{m} \times 120\,\mu\text{m}$) for ultrasonic wire bonding.
* **📉 Smooth Transitions**: Taper the interface pad down to standard $50\,\Omega$ CPW dimensions with a transition angle $< 15^\circ$ to minimize microwave reflections.

---

## 6. Design Rule Checking (DRC) Checklist

| Category | Check Item | Standard Spec (Silicon Substrate) | Consequence of Violation |
| :--- | :--- | :--- | :--- |
| **CPW Lines** | Fillet Radius | $R \ge 3W$ (typically $50 - 99\,\mu\text{m}$) | 🚫 Sharp 90° turns cause impedance mismatch and loss. |
| **Routing** | Trace Overlaps | Strictly $0$ Overlaps | 🚫 Signal-to-signal contact causes an electrical short. |
| **Qubit Capacitor** | Corner Rounding | Radius $\ge 10\,\mu\text{m}$ | 🚫 Sharp corners concentrate electric fields, boosting TLS loss and dropping $T_1$. |
| **Qubit Control** | Flux Line Isolation | 100% Galvanic Isolation (No metal contact) | 🚫 Direct contact injects bias current directly into the qubit. |
| **Resonators** | Meander Spacing | $\ge 150\,\mu\text{m}$ ($> 5 \times (W+2S)$) | 🚫 Self-coupling leads to frequency shifts and $Q$ degradation. |
| **Chip Boundary** | Edge Keep-Out | $\ge 500\,\mu\text{m}$ from dicing lines | 🚫 Cutting defects introduce TLS loss; packaging cavity resonances. |
| **Microwave Ports** | Taper Angle | Transition angle $< 15^\circ$ | 🚫 Abrupt changes cause reflections and control pulse distortion. |

---
*Based on experimental best practices and Qiskit Metal modeling guidelines.*
