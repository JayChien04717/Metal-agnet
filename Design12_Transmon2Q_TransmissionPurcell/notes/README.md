# Design12 notes

## 2026-11-26: topology choice (Stage 1, `r1`)

- **User request:** keep the Design11 transmon architecture and replace the individual Purcell filters with one shared, transmission-type filter, using the Design8 IDC style as the reference.
- **First draft (superseded before review):** a λ/2 filter between two equal IDCs (about 91 fF each, `κ ≈ 146 MHz` per port). It was rejected because a symmetric filter sends half of the readout signal back out of the input port (collection efficiency 50 %).
- **Adopted for `r1` (user, 2026-11-26):** a λ/4 filter with a weak input IDC at the open end (`κ_in ≈ 15 MHz`) and a galvanic output tap near the short (`κ_out ≈ 285 MHz`, about 563 µm from the short).
  - **Ports:** input on P7 and output on P8, so the tap T sits outside both resonator coupling sections.
- **Frequency plan (user decision, 2026-11-26):** keep `f_r = 6.35 / 6.15 GHz` (as Design11), so `κ_PF ≈ 300 MHz`. The user accepted the lower Purcell `T_1` estimate (about 0.14 / 0.19 ms, Sete / Jeffrey bound) over pulling `f_r` closer with a narrower filter.

## 2026-11-26: Stage 2 iterations (all HFSS on the readout-only fixture: PF + RR_Q1 + RR_Q2, no qubits)

| Rev | PF model f / coupling Q1, Q2 | Run (result) | Outcome |
|---|---|---|---|
| r1 | 6.25 GHz / 230, 460 um | S21 3 pass `20261126_163208_..._r1_3p_2b2c127c`, fit `20261126_165813_fit_s21_r1_3p_d6c72cc5` | f_PF 5.999, f_r 6.316 / 5.998 GHz, J 35.9 / 33.8 MHz (Design8 scaling 2x low) |
| r2 | 6.772 GHz (6.5 GHz HFSS target via r1 anchor) / 136, 501 um | eigen 6 pass cancelled; eigen 4 pass `..._purcell_eigen_r2_4p_3040b793` (COM error); S21 cancelled | superseded: user meant 6.5 GHz as the DESIGN value |
| r3 | 6.5 GHz design / 122, 244 um | S21 3 pass `20261126_180157_..._r3_3p_20b82629` | f_PF 6.233, f_r 6.325 / 6.276 GHz (Q2 r1 anchor did not carry: its offset came from the 460 um coupling) |
| r4 | 6.5 GHz design / 106, 216 um | S21 3 pass `20261126_183350_..._r4_3p_d21f9aad`; eigen 4 pass crop `20261126_191628_..._fd067d53` read by `20261126_193420_eigen_extract_..._156b4e08` | f_PF 6.240; eigen f_r 6.385 / 6.155 GHz vs driven 6.401 / 6.072 (driven f_r unreliable) |

- Fit fix: HFSS uses e^{+jwt}; the 3-mode model is fitted to conj(S21) (first r1 fits failed for that reason).
- Eigenmode "failures": the COM call waiting on Analyze raised com_error -2147352567 while HFSS kept solving and
  completed (AEDT message log, re-opened session). `simulation/resonator/eigen_extract.py` reads the saved solution.
  A permanent fix in the shared `run_eigen` is still open.
- r4 selected by the user ("use this parameter"); Q1 is about +35 MHz above 6.35 GHz (eigenmode, +-15 MHz).

## 2026-11-26: E_J and coherence

- Best-E_J study (`readout_ej_opt_r4_67ad6c99`): Q1 14.13 (current 14.10), Q2 12.34 (current 11.93) GHz for
  |2chi| = kappa_r. The Q2 change doubles idle ZZ (-150 -> -324 kHz, `cz_verify_r4_bestEJ_84a5dcee`) and lowers its
  Purcell T1 (148 -> 128 us) for an 8 % readout-rate gain. User kept the current E_J.
- Shot-noise T_phi at |2chi| ~ kappa: about 2 / (2 pi n_th kappa); 106 us at n_th = 1e-3.

## 2026-11-26: Stage 3 (NCKU)

- Single die, NCKU (no markers), junction sites empty, no test junctions, die ID prt2q00 (`layout/export_fab.py`).
- First export attempt (`results/general/20261126_200901_gds_export_r4_ncku_c0ae763b`) was stopped: its checks unioned
  the whole cheesed die in shapely (> 18 min). The checks now use point-in-polygon and windowed booleans.
