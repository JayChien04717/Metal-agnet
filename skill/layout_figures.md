# Skill: Layout figures

How the pictures in `Design<N>/fig/` are generated so every Design has a uniform, reproducible visual record rendered directly from its current GDS.

---

## 1. What `fig/` holds

| File | Content |
|---|---|
| `chip.png` | Whole chip with each part window boxed and labelled with its filename |
| `<part>.png` | One close-up per part (`q1`, `q2`, `coupler`, `readout`, `purcell_resonator_pair`, `hammer`, `cin`, …) |
| `README.md` | Generated index listing each figure's revision, source GDS (SHA256), and window |
| `.gitignore` | Generated; tracks only the standard files above |

- **Multiple chip variants:** use `chip_<variant>.png` and `<variant>_<part>.png`.
- **Not tracked in `fig/`:** design-ID checks, before/after diffs, fab-mask previews, and dated one-off pictures (the generated `.gitignore` keeps local checks untracked).
- **Superseded figures:** move to the ignored `fig/previous/`.

---

## 2. Making the figures

1. Define windows in [`code/general/design_figures.json`](../code/general/design_figures.json): source GDS, revision, optional chip window, and per-part windows `[x0, x1, y0, y1]` in µm in the GDS frame.
2. Point to the **clean companion GDS** of the current revision (`*_sim.gds` or unnotched reference, no cheese holes) so parts read clearly.
3. Render from the repository root in `qmetal`:

   ```powershell
   python -X utf8 -B code/general/plot_design_figs.py              # all Designs
   python -X utf8 -B code/general/plot_design_figs.py <Design>     # one Design
   ```

4. **Fixed script style:** grey = metal (all layer-1 datatypes), white = substrate, red = layer `50/0` junction annotations, scale bar on every figure, title showing Design / revision / part. Very narrow windows are widened to aspect ≥ 0.75 (recorded in `fig/README.md`).

---

## 3. After a geometry change

1. Export the new GDS ([gds_export.md](gds_export.md)).
2. Update the Design's entry in `code/general/design_figures.json` (source GDS, revision, and any moved windows).
3. Re-render, inspect every PNG, and commit GDS and figures together so they correspond to the same revision.

**Adding a view:** add a window key to `design_figures.json` and re-render.

---

## 4. Non-layout figures

Data plots, simulation-geometry overlays, and method figures belong in `Design<N>/notes/fig/` alongside the note that cites them.

For a junction-parameter overlay on a GDS (`code/general/plot_junction_layout.py` or a Design-specific script):
- Annotate junction and array sites from the QGeometry table.
- State the source of L_J, E_J/h, and I_c explicitly, separating HFSS values from independent circuit or fabrication targets.

---

## 5. Checklist before committing figures

- [ ] Rendered by `code/general/plot_design_figs.py` from the current-revision GDS in `design_figures.json`.
- [ ] `fig/` tracks only `chip*.png`, part PNGs, `README.md`, and `.gitignore`.
- [ ] Every part window frames the intended feature after the latest layout change.
- [ ] README and note image links resolve to existing files.
