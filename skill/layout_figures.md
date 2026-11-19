# Skill: Layout figures

How the pictures in `Design<N>/fig/` are made, so that every Design has the same structure and style and anyone can
reproduce them.

Rule change (2026-11-19, from Jay):
- Layout figures are rendered from the current GDS by one shared script.
- They replace the Qiskit Metal GUI screenshots of the 2026-11-13 rule.
- The GDS is what goes on the mask, the render runs headless, and every Design looks the same.

---

## 1. What `fig/` holds

Every Design has the same structure and nothing else:

| File | Content |
|---|---|
| `chip.png` | Whole chip. Each part window is boxed and labelled with its file name. |
| `<part>.png` | One close-up per part: `q1`, `q2`, `coupler`, `readout`, `purcell_resonator_pair`, `hammer`, `cin`, ... |
| `README.md` | Generated. Lists each figure with its revision, source GDS (sha256) and window. |
| `.gitignore` | Generated. Only the files above are tracked. |

- **Several chip variants** (Design8, DesignTestJunction): `chip_<variant>.png` and `<variant>_<part>.png`.
- **Not allowed in `fig/`:** design-ID checks, before/after comparisons, fab-mask previews, dated copies and other
  one-off pictures. A script may still write them there for a local check: the generated `.gitignore` keeps them
  untracked.
- **Superseded figures** go to the ignored `fig/previous/`.

## 2. Making the figures

1. The windows are in [`code/general/design_figures.json`](../code/general/design_figures.json). For each Design
   (and variant) the file gives:
   - the source GDS and revision;
   - an optional chip window;
   - one window per part, `[x0, x1, y0, y1]` in µm in the GDS frame.
2. Use the **clean companion GDS** of the current revision: the simulation GDS or the unnotched reference, with no
   cheese holes. That way the parts read clearly.
3. Render from the repository root in qmetal:

   ```powershell
   python -X utf8 -B code/general/plot_design_figs.py            # all Designs
   python -X utf8 -B code/general/plot_design_figs.py Design6    # one Design
   ```

4. **Style** (fixed in the script):
   - grey = metal (all layer-1 datatypes), white = substrate, red = layer 50/0 junction annotations;
   - a scale bar on every figure;
   - the title gives the Design, revision and part.
5. Very narrow windows are widened to an aspect of at least 0.75, so every close-up reads well. The README records
   the window actually drawn.

## 3. After a geometry change

1. Export the new GDS (see [gds_export.md](gds_export.md)).
2. Point the Design's entry in `design_figures.json` at the new GDS and revision, and adjust any window that moved.
3. Re-render the Design and check every PNG. Commit the GDS and figures together, so they match the same revision.

**Adding a part:** add a window to `design_figures.json` and re-render. The part's file name is its key.

## 4. Figures that are not layout figures

Data plots and method figures belong to the note that cites them, in `Design<N>/notes/fig/`:

| Folder | Made by | Note |
|---|---|---|
| `Design1/notes/fig/mutual/` | `Design1/simulation/flux/plot_mutual_z.py` | Mutual-inductance method |
| `Design1/notes/fig/zz/` | `Design1/simulation/qubit/zz_analysis.py` | ZZ maps |
| `Design2/3/notes/fig/junction_parameters.png` | `code/general/plot_junction_layout.py` | Junction model overlay |
| `Design6/notes/fig/` | `Design6/simulation/qubit/plot_ej_layout.py` and the 11-15 report | Report and S21 evidence |

A junction parameter figure is a model overlay on the GDS:
- Annotate the junction and array sites from the QGeometry table.
- State the source of L_J, E_J/h and I_c.
- Name the HFSS values and any independent circuit or fabrication targets explicitly.
- Use `code/general/plot_junction_layout.py` for Design2 and Design3, and
  `Design6/simulation/qubit/plot_ej_layout.py` for the HFSS-retargeted Design6 junctions.

## 5. Checklist before committing figures

- [ ] Rendered by `plot_design_figs.py` from the GDS named in `design_figures.json` (the current revision).
- [ ] `fig/` holds only `chip*.png`, the part PNGs, `README.md` and `.gitignore`.
- [ ] The part windows show the right parts. Check each PNG after a layout change.
- [ ] README and note image links point to existing files.
