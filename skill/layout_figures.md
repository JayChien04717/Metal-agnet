# Skill: Layout figures ("screenshots")

How to make the layout pictures in `Design<N>/fig/` so that every Design looks the same and anyone
can reproduce them. **Layout screenshots are taken with the Qiskit Metal GUI** (`MetalGUI`), driven
by the layout script (rule from Jay, 2026-11-13). The standard set of views of a Design is in its own notes.

---

## 1. Rules

1. Layout screenshots come from the **Qiskit Metal GUI** (`MetalGUI.screenshot`), called from the
   layout script, so they are reproducible. Do not use Snipping Tool, KLayout or Ansys windows, and
   do not use matplotlib plots of the GDS (`plot_gds`) as layout screenshots.
2. Build the same design that is exported (the junction-free mask version, `jj_metal=False`), so the
   screenshot shows what goes on the mask.
3. Each Design has a fixed, small set of views made by its layout script. Designs that share a
   layout script use the same file names and the same views.
4. Screenshots and GDS are regenerated together and committed together.

## 2. How to take a screenshot

```python
from qiskit_metal import MetalGUI
gui = MetalGUI(design)           # design from build_design(jj_metal=False) / build_chip()
gui.rebuild()
gui.autoscale()                  # whole chip
gui.screenshot(name=os.path.join(FIG_DIR, '<name>_chip.png'), display=False)
gui.zoom_on_components(['Q1'])   # zoomed view of one component (list of component names)
gui.screenshot(name=os.path.join(FIG_DIR, '<name>_q1.png'), display=False)
gui.main_window.close()
```

* `screenshot(name, type_='png', display=True)`: `display=False` saves without showing the image inline.
* `zoom_on_components([...])` frames the listed components; `autoscale()` frames the whole chip.
* `highlight_components([...])` / `clear_highlight()` only when a highlight is wanted in the picture.
* The GUI needs a desktop session (it opens a Qt window); run the layout script locally, not headless.
* Keep the same window size for every run so the pictures of different Designs match.

## 3. Adding a new view

1. Choose the components to frame (or `autoscale()` for the chip).
2. Name it `<layout script name>_<what>.png` (lower case, no spaces).
3. Add the `zoom_on_components` + `screenshot` lines to the `__main__` block of the layout script
   (of every Design that shares it).
4. Rerun, compare the PNGs, and add the view to the Design's figure table.

A one-off picture (for a slide or a question) uses the same calls; save it outside the repo:
`fig/` only holds the Design's standard set.

## 4. Checklist before committing figures

- [ ] Taken with the Qiskit Metal GUI by the layout script from the current design.
- [ ] Mask version (junction-free): junction sites empty.
- [ ] `fig/` holds exactly the Design's standard set; README image links point to existing files.

## 5. Other figures (not layout screenshots)

Data plots are made with matplotlib by the scripts that produce the data:

| folder | made by | note |
|---|---|---|
| `fig/mutual/` | `simulation/plot_mutual_z.py` (impedance) | simulation geometry, not the mask |
| `fig/zz/` | `simulation/zz_analysis.py` | ZZ maps, method in the Design's notes |
| fab checks | `plot_gds` (`pattern design/gds_tools.py`) | checking the exported / cheesed GDS (hole pattern, test junctions); not a layout screenshot |
