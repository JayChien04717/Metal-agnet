# Skill: Layout figures ("screenshots") for Design1 – Design3

How to make the layout pictures in `Design<N>/fig/layout/` so that all Designs look exactly the
same and anyone can reproduce them. Code: `plot_gds` in `SQL HFSS/pattern design/gds_tools.py`.

---

## 1. Rules

1. **Never take a manual screenshot** (KLayout, Qiskit Metal GUI, Ansys window, Snipping Tool).
   Every layout figure is drawn by a script from the **exported, junction-free GDS**, so fillets,
   cut-outs and empty junction sites look exactly as on the mask.
2. Use only the shared function `plot_gds`. Do not copy or modify it inside a Design; do not use
   `matplotlib` styles of your own.
3. Only the four figures of section 3. Same picture in every Design = same file name. A new view
   is added to the `__main__` block of **all three** `Design<N>/pattern design/ftf_layout.py`.
4. Figures and GDS are regenerated together and committed together.

## 2. Format (fixed by `plot_gds`)

| | |
|---|---|
| source | `gds/ftf_layout.gds` (junction-free, no cheesing) |
| figure | 10 × 10 in, white background, no title |
| colours | layer 1 metal `#9aa3ad` (grey), etched = white; junction layers ≥ 50 red (should be none) |
| axes | x, y in **mm**, equal aspect |
| dpi | 120 (140 for the junction zoom) |
| labels | PCB port names (`ports=pcb_ports`, chip view of Design1 / 3; Design2 has no ports) |

## 3. Standard figure set

Run once per Design (Python of the `qmetal` env, see `skill.md` §9):

```
cd "D:\...\SQL HFSS\Design1\pattern design"        # Design2, Design3 the same
& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I ftf_layout.py
```

This writes `gds/ftf_layout.gds`, `gds/ftf_layout_fab.gds` and the four figures below. Windows are
`(x_min, x_max, y_min, y_max)` in mm (Design1 / 3: `Y0 = -0.8`).

| file (`fig/layout/`) | shows | Design1 / Design3 window | Design2 window |
|---|---|---|---|
| `ftf_layout_chip.png` | whole chip | (−5, 5, −5, 5) + port labels | (−5, 5, −5, 5) |
| `ftf_layout_zoom.png` | the FTF set | (−1.0, 1.0, Y0 − 0.55, Y0 + 0.55) | (−3.2, −1.5, −2.9, −0.7) |
| `ftf_layout_fx1.png` | FX1 / Q1 with flux, charge, readout lines | (−0.95, −0.3, Y0 − 0.3, Y0 + 0.3) | (−2.85, −1.75, −2.82, −1.72) |
| `ftf_layout_fxjj.png` | FX1 / Q1 junction site (must be empty) | (−0.48, −0.32, Y0 − 0.045, Y0 + 0.045), dpi 140 | `junction_window(d, 'Q1', half=0.12, half_y=0.045)`, dpi 140 |

`junction_window(design, comp)` centres the view on the junction rows of a component. Call it
**before** `export_gds`, which empties the junction table.

No other layout pictures are made (no SQUID, die-ID or cheesed-fab plots); check those in
KLayout if needed, without saving pictures.

## 4. Adding a new view (only if it is needed in every Design)

1. Decide the window in mm. Open `gds/ftf_layout.gds` in KLayout and read the corner coordinates
   (KLayout shows µm: divide by 1000), or use `junction_window` / a component's `pos_x`, `pos_y`
   so the window follows the layout when it moves.
2. Pick a name `ftf_layout_<what>.png` (lower case, no spaces).
3. Add the line to the `__main__` block of all three `ftf_layout.py` (matching qubit in each
   Design, e.g. FX1 ↔ Q1, TC ↔ Q3):

   ```python
   plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_<what>.png'), window=(x0, x1, y0, y1))
   ```
4. Rerun `ftf_layout.py` in each Design, compare the three PNGs side by side, add the file to the
   table in section 3.

One-off picture (for a slide or a question), without changing the scripts:

```python
import sys; sys.path.insert(0, r'..\..\pattern design')
from gds_tools import plot_gds
plot_gds(r'..\gds\ftf_layout.gds', r'C:\temp\my_view.png', window=(-0.5, -0.2, -0.95, -0.65))
```
Save it outside the repo: `fig/` only holds the files of section 3.

## 5. Checklist before committing figures

- [ ] Made by `ftf_layout.py` from the current GDS (time stamps of GDS and PNG agree).
- [ ] Junction site empty in `ftf_layout_fxjj.png`; no red in any plot.
- [ ] `Design1/fig/layout`, `Design2/fig/layout`, `Design3/fig/layout` hold exactly the four files.
- [ ] README image links point to existing files (`Design1/3`: `ftf_layout_chip.png`, `Design2`: `ftf_layout_zoom.png`).

## 6. Other figure folders (not layout)

| folder | made by | note |
|---|---|---|
| `fig/mutual/` | `simulation/plot_mutual_z.py` (impedance); `fluxonium_arrayjj*.png` from the `__main__` of the shared `pattern design/fluxonium.py` (`plot_design`, qgeometry tables) | simulation geometry (lumped-L array), not the mask |
| `fig/zz/` | `simulation/zz_analysis.py` | ZZ maps, `../note/zz_method.md` |
