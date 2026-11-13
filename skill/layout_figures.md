# Skill: Layout figures ("screenshots")

How to make the layout pictures in `Design<N>/fig/` so that every Design looks the same and anyone
can reproduce them. Code: `plot_gds` in `pattern design/gds_tools.py`. The standard figure set and
windows of a Design are in its own notes.

---

## 1. Rules

1. **Never take a manual screenshot** (KLayout, Qiskit Metal GUI, Ansys window, Snipping Tool).
   Every layout figure is drawn by a script from the **exported, junction-free GDS**, so fillets,
   cut-outs and empty junction sites look exactly as on the mask.
2. Use only the shared function `plot_gds`. Do not copy or modify it inside a Design; do not use
   `matplotlib` styles of your own.
3. Each Design has a fixed, small set of figures made by its layout script. Designs that share a
   layout script use the same file names and the same views.
4. Figures and GDS are regenerated together and committed together.

## 2. Format (fixed by `plot_gds`)

| | |
|---|---|
| source | the junction-free GDS without cheesing |
| figure | 10 × 10 in, white background, no title |
| colours | layer 1 metal `#9aa3ad` (grey), etched = white; junction layers ≥ 50 red (should be none) |
| axes | x, y in **mm**, equal aspect |
| dpi | 120 (140 for a junction zoom) |
| labels | PCB port names (`ports=pcb_ports`) on the chip view, if the Design has PCB ports |

`junction_window(design, comp)` centres a view on the junction rows of a component. Call it
**before** `export_gds`, which empties the junction table.

## 3. Adding a new view

1. Decide the window `(x_min, x_max, y_min, y_max)` in mm: read the corners in KLayout (µm / 1000),
   or derive it from `junction_window` / a component's `pos_x`, `pos_y` so it follows the layout.
2. Name it `<layout script name>_<what>.png` (lower case, no spaces).
3. Add the `plot_gds` line to the `__main__` block of the layout script (of every Design that shares it):

   ```python
   plot_gds(gds_file, os.path.join(FIG_DIR, '<name>_<what>.png'), window=(x0, x1, y0, y1))
   ```
4. Rerun, compare the PNGs, and add the view to the Design's figure table.

One-off picture (for a slide or a question), without changing the scripts:

```python
import sys; sys.path.insert(0, r'..\..\pattern design')
from gds_tools import plot_gds
plot_gds(r'..\gds\<file>.gds', r'C:\temp\my_view.png', window=(-0.5, -0.2, -0.95, -0.65))
```
Save it outside the repo: `fig/` only holds the files of section 3.

## 4. Checklist before committing figures

- [ ] Made by the layout script from the current GDS (time stamps of GDS and PNG agree).
- [ ] Junction sites empty; no red in any plot.
- [ ] `fig/` holds exactly the Design's standard set; README image links point to existing files.

## 5. Other figure folders (not layout)

| folder | made by | note |
|---|---|---|
| `fig/mutual/` | `simulation/plot_mutual_z.py` (impedance); `fluxonium_arrayjj*.png` from the `__main__` of the shared `pattern design/fluxonium.py` (`plot_design`, qgeometry tables) | simulation geometry (lumped-L array), not the mask |
| `fig/zz/` | `simulation/zz_analysis.py` | ZZ maps, method in the Design's notes |

