"""Junction-free GDS export and GDS plots, shared by Design1 - Design3 (skill/gds_export.md).

    from gds_tools import export_gds, plot_gds
    export_gds(build_design(jj_metal=False), 'gds/ftf_layout.gds')             # no cheesing
    export_gds(build_design(jj_metal=False), 'gds/ftf_layout_fab.gds', fab=True)  # cheesed

All Designs use the same plot function, so their figures look identical: layer 1 grey on
white, square 10 x 10 in figure, axes in mm (skill/layout_figures.md).
"""
import os


# Lab fabrication / cheesing settings for the GDS renderer (skill/gds_export.md)
# design units are mm: 0.005 = 5 um holes, 0.05 = 50 um pitch
def apply_fab_gds_options(gds_render):
    gds_render.options.short_segments_to_not_fillet = False
    gds_render.options.tolerance = 1e-6
    gds_render.options.fabricate = True
    gds_render.options.negative_mask = dict(main=[2])
    gds_render.options.no_cheese['view_in_file']['main'] = {1: False, 2: False}
    gds_render.options.no_cheese['buffer'] = '50um'
    gds_render.options.cheese['edge_nocheese'] = '0.5mm'
    gds_render.options.cheese['view_in_file']['main'] = {1: True, 2: False}
    gds_render.options.cheese.update(dict(cheese_0_x=0.005, cheese_0_y=0.005, delta_x=0.05, delta_y=0.05))


def export_gds(design, filename, fab=False):
    """Write the junction-free GDS (routes get their real fillets here).

    Build the design with jj_metal=False: then the pads carry no junction-related metal.
    The junction table (JJ, array arms, SQUID junctions) is dropped here, so no junction
    cell is placed and every junction site is empty; the fab draws the junctions.
    fab=True applies the lab cheesing settings (apply_fab_gds_options).
    Note: this edits `design` in place; build a new design for the next file."""
    design.qgeometry.tables['junction'] = design.qgeometry.tables['junction'].iloc[0:0]
    gds = design.renderers.gds
    if fab:
        apply_fab_gds_options(gds)
    else:
        gds.options['no_cheese']['view_in_file'] = {'main': {1: False}}
        gds.options['cheese']['view_in_file'] = {'main': {1: False}}
    gds.export_to_gds(filename)
    _prune_gds(filename)


def _prune_gds(filename):
    """Convert to um and drop cells not referenced from TOP, and references to cells that
    no longer exist (qiskit's fabricate mode deletes its intermediate Cheese_diff cell but
    leaves the reference to it)."""
    import gdspy
    # qiskit writes mm (unit 1e-3); convert to um (unit 1e-6) like the lab's GDS files.
    # (read_gds' default units='skip' would keep the mm numbers but label them um)
    lib = gdspy.GdsLibrary(infile=filename, unit=1e-6, units='convert')
    top = lib.cells['TOP']
    keep = {top.name: top}
    stack = [top]
    while stack:
        cell = stack.pop()
        cell.references = [r for r in cell.references if not isinstance(r.ref_cell, str)]
        for ref in cell.references:
            c = ref.ref_cell
            if c.name not in keep:
                keep[c.name] = c
                stack.append(c)
    out = gdspy.GdsLibrary(name=lib.name, unit=lib.unit, precision=lib.precision)
    out.add(list(keep.values()))
    out.write_gds(filename)


def plot_gds(filename, png, window=(-5, 5, -5, 5), dpi=120, ports=None):
    """Plot layer 1 (metal, grey) and any junction layers (red) of an exported GDS, in mm.
    ports: {name: options with pos_x / pos_y} labels (e.g. pcb_ports), optional."""
    import gdspy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    lib = gdspy.GdsLibrary(infile=filename, unit=1e-3, units='convert')     # plot in mm
    top = [c for c in lib.top_level() if c.name == 'TOP'][0]
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor('white')
    for (layer, _), ps in sorted(top.get_polygons(by_spec=True).items()):
        if layer == 1:
            ax.add_collection(PolyCollection(ps, facecolor='#9aa3ad', edgecolor='none'))
        elif layer >= 50:
            ax.add_collection(PolyCollection(ps, facecolor='red', edgecolor='none'))
    for name, opt in (ports or {}).items():
        x, y = float(opt['pos_x'][:-2]), float(opt['pos_y'][:-2])
        if window[0] < x < window[1] and window[2] < y < window[3]:
            ax.annotate(name, (x, y), fontsize=11, color='darkred', weight='bold', ha='center', va='center')
    ax.set_xlim(window[0], window[1]); ax.set_ylim(window[2], window[3])
    ax.set_aspect('equal'); ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    fig.tight_layout(); fig.savefig(png, dpi=dpi); plt.close(fig)


def junction_window(design, comp, half=0.06, half_y=None):
    """Plot window (mm) centred on the junction site(s) of component `comp` (centre of its
    junction-table rows). Call before export_gds, which empties the junction table."""
    jt = design.qgeometry.tables['junction']
    rows = jt[jt['component'] == design.components[comp].id]
    x0, y0, x1, y1 = rows.geometry.total_bounds
    x, y = (x0 + x1) / 2, (y0 + y1) / 2
    hy = half if half_y is None else half_y
    return (x - half, x + half, y - hy, y + hy)
