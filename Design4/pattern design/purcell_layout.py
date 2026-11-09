"""Design4: transmission Purcell filter at 6.7 GHz on a 10 x 10 mm chip.

    P1 --- feed --- [interdigital C] --- lambda/2 meander --- [interdigital C] --- feed --- P6

A lambda/2 CPW resonator (open at both ends) is coupled to the input and output line through
two identical interdigital finger capacitors, so the device is a band-pass in transmission
(f0 = 6.7 GHz, bandwidth kappa_F). Launch pads sit at the PCB ports P1 / P6 (pcbport.py).
All CPW: 10 um centre / 6 um gap (50.9 Ohm), 99 um fillets.

Method (note/purcell_filter_method.md):
    Q_F = f0 / kappa_F;  each port gives Q_c = 2 Q_F;  Cc = sqrt(pi / (Q_c Z0 Zr)) / w
    interdigital C from Alley's formula (first estimate, to be replaced by Q3D);
    lambda/2 length with both ends loaded by Cc.

Usage:  python purcell_layout.py     (builds, exports gds/, plots fig/layout/)
"""
import os
import sys
from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GDS_DIR = os.path.join(ROOT, 'gds')
FIG_DIR = os.path.join(ROOT, 'fig', 'layout')
sys.path.insert(0, HERE)
from pcbport import pcb_ports
from interdigital_cap import InterdigitalCap

# ---- specification ----------------------------------------------------------
F0 = 6.7                # GHz, filter centre frequency
KAPPA_F = 200.0         # MHz, filter bandwidth (Q_F = 33.5; typical Purcell filter Q_F ~ 30)
CPW_WIDTH, CPW_GAP = '10 um', '6 um'
CPW_Z0 = 50.9           # Ohm
CPW_C = 163.4e-12       # F/m
CPW_V = 1.2016e8        # m/s
EPS_R = 11.45           # silicon
FILLET = '99um'

# ---- interdigital capacitor (N fingers, 6 um fingers / gaps) ----------------
FINGER_COUNT = 10
FINGER_W, FINGER_G, END_GAP = 6.0, 6.0, 6.0       # um
BUS_W, LEAD_L = 20.0, 40.0                        # um
TAPER_L = 100.0                                   # um, sloped taper CPW -> comb height (impedance match)

# ---- chip placement (mm) ----------------------------------------------------
Y_LINE = 1.6            # P1 / P6 height
CAP_X = 1.2             # interdigital caps at x = -CAP_X, +CAP_X
LP = dict(chip='main', lead_length='30um', pad_width='200um', pad_height='200um', pad_gap='76.5um')


def alley_cap_fF(N, overlap_um, w_um=FINGER_W, h_um=525.0, eps_r=EPS_R):
    """Interdigital capacitance [fF], Alley / Bahl formula (substrate thickness h):
    C = (eps_r + 1) l [(N - 3) A1 + A2],  A1 = 4.409 tanh(0.55 (h/W)^0.45) e-6 pF/um,
    A2 = 9.92 tanh(0.52 (h/W)^0.5) e-6 pF/um,  W = finger width."""
    A1 = 4.409e-6 * np.tanh(0.55 * (h_um / w_um) ** 0.45)
    A2 = 9.92e-6 * np.tanh(0.52 * (h_um / w_um) ** 0.5)
    return (eps_r + 1) * overlap_um * ((N - 3) * A1 + A2) * 1e3


def overlap_for(C_fF, N=FINGER_COUNT):
    """Finger overlap [um] giving C_fF (inverse of alley_cap_fF, linear in l)."""
    return C_fF / alley_cap_fF(N, 1.0)


def coupling_cap(f_ghz, Q_c, Z0=CPW_Z0, Zr=CPW_Z0):
    """Cc [fF] for an external Q_c (lab formula): Q_c = pi / (Z0 Zr (Cc w)^2)."""
    w = 2 * np.pi * f_ghz * 1e9
    return np.sqrt(np.pi / (Q_c * Z0 * Zr)) / w * 1e15


def res_length(f_ghz, C_load_ff):
    """lambda/2 length [mm], total end loading C_load: f = v / 2L (1 - C_load / C_line)."""
    f = f_ghz * 1e9
    L = CPW_V / (2 * f)
    for _ in range(50):
        L = CPW_V / (2 * f) * (1 - C_load_ff * 1e-15 / (CPW_C * L))
    return L * 1e3


# Calibration from the first runs (note/purcell_filter_method.md):
#  - HFSS S21 (s21_v1): 3 dB bandwidth 390 MHz where the lab formula (Q_c = 2 Q_F per port) predicted
#    188 MHz from the Q3D Cc -> the formula gives kappa a factor KAPPA_CORR too small for this
#    two-port filter, so the formula is evaluated for kappa_F / KAPPA_CORR.
#  - Q3D Cc (97.99 fF) / Alley (101.0 fF) = Q3D_OVER_ALLEY.
KAPPA_CORR = 390.0 / 188.1
Q3D_OVER_ALLEY = 97.99 / 101.0
Q_F = F0 * 1e3 / KAPPA_F
CC = coupling_cap(F0, 2 * Q_F * KAPPA_CORR)                    # fF per port (Q3D value wanted)
OVERLAP = round(overlap_for(CC / Q3D_OVER_ALLEY), 1)           # um

CAP_FF = alley_cap_fF(FINGER_COUNT, OVERLAP) * Q3D_OVER_ALLEY
# filter length, measured between the two cap pins (the caps' own length adds a bit of line):
CAP_LEN = (OVERLAP + 2 * END_GAP + 2 * BUS_W + 2 * (LEAD_L + TAPER_L)) / 1000          # mm, whole cap
FILTER_LENGTH = round(res_length(F0, 2 * CAP_FF) * (1 + 0.0), 3)           # mm; hfss_scale later
HFSS_SCALE = 1.0        # (unused; replaced by the empirical fit below)
# Empirical electrical length from two HFSS runs (s21_v1: path 7.206 mm, Cc 98 fF, f0 6.48 GHz;
# s21_v2: path 7.462 mm, Cc 70 fF, f0 6.45 GHz): v / (2 f0) = L_path + L0 + GAMMA * Cc,
# L_path = filleted centre line between the two cap pins (caps' own metal / taper included in L0)
L0_MM, GAMMA_MM_FF = 1.320, 0.00761
PIN_LENGTH = round(CPW_V / (2 * F0 * 1e9) * 1e3 - L0_MM - GAMMA_MM_FF * CAP_FF, 4)   # mm


def filleted_length(coords, r):
    """Length of a polyline after rounding every corner with radius r."""
    c = np.asarray(coords, float)
    L = np.linalg.norm(np.diff(c, axis=0), axis=1).sum()
    for p0, p1, p2 in zip(c[:-2], c[1:-1], c[2:]):
        a, b = p1 - p0, p2 - p1
        t = np.arccos(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))
        L -= 2 * r * np.tan(t / 2) - r * t
    return L


def path_length(design, comp, name=None):
    """Filleted centre-line length [mm] of a route in the design."""
    t = design.qgeometry.tables['path']
    rows = t[(t['component'] == design.components[comp].id) & (~t['subtract'])]
    if name:
        rows = rows[rows['name'] == name]
    r = rows.iloc[0]
    fil = r['fillet'] if r['fillet'] == r['fillet'] and r['fillet'] else 0.0
    return filleted_length(r.geometry.coords, float(fil))


def build_design(overlap_um=None, filter_length=None, finger_count=None):
    ov = OVERLAP if overlap_um is None else overlap_um
    L_pin = PIN_LENGTH if filter_length is None else filter_length
    N = FINGER_COUNT if finger_count is None else finger_count
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    design.chips.main.size.size_x = '10mm'
    design.chips.main.size.size_y = '10mm'
    design.chips.main.size.size_z = '-525um'

    for name in ('P1', 'P6'):
        LaunchpadWirebond(design, 'LP' + name, options=dict(LP, **pcb_ports[name]))

    cap = dict(chip='main', finger_count=str(N), finger_width=f'{FINGER_W}um',
               finger_gap=f'{FINGER_G}um', overlap=f'{ov}um', end_gap=f'{END_GAP}um',
               bus_width=f'{BUS_W}um', taper_length=f'{TAPER_L}um', lead_length=f'{LEAD_L}um',
               trace_width=CPW_WIDTH, trace_gap=CPW_GAP, ground_gap='6um')
    InterdigitalCap(design, 'cap_in', options=dict(pos_x=f'{-CAP_X}mm', pos_y=f'{Y_LINE}mm', **cap))
    InterdigitalCap(design, 'cap_out', options=dict(pos_x=f'{CAP_X}mm', pos_y=f'{Y_LINE}mm', **cap))

    # input / output feed lines
    for nm, lp, c, pin in (('feed_in', 'LPP1', 'cap_in', 'a'), ('feed_out', 'LPP6', 'cap_out', 'b')):
        RouteStraight(design, nm, Dict(chip='main', fillet=FILLET, pin_inputs=Dict(
            start_pin=Dict(component=lp, pin='tie'), end_pin=Dict(component=c, pin=pin))))

    # lambda/2 filter between the inner pins of the two caps
    RouteMeander(design, 'filter', Dict(
        chip='main', total_length=f'{L_pin:.4f}mm',
        fillet=FILLET, hfss_wire_bonds=False,
        lead=Dict(start_straight='150um', end_straight='150um'),
        meander=Dict(spacing='200um', asymmetry='0um'),
        pin_inputs=Dict(start_pin=Dict(component='cap_in', pin='b'),
                        end_pin=Dict(component='cap_out', pin='a'))))
    # fillets shorten the line: iterate on the filleted centre-line length. The two caps
    # contribute their own metal (bus + fingers) to the resonator, which is not counted here.
    res = design.components['filter']
    for _ in range(10):
        err = L_pin - path_length(design, 'filter')
        if abs(err) < 0.002:
            break
        res.options.total_length = f'{float(res.options.total_length[:-2]) + err:.4f}mm'
        res.rebuild()
    return design


def plot_gds(filename, png, window=(-5, 5, -5, 5), dpi=120):
    """Plot layer 1 (metal) of an exported GDS."""
    import gdspy
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    lib = gdspy.GdsLibrary(infile=filename, unit=1e-3, units='convert')
    top = [c for c in lib.top_level() if c.name == 'TOP'][0]
    fig, ax = plt.subplots(figsize=(10, 10 * (window[3] - window[2]) / (window[1] - window[0])))
    ax.set_facecolor('white')
    for (layer, _), ps in sorted(top.get_polygons(by_spec=True).items()):
        if layer == 1:
            ax.add_collection(PolyCollection(ps, facecolor='#9aa3ad', edgecolor='none'))
    for name, opt in pcb_ports.items():
        x, y = float(opt['pos_x'][:-2]), float(opt['pos_y'][:-2])
        if window[0] < x < window[1] and window[2] < y < window[3]:
            ax.annotate(name, (x, y), fontsize=11, color='darkred', weight='bold', ha='center', va='center')
    ax.set_xlim(window[0], window[1]); ax.set_ylim(window[2], window[3])
    ax.set_aspect('equal'); ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    fig.tight_layout(); fig.savefig(png, dpi=dpi); plt.close(fig)


def _prune_gds(filename):
    """Rewrite the GDS in um, keeping only cells referenced from TOP."""
    import gdspy
    lib = gdspy.GdsLibrary(infile=filename, unit=1e-6, units='convert')
    top = lib.cells['TOP']
    keep, stack = {top.name: top}, [top]
    while stack:
        cell = stack.pop()
        cell.references = [r for r in cell.references if not isinstance(r.ref_cell, str)]
        for ref in cell.references:
            if ref.ref_cell.name not in keep:
                keep[ref.ref_cell.name] = ref.ref_cell
                stack.append(ref.ref_cell)
    out = gdspy.GdsLibrary(name=lib.name, unit=lib.unit, precision=lib.precision)
    out.add(list(keep.values()))
    out.write_gds(filename)


def export_gds(design, filename):
    design.renderers.gds.options['no_cheese']['view_in_file'] = {'main': {1: False}}
    design.renderers.gds.options['cheese']['view_in_file'] = {'main': {1: False}}
    design.renderers.gds.export_to_gds(filename)
    _prune_gds(filename)


def summary():
    print(f'f0 {F0} GHz, kappa_F {KAPPA_F} MHz -> Q_F {Q_F:.1f}, Q_c per port {2 * Q_F:.1f}')
    print(f'Cc per port {CC:.1f} fF -> N={FINGER_COUNT} fingers, overlap {OVERLAP} um (Alley {CAP_FF:.1f} fF)')
    print(f'analytic lambda/2 length {FILTER_LENGTH} mm; empirical pin-to-pin path {PIN_LENGTH} mm')


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    summary()
    d = build_design()
    for c in ('cap_in', 'cap_out', 'filter'):
        print(c, {k: tuple(np.round(v.middle, 4)) + tuple(np.round(v.normal, 2)) for k, v in d.components[c].pins.items()})
    print('filter filleted length', round(path_length(d, 'filter'), 4), 'mm')
    gds_file = os.path.join(GDS_DIR, 'purcell_filter.gds')
    export_gds(d, gds_file)
    plot_gds(gds_file, os.path.join(FIG_DIR, 'purcell_chip.png'))
    plot_gds(gds_file, os.path.join(FIG_DIR, 'purcell_zoom.png'), window=(-2.2, 2.2, Y_LINE - 1.5, Y_LINE + 1.5))
    plot_gds(gds_file, os.path.join(FIG_DIR, 'purcell_cap.png'), window=(-CAP_X - 0.18, -CAP_X + 0.18, Y_LINE - 0.1, Y_LINE + 0.1), dpi=140)
