"""Design3's readout transmission line + Purcell filter, without qubits.

    P8 -> interdigital C -> lambda/2 Purcell filter line [meander + 3 coupled-line tees carrying the
          hanger readout resonators, open at the far end, no qubit] -> interdigital C -> P7

The feedline, tees, resonator routing, CPW and fillets are Design3's (ftf_layout.py of
Design3 is imported and its RESONATORS table is re-evaluated). The resonators are
retuned to sit inside the filter pass band (6.6 / 6.7 / 6.8 GHz) and have no qubit pad
(C_pad = 0), so each ends open at the old qubit readout pin.

Usage: python design3_purcell.py   (builds, exports gds/design3_purcell.gds, plots)
"""
import os
import sys
import warnings
from collections import OrderedDict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
D3 = os.path.join(os.path.dirname(ROOT), 'Design3', 'pattern design')
sys.path.insert(0, HERE)
import purcell_layout as pl  # noqa: E402  (Design4: interdigital cap, filter numbers, gds helpers)
sys.path.insert(1, D3)
import ftf_layout as f3  # noqa: E402  (Design3)
from qiskit_metal import designs, Dict  # noqa: E402
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond  # noqa: E402
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee  # noqa: E402
from qiskit_metal.qlibrary.tlines.mixed_path import RouteMixed  # noqa: E402
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander  # noqa: E402
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402
from interdigital_cap import InterdigitalCap  # noqa: E402
from res_end import ResEnd  # noqa: E402
from tap_junction import TapJunction  # noqa: E402
from pcbport import pcb_ports  # noqa: E402

GDS_FILE = os.path.join(ROOT, 'gds', 'design3_purcell.gds')

# readout resonators inside the filter band; no qubit pad
RES_F = OrderedDict([('FX1', 6.6), ('TC', 6.7), ('FX2', 6.8)])        # GHz
# hfss_scale = 1: the analytic lengths are used. (An earlier scale of 0.969 / 0.977 / 0.980 came from HFSS
# eigenmode runs of 10-12 passes that had NOT converged: the mode keeps climbing until ~pass 14 and settles
# at 6.79 GHz for the scaled FX1 resonator, i.e. the analytic length is within ~0.3 %.)
HFSS_SCALE = {'FX1': 1.0, 'TC': 1.0, 'FX2': 1.0}
RESONATORS = OrderedDict()
# Q3D of the 10/6 um tee (q3d_tee_sweep.py, 1.0/1.3/1.6 mm): Cc[fF] = 20.505 * L[mm] - 2.090
# (Design3's 23.985 / -1.082 was for a different CPW). Target kappa = 3 MHz for every resonator.
CC_SLOPE, CC_OFFSET = 20.505, -2.090
KAPPA_TARGET = 3.0   # MHz
for _k, _r in f3.RESONATORS.items():
    r = dict(_r, f=RES_F[_k], C_pad=0.0, hfss_scale=HFSS_SCALE[_k], kappa=KAPPA_TARGET)
    r['Cc'] = f3.coupling_cap(r['f'], r['kappa'])
    r['coupling_length'] = round((r['Cc'] - r['C_extra'] - CC_OFFSET) / CC_SLOPE, 3)
    r['length'] = round(f3.res_length(r['f'], r['C_pad'] + r['Cc']) * r['hfss_scale'], 3)
    RESONATORS[_k] = r

# readout-resonator coupler: 'tee' = Design3 coupled-line tee (long parallel coupling), 'finger' = short T-tap
# on the feedline + interdigital finger capacitor (point-like coupling, as in the literature)
COUPLER = 'finger'
TAP_X = {'FX1': -1.62, 'TC': -0.43, 'FX2': 1.68}                        # mm, feedline positions of the taps
FC = dict(N=6, w=6.0, g=6.0, eg=6.0, bw=20.0, taper=100.0, lead=40.0)   # um, finger-coupler geometry
TAP_PRIME, TAP_STUB, TAP_GAP = 0.30, 0.12, 0.05                         # mm


def finger_overlap(cc_ff):
    """Finger overlap [um] giving cc_ff (Alley estimate scaled by the Q3D calibration)."""
    return round(pl.overlap_for(cc_ff / pl.Q3D_OVER_ALLEY, FC['N']), 1)


def finger_xpin(overlap_um):
    """Centre -> pin distance of the finger coupler [um]."""
    return (overlap_um + 2 * FC['eg']) / 2 + FC['bw'] + FC['taper'] + FC['lead']


def finger_inner(overlap_um):
    """Centre -> resonator-side comb face [um] (the electrical open end of the resonator)."""
    return (overlap_um + 2 * FC['eg']) / 2


# interdigital caps on the feedline height, left / right of the tee chain
CAP_IN_X, CAP_OUT_X = -3.8, 3.8


def qubit_pins():
    """Readout-pin positions / directions of Design3 (read once from its layout)."""
    ref = f3.build_design(for_process=True)
    out = {}
    for n in RESONATORS:
        pin = ref.components[n].pins['readout_line']
        out[n] = (np.array(pin.middle, float), np.array(pin.normal, float))
    return out


def build_design():
    warnings.filterwarnings('ignore')
    qpins = qubit_pins()
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = f3.CPW_WIDTH
    design.variables['cpw_gap'] = f3.CPW_GAP
    design.chips.main.size.size_x = '10mm'
    design.chips.main.size.size_y = '10mm'
    design.chips.main.size.size_z = '-525um'

    for name in ('P8', 'P7'):
        LaunchpadWirebond(design, 'LP' + name, options=dict(f3.LP, **pcb_ports[name]))

    # ---- feedline with three tees (as in Design3) ----
    tee = dict(chip='main', prime_width='cpw_width', prime_gap='cpw_gap',
               second_width='cpw_width', second_gap='cpw_gap', coupling_space='3um',
               down_length=f'{f3.TEE_DOWN}mm', fillet=f3.FILLET, open_termination=True)
    drops = {}
    if COUPLER == 'tee':
        widths = [2 * r['coupling_length'] for r in RESONATORS.values()]
        x = -(sum(widths) + f3.TEE_GAP * (len(widths) - 1)) / 2
        for (name, r), wdt in zip(RESONATORS.items(), widths):
            cx, cl = x + wdt / 2, r['coupling_length']
            mirror = abs(cx - cl / 2) < abs(cx + cl / 2)
            drops[name] = cx - cl / 2 if mirror else cx + cl / 2
            CoupledLineTee(design, f'tee_{name}', options=dict(
                pos_x=f'{cx:.4f}mm', pos_y=f'{f3.FEED_Y}mm', mirror=mirror, coupling_length=f'{cl}mm', **tee))
            x += wdt + f3.TEE_GAP
    else:
        for name, r in RESONATORS.items():
            ov = finger_overlap(r['Cc'])
            r['overlap'] = ov
            cx = TAP_X[name]
            drops[name] = cx
            TapJunction(design, f'tee_{name}', options=dict(
                pos_x=f'{cx:.4f}mm', pos_y=f'{f3.FEED_Y}mm', prime_length=f'{TAP_PRIME}mm',
                stub_length=f'{TAP_STUB}mm', trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP))
            yc = f3.FEED_Y - TAP_STUB - TAP_GAP - finger_xpin(ov) / 1000
            InterdigitalCap(design, f'cc_{name}', options=dict(
                chip='main', pos_x=f'{cx:.4f}mm', pos_y=f'{yc:.5f}mm', orientation='-90',
                finger_count=str(FC['N']), finger_width=f"{FC['w']}um", finger_gap=f"{FC['g']}um",
                overlap=f'{ov}um', end_gap=f"{FC['eg']}um", bus_width=f"{FC['bw']}um",
                taper_length=f"{FC['taper']}um", lead_length=f"{FC['lead']}um",
                trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP, ground_gap='6um'))
            RouteStraight(design, f'stub_{name}', Dict(chip='main', pin_inputs=Dict(
                start_pin=Dict(component=f'tee_{name}', pin='second_end'),
                end_pin=Dict(component=f'cc_{name}', pin='a'))))
    names = list(RESONATORS)
    # ---- Purcell filter = the line between the two interdigital caps; the tees (and so the
    #      readout resonators) sit on it. cap_in -> meander -> tee chain -> cap_out ----
    cap = dict(chip='main', finger_count=str(pl.FINGER_COUNT), finger_width=f'{pl.FINGER_W}um',
               finger_gap=f'{pl.FINGER_G}um', overlap=f'{pl.OVERLAP}um', end_gap=f'{pl.END_GAP}um',
               bus_width=f'{pl.BUS_W}um', taper_length=f'{pl.TAPER_L}um', lead_length=f'{pl.LEAD_L}um',
               trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP, ground_gap='6um', orientation='0')
    InterdigitalCap(design, 'cap_in', options=dict(pos_x=f'{CAP_IN_X}mm', pos_y=f'{f3.FEED_Y}mm', **cap))
    InterdigitalCap(design, 'cap_out', options=dict(pos_x=f'{CAP_OUT_X}mm', pos_y=f'{f3.FEED_Y}mm', **cap))
    f3.route(design, 'feed_in', ('LPP8', 'tie'), ('cap_in', 'a'),
             [(-1.1, 3.9), (CAP_IN_X - 0.6, 3.9), (CAP_IN_X - 0.6, f3.FEED_Y)])
    if COUPLER == 'tee':
        RouteMeander(design, 'filter_in', Dict(
            chip='main', total_length='3mm', fillet=f3.FILLET, hfss_wire_bonds=False,
            lead=Dict(start_straight='150um', end_straight='150um'),
            meander=Dict(spacing='200um', asymmetry='0um'),
            pin_inputs=Dict(start_pin=Dict(component='cap_in', pin='b'),
                            end_pin=Dict(component=f'tee_{names[0]}', pin='prime_start'))))
    else:      # the caps are placed so that the straight line has the filter length
        RouteStraight(design, 'filter_in', Dict(chip='main', pin_inputs=Dict(
            start_pin=Dict(component='cap_in', pin='b'),
            end_pin=Dict(component=f'tee_{names[0]}', pin='prime_start'))))
    for a_, b_ in zip(names[:-1], names[1:]):
        RouteStraight(design, f'feed_{a_}_{b_}', Dict(chip='main', pin_inputs=Dict(
            start_pin=Dict(component=f'tee_{a_}', pin='prime_end'),
            end_pin=Dict(component=f'tee_{b_}', pin='prime_start'))))
    RouteStraight(design, 'feed_out', Dict(chip='main', pin_inputs=Dict(
        start_pin=Dict(component=f'tee_{names[-1]}', pin='prime_end'),
        end_pin=Dict(component='cap_out', pin='a'))))
    f3.route(design, 'out_line', ('cap_out', 'b'), ('LPP7', 'tie'),
             [(CAP_OUT_X + 0.6, f3.FEED_Y), (CAP_OUT_X + 0.6, 3.9), (1.1, 3.9)])
    # total electrical path between the caps' inner pins = pl.PIN_LENGTH (empirical fit)
    fixed = (sum(pl.path_length(design, f'tee_{n}', 'prime_cpw') for n in names)
             + sum(pl.path_length(design, c) for c in ('feed_FX1_TC', 'feed_TC_FX2', 'feed_out')))
    mend = design.components['filter_in']
    for _ in range(10 if COUPLER == 'tee' else 0):
        err = pl.PIN_LENGTH - (fixed + pl.path_length(design, 'filter_in'))
        if abs(err) < 0.002:
            break
        mend.options.total_length = f'{float(mend.options.total_length[:-2]) + err:.4f}mm'
        mend.rebuild()
    design.filter_total = fixed + pl.path_length(design, 'filter_in')

    # ---- hanger resonators: tee -> meander -> open end at the old qubit readout pin ----
    jog_up = {'FX1': 0.30, 'FX2': 0.30, 'TC': 0.75}
    start_comp = {n: ((f'tee_{n}', 'second_end') if COUPLER == 'tee' else (f'cc_{n}', 'b')) for n in RESONATORS}

    def res_len0(name, r):
        if COUPLER == 'tee':
            return r['length'] - r['coupling_length'] - f3.TEE_DOWN
        return r['length'] - (finger_xpin(r['overlap']) - finger_inner(r['overlap'])) / 1000
    for name, r in RESONATORS.items():
        pos, nrm = qpins[name]
        ang = float(np.degrees(np.arctan2(nrm[1], nrm[0])))
        ResEnd(design, f'end_{name}', options=dict(pos_x=f'{pos[0]:.5f}mm', pos_y=f'{pos[1]:.5f}mm',
                                                  orientation=f'{ang:.3f}'))
        y_jog = pos[1] + jog_up[name]
        RouteMixed(design, f'res_{name}', Dict(
            chip='main', total_length=f'{res_len0(name, r):.4f}mm',
            fillet=f3.FILLET, hfss_wire_bonds=False,
            meander=Dict(spacing='230um', asymmetry='0um'),
            lead=Dict(start_straight='200um', end_straight='100um'),
            anchors=OrderedDict({0: np.array([drops[name], y_jog]), 1: np.array([pos[0], y_jog])}),
            between_anchors=OrderedDict({0: 'M', 1: 'S', 2: 'S'}),
            pin_inputs=Dict(start_pin=Dict(component=start_comp[name][0], pin=start_comp[name][1]),
                            end_pin=Dict(component=f'end_{name}', pin='end'))))
        rc = design.components[f'res_{name}']
        tee_len = (f3.path_length(design, f'tee_{name}', 'second_cpw') if COUPLER == 'tee'
                   else (finger_xpin(r['overlap']) - finger_inner(r['overlap'])) / 1000)
        for _ in range(8):
            err = r['length'] - (tee_len + f3.path_length(design, f'res_{name}'))
            if abs(err) < 0.002:
                break
            rc.options.total_length = f'{float(rc.options.total_length[:-2]) + err:.4f}mm'
            rc.rebuild()
    return design


if __name__ == '__main__':
    d = build_design()
    for n, r in RESONATORS.items():
        print(n, 'f', r['f'], 'GHz  length', r['length'], 'mm  tee coupling length', r['coupling_length'], 'mm  Cc', round(r['Cc'], 2), 'fF')
    print('filter line between caps', round(d.filter_total, 4), 'mm (target', pl.PIN_LENGTH, ')')
    pl.export_gds(d, GDS_FILE)
    pl.plot_gds(GDS_FILE, os.path.join(ROOT, 'fig', 'layout', 'design3_purcell_chip.png'))
    pl.plot_gds(GDS_FILE, os.path.join(ROOT, 'fig', 'layout', 'design3_purcell_zoom.png'), window=(-4.6, 4.6, -1.2, 4.7))
