"""Single-channel Design7 comparison with identical feed and hammer metal."""

from argparse import ArgumentParser
import json
from math import pi
from pathlib import Path
import sys

from shapely.geometry import box
from shapely.ops import unary_union

HOME = Path(__file__).resolve().parents[1]
ROOT = HOME.parent
sys.path[:0] = [str(ROOT / path) for path in ('code/general', 'Design4/layout')]
import purcell_6500 as P

SETTINGS = HOME / 'parameters/geometry.json'
VARIANTS = ('parallel', 'opposing')
NET_COMPONENTS = {
    'feed': ('feed_in', 'n0', 'stub_p', 'feed_out'),
    'filter': ('feed_jf', 'jf', 'stub_f', 'filter_tail'),
    'resonator': ('stub_r', 'jr', 'res_open_run', 'res_tail'),
}


def tap(design, name, x, y, **options):
    P.d5.TapJunction(design, name, options=dict(P.CPW, **P.position(x, y), **options))


def end(design, name, x, y, direction):
    P.d5.ResEnd(design, name, options=dict(P.CPW, **P.position(x, y), orientation=str(direction)))


def folded_u_tail(design, name, start, length, fold_x, fold_y, return_x):
    """Route a compact U and set its terminal leg from the filleted length."""
    x, y = design.components[start[0]].pins[start[1]].middle
    radius = design.parse_value(P.d5.f3.FILLET)
    corners = 3
    straight = abs(fold_x - x) + abs(fold_y - y) + abs(return_x - fold_x)
    last = length - straight + corners * radius * (2 - pi / 2)
    if last <= radius:
        raise ValueError(f'{name}: final leg is too short for the fillet')
    stop = name + '_short'
    end(design, stop, return_x, fold_y + last, 270)
    P.d5.f3.route(design, name, start, (stop, 'end'),
                  [(fold_x, y), (fold_x, fold_y), (return_x, fold_y)])


def build(variant, settings=SETTINGS):
    if variant not in VARIANTS:
        raise ValueError(f'Unknown layout {variant!r}')
    cfg = json.loads(settings.read_text(encoding='utf-8'))
    channel = cfg['channel']
    hammer = cfg['hammer_um']
    layout = cfg['variants'][variant]
    length_f = channel['filter_length_mm']
    length_r = channel['resonator_length_mm']
    tap_f, tap_r = channel['filter_tap_mm'], channel['resonator_tap_mm']
    pin_offset = (hammer['gap'] / 2 + hammer['width'] + hammer['lead']) / 1000
    design = P.d5.designs.DesignPlanar()
    design.chips.main.size.update(P.WAFER_CHIP)
    design.chips.main.material = P.d5.SI_MATERIAL
    design.variables.update(cpw_width=P.CPW['trace_width'], cpw_gap=P.CPW['trace_gap'])
    tap(design, 'n0', 0, 0, stub_length='100um')
    P.enclosing(design, 'cp', 0, -P.CP_Y, -90)
    port_x = channel['feed_halfspan_mm']
    end(design, 'input', -port_x, 0, 0)
    end(design, 'output', port_x, 0, 180)
    P.straight(design, 'feed_in', ('input', 'end'), ('n0', 'prime_start'))
    P.straight(design, 'feed_out', ('output', 'end'), ('n0', 'prime_end'))
    P.straight(design, 'stub_p', ('n0', 'second_end'), ('cp', 'a'))

    cp_y = float(design.components['cp'].pins['b'].middle[1])
    half = layout['tee_prime_mm'] / 2
    tee_options = dict(orientation='90', prime_length=f'{2 * half}mm')
    if variant == 'parallel':
        spacing = layout['tee_spacing_mm']
        y_f = y_r = cp_y - tap_f
        x_r = spacing
        cj_y = y_f
        # The 20 um straight link matches the Design6 tap-to-hammer connection.
        stub = spacing / 2 - pin_offset - 0.02
        tap(design, 'jf', 0, y_f, stub_up='False', stub_length=f'{stub:.9f}mm', **tee_options)
        tap(design, 'jr', x_r, y_r, stub_up='True', stub_length=f'{stub:.9f}mm', **tee_options)
        pf, pr = 'second_end', 'second_end'
        short_f, short_r = 'prime_start', 'prime_end'
        open_r, open_sign, orient = 'prime_start', -1, 0
        through_f = through_r = half
    else:
        y_f, y_r = layout['filter_tee_y_mm'], layout['resonator_tee_y_mm']
        x_r, cj_y = 0, layout['hammer_y_mm']
        tap_f = cp_y - y_f
        tap_r = y_r - layout['resonator_open_y_mm']
        stub = layout['tail_tee_stub_mm']
        tap(design, 'jf', 0, y_f, stub_up='False', stub_length=f'{stub}mm', **tee_options)
        tap(design, 'jr', x_r, y_r, stub_up='True', stub_length=f'{stub}mm', **tee_options)
        pf, pr = 'prime_start', 'prime_end'
        short_f = short_r = 'second_end'
        open_r, open_sign, orient = 'prime_start', -1, -90
        through_f = through_r = stub
    P.d5.HammerCap(design, 'cj', options=dict(
        P.CPW, **P.position(x_r / 2, cj_y), orientation=str(orient),
        head_width=f"{hammer['width']}um", head_height=f"{hammer['height']}um",
        head_gap=f"{hammer['gap']}um", ground_gap=f"{hammer['ground']}um",
        lead_length=f"{hammer['lead']}um"))
    P.straight(design, 'stub_f', ('jf', pf), ('cj', 'a'))
    P.straight(design, 'stub_r', ('cj', 'b'), ('jr', pr))
    P.straight(design, 'feed_jf', ('cp', 'b'), ('jf', 'prime_end'))
    P.OpenToGround(design, 'readout_open', options=dict(
        **P.position(x_r, y_r + open_sign * tap_r), orientation='270',
        width=P.CPW['trace_width'], gap=P.CPW['trace_gap'], termination_gap=P.CPW['trace_gap']))
    P.straight(design, 'res_open_run', ('jr', open_r), ('readout_open', 'open'))
    tail_f, tail_r = length_f - tap_f - through_f, length_r - tap_r - through_r
    if variant == 'parallel':
        top = layout['upper_turn_y_mm']
        P.serpentine_tail(design, 'filter_tail', ('jf', short_f), tail_f, layout['filter_lead_mm'],
                          -1, -1, height=top - (y_f - half - layout['filter_lead_mm']))
        P.serpentine_tail(design, 'res_tail', ('jr', short_r), tail_r, top - (y_r + half),
                          1, 1, height=layout['resonator_fold_height_mm'])
    else:
        for name, tee, pin, length, fold in (
            ('filter_tail', 'jf', short_f, tail_f, layout['filter_fold_mm']),
            ('res_tail', 'jr', short_r, tail_r, layout['resonator_fold_mm']),
        ):
            folded_u_tail(design, name, (tee, pin), length, **fold)

    measured = {}
    for label, routes, tee_length, target in (
        ('filter', ('feed_jf', 'filter_tail'), half + through_f, length_f),
        ('resonator', ('res_open_run', 'res_tail'), half + through_r, length_r),
    ):
        actual = sum(P.d5.f3.path_length(design, route) for route in routes) + tee_length
        if abs(actual - target) > 0.001:
            raise ValueError(f'{label}: drawn {actual:.6f} mm differs from {target:.6f} mm')
        measured[label] = actual
    return design, measured


def checks(design, measured):
    """Check connected signal nets and report inter-net separation outside Cj."""
    groups = {name: net for net, names in NET_COMPONENTS.items() for name in names}
    pieces = {net: [] for net in NET_COMPONENTS}
    for kind in ('poly', 'path'):
        for row in design.qgeometry.tables[kind].itertuples():
            geom = row.geometry
            if not geom.is_valid or geom.is_empty:
                raise ValueError(f'Invalid {kind} geometry {row.name}')
            if row.subtract:
                continue
            comp = design._components[row.component].name
            net = groups.get(comp)
            if comp in ('cp', 'cj'):
                side_a = row.name.startswith(('frame', 'head_a', 'lead_a'))
                net = ('feed' if side_a else 'filter') if comp == 'cp' else (
                    'filter' if side_a else 'resonator')
            pieces[net].append(geom.buffer(row.width / 2, cap_style=2) if kind == 'path' else geom)
    nets = {net: unary_union([geom.buffer(0.000001) for geom in parts]) for net, parts in pieces.items()}
    for index, (name, net) in enumerate(nets.items()):
        if net.geom_type != 'Polygon':
            raise ValueError(f'{name} has disconnected metal')
        if any(net.intersects(other) for other in list(nets.values())[index + 1:]):
            raise ValueError(f'{name} has an inter-net short')
    cj = design.components['cj']
    x, y = (design.parse_value(cj.options[key]) for key in ('pos_x', 'pos_y'))
    local = box(x - .30, y - .30, x + .30, y + .30)
    gap = nets['filter'].difference(local).distance(nets['resonator'].difference(local))
    min_x, min_y, max_x, max_y = unary_union([nets['filter'], nets['resonator']]).bounds
    return dict(lengths_mm=measured, filter_res_clearance_outside_hammer_box_um=round(1000 * gap, 3),
                branch_span_mm=round(max_x - min_x, 3), branch_height_mm=round(max_y - min_y, 3),
                hammer_box_um=[600, 600], connected_signal_nets=list(nets))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figures', action='store_true', help='Save Qiskit Metal GUI previews')
    parser.add_argument('--geometry', type=Path, default=SETTINGS)
    parser.add_argument('--variant', choices=VARIANTS)
    args = parser.parse_args()
    available = tuple(json.loads(args.geometry.read_text(encoding='utf-8'))['variants'])
    if args.variant and args.variant not in available:
        parser.error(f'{args.geometry.name} does not define {args.variant}')
    variants = (args.variant,) if args.variant else available
    suffix = '' if args.geometry.resolve() == SETTINGS else '_' + args.geometry.stem
    for variant in variants:
        design, measured = build(variant, args.geometry)
        print(variant, json.dumps(checks(design, measured)), flush=True)
        if args.figures:
            from qiskit_metal import MetalGUI
            gui = MetalGUI(design)
            gui.rebuild()
            gui.autoscale()
            gui.screenshot(name=str(HOME / 'fig' / f'design7_layout_{variant}{suffix}'), display=False)
            if variant == 'opposing':
                gui.zoom_on_components(['jf', 'cj', 'jr'])
                gui.screenshot(name=str(HOME / 'fig' / f'design7_layout_hammer{suffix}'), display=False)
            gui.main_window.close()


if __name__ == '__main__':
    main()
