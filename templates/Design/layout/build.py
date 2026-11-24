# -*- coding: utf-8 -*-
r"""Template single layout entry point (see ``skill/layout_structure.md``).

Reads ``placement.json`` and ``routes.json`` directly. Register any shared or
Design-specific ``QComponent``s in ``COMPONENT_REGISTRY``; never import
simulation modules or mutate/delete components after building.
"""
from collections import OrderedDict
import copy
import json
from pathlib import Path
import sys

_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'code' / 'general' / 'project_paths.py').is_file())
sys.path.insert(0, str(_ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)

DESIGN_DIR = project_paths.design_dir(__file__)
LAYOUT_DIR = Path(__file__).resolve().parent
COMPONENTS_DIR = LAYOUT_DIR / 'components'
if COMPONENTS_DIR.is_dir() and str(COMPONENTS_DIR) not in sys.path:
    sys.path.insert(0, str(COMPONENTS_DIR))

import numpy as np
from qiskit_metal import Dict, designs
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.mixed_path import RouteMixed
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

from bondpad import BondPad, add_bond_pad
from check_keepouts import check_keepouts, collect_keepouts
from junction_test_coupon import HookedTestJunction
from nthu_tiau_marker import TiAuChipMarker
from route_length_solver import path_length, solve_route_length
from tsri_a12_mark import TSRIA12Mark
from xmark import XMark

COMPONENT_REGISTRY = {
    'BondPad': BondPad,
    'XMark': XMark,
    'TSRIA12Mark': TSRIA12Mark,
    'TiAuChipMarker': TiAuChipMarker,
    'HookedTestJunction': HookedTestJunction,
}


def _deep_merge(base, patch):
    out = copy.deepcopy(base)
    for k, v in (patch or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _to_metal_dict(obj):
    if isinstance(obj, dict):
        return Dict({k: _to_metal_dict(v) for k, v in obj.items()})
    return obj


def _resolve_coord(design, val):
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, dict) and 'pin' in val:
        comp_name, pin_name = val['pin']
        axis = 0 if val['axis'] == 'x' else 1
        return float(design.components[comp_name].pins[pin_name].middle[axis]) + float(val.get('offset', 0.0))
    raise ValueError(f'Unsupported anchor coordinate specification: {val!r}')


def build_design(variant='sim', overrides=None, placement_path=None, routes_path=None,
                 check_keepout_um=100.0):
    """Build the chip design from ``placement.json`` and ``routes.json``."""
    if variant not in ('sim', 'process', 'fab'):
        raise ValueError(f"variant must be 'sim', 'process', or 'fab'; got {variant!r}")

    placement = json.loads(Path(placement_path or (LAYOUT_DIR / 'placement.json')).read_text(encoding='utf-8'))
    routes_doc = json.loads(Path(routes_path or (LAYOUT_DIR / 'routes.json')).read_text(encoding='utf-8'))
    overrides = overrides or {}
    comp_overrides = overrides.get('components', {})
    route_overrides = overrides.get('routes', {})

    chip = placement['chip']
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = chip['material']
    design.variables['cpw_width'] = chip['cpw_width']
    design.variables['cpw_gap'] = chip['cpw_gap']
    design.chips.main.size.size_x = chip['size_x']
    design.chips.main.size.size_y = chip['size_y']
    design.chips.main.size.size_z = chip['size_z']

    for name, spec in placement['components'].items():
        allowed = spec.get('variants')
        if allowed is not None and variant not in allowed:
            continue
        if 'pcb_anchor' in spec:
            add_bond_pad(design, name, _deep_merge(spec['pcb_anchor'], comp_overrides.get(name)))
            continue
        cls_name = spec.get('variant_classes', {}).get(variant, spec['class'])
        opts = _deep_merge(spec.get('options', {}), spec.get(f'{variant}_options'))
        opts = _deep_merge(opts, comp_overrides.get(name))
        COMPONENT_REGISTRY[cls_name](design, name, options=_to_metal_dict(opts))

    length_report = {}
    default_fillet = routes_doc.get('fillet', '99um')
    for name, spec in routes_doc['routes'].items():
        cfg = _deep_merge(spec, route_overrides.get(name))
        pin_inputs = Dict(
            start_pin=Dict(component=cfg['start'][0], pin=cfg['start'][1]),
            end_pin=Dict(component=cfg['end'][0], pin=cfg['end'][1]),
        )
        rclass = cfg['class']
        if rclass == 'RouteStraight':
            RouteStraight(design, name, Dict(chip='main', pin_inputs=pin_inputs))
        elif rclass == 'RouteAnchors':
            anchors = OrderedDict(
                (i, np.array([_resolve_coord(design, a[0]), _resolve_coord(design, a[1])], dtype=float))
                for i, a in enumerate(cfg.get('anchors', []))
            )
            RouteAnchors(design, name, Dict(
                chip='main', fillet=cfg.get('fillet', default_fillet),
                lead=_to_metal_dict(cfg.get('lead', {'start_straight': '0um', 'end_straight': '0um'})),
                anchors=anchors, pin_inputs=pin_inputs,
            ))
        elif rclass == 'RouteMixed':
            anchors = OrderedDict(
                (i, np.array([_resolve_coord(design, a[0]), _resolve_coord(design, a[1])], dtype=float))
                for i, a in enumerate(cfg.get('anchors', []))
            )
            between = OrderedDict((int(k), v) for k, v in cfg['between_anchors'].items())
            init_len = cfg.get('initial_total_length_mm', cfg.get('target_length_mm'))
            RouteMixed(design, name, Dict(
                chip='main', total_length=f'{float(init_len):.4f}mm',
                fillet=cfg.get('fillet', default_fillet), hfss_wire_bonds=False,
                meander=_to_metal_dict(cfg.get('meander', {})),
                lead=_to_metal_dict(cfg.get('lead', {})),
                anchors=anchors, between_anchors=between, pin_inputs=pin_inputs,
            ))
            if 'target_length_mm' in cfg:
                extra_fn = None
                if 'extra_length_path' in cfg:
                    ec, en = cfg['extra_length_path']
                    extra_fn = lambda d, c=ec, n=en: path_length(d, c, n)
                length_report[name] = solve_route_length(
                    design, name, float(cfg['target_length_mm']),
                    extra_length_fn=extra_fn,
                    tol_mm=float(cfg.get('tol_mm', 0.001)),
                    check_segments=bool(cfg.get('check_segments', True)),
                )
        else:
            raise ValueError(f'Unsupported route class {rclass!r} for route {name!r}')

    design.metadata['layout_revision'] = placement.get('revision', 'unknown')
    design.metadata['variant'] = variant
    design.metadata['route_lengths'] = length_report
    if check_keepout_um is not None and collect_keepouts(design):
        design.metadata['keepout_check'] = check_keepouts(design, min_um=float(check_keepout_um))
    return design
