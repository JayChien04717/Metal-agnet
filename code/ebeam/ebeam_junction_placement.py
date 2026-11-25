"""Place e-beam junction drawings on the contact hooks of fab QComponents (skill/ebeam_junction.md).

Used inside a Design's ``build_design()`` for the fab render (gds_export.md section 2.1), so the
junction layers come from the same layout code as the base metal; no GDS is read or measured.

    from ebeam_junction_placement import add_ebeam_junctions
    add_ebeam_junctions(design, fluxonium={'FX1': 4.0, 'FX2': 4.0},
                        squid={'TC': (11.528, 11.528)})

* fluxonium: {component name: small-JJ EJ/h in GHz}. The component must provide
  ``contact_hooks()`` (``MyFluxoniumDolanFab``). The drawing is centred between the two hook mouths;
  ``side`` follows the component's array direction; the contact body runs ``overlap`` past the
  hook bottom.
* squid: {component name: (EJ1, EJ2) in GHz}, JJ1 = left branch. The component must provide
  ``contact_hooks()`` (``TmonCouplerDolanNthuFab``). The island hooks must lie 5 um above the ground
  hooks' mouths (TransmonJunctionDolan frame: mouths at +/-2.5 um); ``lower_contact_offset`` is the
  measured ground-hook x offset, and the stem is moved inside the body when needed to keep the
  shadow clearance.

Every small JJ keeps the finger-left orientation, so only components whose hooks are vertical
(stems along y, island above ground for a SQUID) are accepted. Junction sizes come from
code/qubit/jj_resistance_model.py; pass ``model_kw`` to use a recalibrated model. Returns a record
(placements, options, EJ inputs) for the export manifest.
"""
from pathlib import Path
import sys

sys.path[:0] = [str(Path(__file__).resolve().parents[2] / p) for p in
                ('pattern design/ebeam_junction', 'code/qubit')]
from fluxonium_junction_dolan import FluxoniumJunctionDolan  # noqa: E402
from transmon_junction_dolan import TransmonJunctionDolan  # noqa: E402
from transmon_junction_manhattan import TransmonJunctionManhattan  # noqa: E402
import jj_resistance_model as jjm  # noqa: E402

LAYERS_ON_BASE_1 = dict(layer_jj='11', layer_jj_shadow='12', layer_array='13', layer_array_shadow='14')
MANHATTAN_LAYERS_ON_BASE_1 = dict(layer_jj='11', layer_jj_undercut='12')
TOL = 1e-4   # um


def _um(design, value_mm):
    return value_mm / design.parse_value('1um')


def add_ebeam_junctions(design, *, fluxonium=None, squid=None, manhattan=None, layers=None,
                        manhattan_layers=None, margin_um=1.0, overlap_um=1.0, suffix='_EBJ', model_kw=None):
    layers = dict(LAYERS_ON_BASE_1 if layers is None else layers)
    m_layers = dict(MANHATTAN_LAYERS_ON_BASE_1 if manhattan_layers is None else manhattan_layers)
    model_kw = model_kw or {}
    record = dict(layers=layers, manhattan_layers=m_layers, model='code/qubit/jj_resistance_model.py',
                  model_kw=model_kw, ej_inputs_ghz={}, placements={})
    for name, ej in (fluxonium or {}).items():
        info = design.components[name].contact_hooks()
        hooks = {h['name']: h for h in info['hooks']}
        (x1, y1), (x2, y2) = (tuple(_um(design, v) for v in hooks[k]['mouth']) for k in ('top', 'bottom'))
        if abs(x1 - x2) > TOL:
            raise ValueError(f'{name}: hooks are not vertical; the finger-left rule needs stems along y')
        top, bottom = (y1, y2) if y1 > y2 else (y2, y1)
        width = _um(design, hooks['top']['width'])
        depth = abs(_um(design, hooks['top']['end'][1]) - y1)
        cx, cy = x1, (top + bottom) / 2
        ax = info['array_dir'][0]
        if abs(abs(ax) - 1) > 1e-6:
            raise ValueError(f'{name}: array direction is not along x')
        size = jjm.size_from_ej(ej, **model_kw)
        options = dict(pos_x=f'{cx:.4f}um', pos_y=f'{cy:.4f}um', side='left' if ax < 0 else 'right',
                       array_contact_body_width=f'{width + 2 * margin_um:g}um',
                       array_contact_end=f'{(top - cy) + depth + overlap_um:.4f}um',
                       **size.options(), **{k: v for k, v in layers.items()})
        FluxoniumJunctionDolan(design, name + suffix, options=options)
        record['ej_inputs_ghz'][name] = ej
        record['placements'][name + suffix] = dict(component='FluxoniumJunctionDolan', options=options,
                                                   hook_mouths_um=[[x1, top], [x1, bottom]])
    for name, (ej1, ej2) in (squid or {}).items():
        hooks = {h['name']: h for h in design.components[name].contact_hooks()}
        isl = sorted((tuple(_um(design, v) for v in hooks[f'island_{t}']['mouth']), t) for t in ('jj1', 'jj2'))
        gnd = {t: tuple(_um(design, v) for v in hooks[f'ground_{t}']['mouth']) for t in ('jj1', 'jj2')}
        (xl, yl), tl = isl[0]
        (xr, yr), tr = isl[1]
        if abs(yl - yr) > TOL or abs(gnd[tl][1] - gnd[tr][1]) > TOL:
            raise ValueError(f'{name}: SQUID hooks are not level')
        gap = yl - gnd[tl][1]
        if abs(gap - 5.0) > TOL:
            raise ValueError(f'{name}: island-to-ground mouth gap {gap:.4f} um; TransmonJunctionDolan needs 5 um '
                             '(island hook above ground hook)')
        offsets = {gnd[t][0] - x for (x, _), t in isl}
        if max(offsets) - min(offsets) > TOL:
            raise ValueError(f'{name}: the two ground hooks have different x offsets {offsets}')
        offset = offsets.pop()
        width = _um(design, hooks['island_jj1']['width'])
        depth = abs(_um(design, hooks['island_jj1']['end'][1]) - yl)
        if abs(depth - 4.0) > TOL or abs(overlap_um - 1.0) > TOL:
            raise ValueError(f'{name}: TransmonJunctionDolan contacts assume 4 um deep hooks and 1 um overlap '
                             f'(got {depth:.4f} um, {overlap_um} um)')
        # Stem centre relative to the body: centred when it clears the shadow openings, else flush right.
        lead_offset = 0.0 if offset - 1.0 >= 0.338 + 0.648 + 0.075 else 1.0
        sq = jjm.squid_options(ej1, ej2, **model_kw)
        options = dict(pos_x=f'{(xl + xr) / 2:.4f}um', pos_y=f'{yl - 2.5:.4f}um',
                       squid_spacing=f'{xr - xl:.4f}um', upper_contact_body_width=f'{width + 2 * margin_um:g}um',
                       lower_contact_body_width=f'{width + 2 * margin_um:g}um',
                       lower_contact_offset=f'{offset:.4f}um', lower_contact_lead_offset=f'{lead_offset:g}um',
                       layer_jj=layers['layer_jj'], layer_jj_shadow=layers['layer_jj_shadow'], **sq)
        TransmonJunctionDolan(design, name + suffix, options=options)
        record['ej_inputs_ghz'][name] = [ej1, ej2]
        record['placements'][name + suffix] = dict(component='TransmonJunctionDolan', options=options,
                                                   island_mouths_um=[[xl, yl], [xr, yr]],
                                                   ground_offset_um=offset)
    for name, spec in (manhattan or {}).items():
        info = design.components[name].manhattan_sites()
        ori = info['orientation'] % 360.0
        if min(ori, 360.0 - ori) > 1e-6:
            raise ValueError(f'{name}: Manhattan sites have orientation {ori} deg; keep orientation=0 '
                             'so all Manhattan fingers stay along chip x and y')
        cx, cy = (_um(design, v) for v in info['center'])
        options = dict(
            pos_x=f'{cx:.4f}um',
            pos_y=f'{cy:.4f}um',
            orientation='0',
            squid='True' if info['squid'] else 'False',
            layer_jj=m_layers['layer_jj'],
            layer_jj_undercut=m_layers['layer_jj_undercut'],
        )
        if info['squid']:
            options['squid_spacing'] = f"{_um(design, info['squid_spacing']):.4f}um"
        if isinstance(spec, dict):
            options.update(spec)
        TransmonJunctionManhattan(design, name + suffix, options=options)
        record['placements'][name + suffix] = dict(
            component='TransmonJunctionManhattan',
            options=options,
            sites_um=[[_um(design, s['center'][0]), _um(design, s['center'][1])] for s in info['sites']],
        )
    return record
