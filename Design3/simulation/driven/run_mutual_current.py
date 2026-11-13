"""Reproduce the Design3 low-frequency flux-line to loop mutual-inductance runs.

Run in the qmetal environment with AEDT open:
    python run_mutual_current.py tc-left tc-right fx1 fx2 tc-hex
"""
import argparse
import math
from pathlib import Path
import sys

from shapely.geometry import LineString, Point

DESIGN3 = Path(__file__).resolve().parents[2]
ROOT = DESIGN3.parent
sys.path[:0] = [
    str(DESIGN3 / 'pattern design'),
    str(ROOT / 'pattern design'),
    str(ROOT / 'Design5' / 'simulation'),
]
from ftf_layout import build_design
from hfss_mutual_fx import run_design
import save_aedt

TC_CASES = {
    'tc-left': ('d3_tc_rect20x14_m_left', 'rect_jj1', False),
    'tc-right': ('d3_tc_rect20x14_m_right', 'rect_jj2', False),
    'tc-hex': ('d3_tc_hex16x24_m_left', 'rect_jj1', True),
}


def build_tc_fixture(open_jj, legacy):
    """Replace the other junction with a PEC sheet for a single open loop port."""
    options = dict(squid_loop_shape='hexagon', squid_loop_width='16um',
                   squid_loop_height='24um', squid_contact_overlap='0um') if legacy else None
    design = build_design(marks=False, tmon_options=options)
    tc = design.components['TC']
    island, ground, hole, *_ = tc._place(list(tc._squid_geometry()))
    assert math.isclose(hole.area, (288 if legacy else 280)*1e-6, abs_tol=1e-12)

    junctions = design.qgeometry.tables['junction']
    return_jj = junctions.loc[junctions.component.eq(tc.id) & junctions.name.ne(open_jj)]
    assert len(return_jj) == 1
    row, = return_jj.itertuples()
    a, b = [Point(coord) for coord in row.geometry.coords]
    assert a.distance(island) < 1e-12 and b.distance(ground) < 1e-12
    ux, uy = (b.x-a.x)/row.geometry.length, (b.y-a.y)/row.geometry.length
    land = .00025  # mm overlap of the simulation-only closure onto each electrode
    wire = LineString([(a.x-land*ux, a.y-land*uy), (b.x+land*ux, b.y+land*uy)])
    bridge = wire.buffer(row.width/2, cap_style=2)
    assert bridge.intersection(island).area > 0 and bridge.intersection(ground).area > 0
    tc.add_qgeometry('poly', {'mutual_loop_closure': bridge})
    design.qgeometry.tables['junction'] = junctions.drop(index=return_jj.index)
    return design, row.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', nargs='+', choices=(*TC_CASES, 'fx1', 'fx2'))
    args = parser.parse_args()
    save_aedt.use_project(DESIGN3 / 'simulation/aedt/Design3.aedt')
    for case in args.case:
        info = {'ground': 'local box plus buffer'}
        if case in TC_CASES:
            tag, junction, legacy = TC_CASES[case]
            component, flux_pin = 'TC', 'flux_bias_line'
            design, closure = build_tc_fixture(junction, legacy)
            info.update(open_jj=junction, pec_return_jj=closure,
                        opening_um=[16, 24] if legacy else [20, 14])
        else:
            tag, component, flux_pin, junction = 'd3_current_m_'+case, case.upper(), 'flux_line', 'rect_jj'
            design = build_design(for_process=True, marks=False, fx_options={'flux_as_ground': 'False'})
            q = design.components[component]
            polygons = design.qgeometry.tables['poly']
            assert 'array_end' in set(polygons.loc[polygons.component.eq(q.id), 'name'])
            info.update(flux_height_um=80, flux_tip_gap_um=11, stem_width_um=10,
                        array_model='two 100 nH arms and metal end')
        base_tag = tag
        suffix = 2
        results = DESIGN3 / 'simulation/results/runs'
        project = save_aedt._desktop().get_active_project()
        existing_designs = set(project.get_design_names())
        while (results / f'{tag}_summary.json').exists() or f'{tag}_hfss' in existing_designs:
            tag = f'{base_tag}_r{suffix}'
            suffix += 1
        run_design(tag, design, component, flux_pin, junction,
                   [(component, 'readout_line'), (component, 'charge_line')], info=info)
        save_aedt._desktop().get_active_project()._project.Save()


if __name__ == '__main__':
    main()

