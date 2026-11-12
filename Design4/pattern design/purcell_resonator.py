"""Reusable Design4 Purcell filter and readout resonator without hammer electrodes."""

import warnings

import numpy as np
from qiskit_metal.qlibrary.core import QComponent
from shapely.affinity import rotate, translate
from shapely.geometry import LineString

from purcell_6500 import build_no_hammer

RECOMMENDED_CPW_SPACING_UM = 250


class PurcellResonator(QComponent):
    """A fixed 10/6 µm CPW pattern, with two feedline pins and HFSS wirebonds.

    Use on a 500 µm silicon chip for comparison with the spacing study.
    The CPW spacing option is measured between the two straight centerlines.
    """

    component_metadata = dict(
        short_name='purcell_res',
        _qgeometry_table_poly='True',
        _qgeometry_table_path='True',
    )
    default_options = dict(
        cpw_spacing='400um',
        filter_length='4.191257715mm',
        resonator_length='4.505973012mm',
        resonator_trim='0um',
        ground_wire_bonds=True,
    )

    def make(self):
        p = self.p
        spacing_um = p.cpw_spacing * 1000
        if spacing_um < RECOMMENDED_CPW_SPACING_UM:
            warnings.warn(
                f'{self.name}: {spacing_um:g} µm CPW centerline spacing is below the '
                f'{RECOMMENDED_CPW_SPACING_UM} µm reference for reduced parasitic coupling. '
                'Verify the minimum mode splitting for this layout.',
                UserWarning,
                stacklevel=2,
            )
        source = build_no_hammer(
            (p.filter_length, p.resonator_length),
            spacing_mm=p.cpw_spacing,
            res_trim=p.resonator_trim,
        )
        excluded = {'LPPin', 'LPPout'}
        names = {component.id: name for name, component in source.components.items()}

        def place(geometry):
            return translate(rotate(geometry, p.orientation, origin=(0, 0)),
                             xoff=p.pos_x, yoff=p.pos_y)

        for kind in ('poly', 'path'):
            for row in source.qgeometry.tables[kind].itertuples():
                name = names[row.component]
                if name in excluded:
                    continue
                options = {}
                if kind == 'path':
                    options.update(width=row.width, fillet=row.fillet,
                                   hfss_wire_bonds=bool(row.subtract and p.ground_wire_bonds))
                self.add_qgeometry(kind, {f'{name}_{row.name}': place(row.geometry)},
                                   subtract=bool(row.subtract), helper=bool(row.helper),
                                   **options)

        for name, component, endpoint in (
                ('feed_in', 'feed_in', 'start'), ('feed_out', 'feed_out', 'end')):
            pin = source.components[component].pins[endpoint]
            points = np.asarray(place(LineString(pin.points)).coords)
            self.add_pin(name, points, width=pin.width, gap=pin.gap)

