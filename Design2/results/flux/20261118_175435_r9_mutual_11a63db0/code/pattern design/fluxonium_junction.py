"""Parameterized junction-mask geometry referenced to ``layout3 ftf.gds``.

This is a multilayer fabrication drawing, not a lumped HFSS junction model.
All dimensions are parsed in the design's units (normally mm in Qiskit Metal).
See fluxonium_junction.md for source selection and validation boundaries.
"""
import json
from pathlib import Path
from math import isfinite

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class FluxoniumJunction(QComponent):
    """Fixed-size array cells, variable count, and a shared main-junction gap.

    ``num_cells`` counts rectangles PER ARM, not electrical junctions. Changing
    it moves the left return and extends the two array bars; the junction and
    right-hand contacts stay fixed. ``cell_length`` is along x, ``cell_width``
    along y, and ``cell_gap`` is the space between cells on layer_array.
    Layer numbers retain the source drawing's mask meaning, not inferred nets.
    """

    default_options = Dict(
        num_cells='43', cell_length='0.998um', cell_width='2um',
        cell_gap='0.202um', junction_gap='0.075um',
        layer='1', layer_recess='2', layer_array='3', layer_array_bar='4',
        jj_metal='True',
    )
    component_metadata = Dict(short_name='fluxonium_junction',
                              _qgeometry_table_poly='True')

    @staticmethod
    def reference():
        """Return the immutable, device-local source polygon snapshot (um)."""
        path = Path(__file__).with_name('fluxonium_junction_reference.json')
        return json.loads(path.read_text(encoding='utf-8'))

    def make(self):
        p = self.p
        n = float(p.num_cells)
        if not isfinite(n) or not n.is_integer() or n < 1:
            raise ValueError('num_cells must be a positive integer per arm')
        n = int(n)
        u = self.parse_value('1um')
        length, width, spacing, gap = (
            float(v) / u for v in
            (p.cell_length, p.cell_width, p.cell_gap, p.junction_gap))
        if not all(isfinite(v) and v > 0 for v in (length, width, spacing, gap)):
            raise ValueError('Cell dimensions, cell_gap and junction_gap must be positive')
        if width > 2:
            raise ValueError('cell_width must be <= 2um for the reference 8um arm spacing')
        if gap >= 0.25:
            raise ValueError('junction_gap must be < 0.25um for the reference junction')
        enabled = str(p.jj_metal).lower()
        if enabled not in ('true', 'false', '1', '0'):
            raise ValueError('jj_metal must be True or False')
        if enabled in ('false', '0'):
            return

        source = self.reference()['polygons_um']
        groups = {int(k): {f'reference_{i}': draw.Polygon(x)
                           for i, x in enumerate(v)} for k, v in source.items()}
        # One common gap controls both sides of the finger, the top landing,
        # and all three clearances around the right-hand horizontal bridge.
        recess = groups[2]
        recess['reference_0'] = draw.box(2.713, 1.759, 3.384-gap, 3-gap)
        recess['reference_1'] = draw.box(3.384-gap, 1.759, 3.484+gap, 1.934)
        recess['reference_9'] = draw.box(3.484+gap, 1.759, 4.131, 2.446-gap)
        recess['reference_10'] = draw.box(3.484+gap, 2.446-gap,
                                        3.732-gap, 2.696+gap)
        recess['reference_11'] = draw.box(3.484+gap, 2.696+gap, 4.131, 3.812)

        pitch = length + spacing
        # Preserve the 0.202um terminal clearance at the left return. This
        # also lets cell_length change without moving the rightmost contact.
        last_left = -0.501 - length - (n-1)*pitch
        return_right = last_left - 0.202
        delta = return_right - (-52.1)
        for i in (12, 13):
            recess[f'reference_{i}'] = draw.translate(recess[f'reference_{i}'], delta, 0)

        array = {}
        for i in range(3):
            array[f'return_{i}'] = draw.translate(groups[3][f'reference_{i}'], delta, 0)
        # Source end contacts, independent of cell count.
        for i, coords in enumerate(source['3'][-2:]):
            array[f'contact_{i}'] = draw.Polygon(coords)
        for arm, y in (('top', 4), ('bottom', -4)):
            for i in range(n):
                right = -0.501 - i*pitch
                array[f'{arm}_cell_{i:04d}'] = draw.box(
                    right-length, y-width/2, right, y+width/2)
        groups[3] = array
        for i, y in enumerate((4, -4)):
            groups[4][f'reference_{i}'] = draw.box(
                return_right-2, y-width/2-0.05, -1.069, y+width/2+0.05)
        if return_right-2 >= -1.069:
            raise ValueError('Array is too short for the reference bar termination')

        layers = {1: p.layer, 2: p.layer_recess, 3: p.layer_array, 4: p.layer_array_bar}
        for source_layer, shapes in groups.items():
            placed = {}
            for name, shape in shapes.items():
                if not shape.is_valid or shape.is_empty or shape.area <= 0:
                    raise ValueError(f'Invalid mask geometry: {source_layer}/{name}')
                shape = draw.scale(shape, u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layers[source_layer],
                               chip=p.chip, subtract=False)

        for name, y, direction in (('top', 16, 1), ('bottom', -16, -1)):
            # input_as_norm uses an inner point followed by the terminal point.
            line = draw.LineString([(0, y-direction), (0, y)])
            line = draw.scale(line, u, u, origin=(0, 0))
            line = draw.rotate(line, p.orientation, origin=(0, 0))
            line = draw.translate(line, p.pos_x, p.pos_y)
            self.add_pin(name, np.asarray(line.coords), width=1.002*u,
                         input_as_norm=True, chip=p.chip)


# Accept the spelling originally requested, while keeping a conventional class name.
fluxonium_junction = FluxoniumJunction
fluxonium_juncion = FluxoniumJunction
