"""Parameterized junction-mask geometry referenced to ``layout3 ftf.gds``.

This is a multilayer fabrication drawing, not a lumped HFSS junction model.

Sketch (not to scale; layers keep the source mask meaning)::

      return    num_array per arm (array_length x array_width, array_gap apart)        contact
        _    ________________________________________________     ____
       | |  |_||_||_||_||_||_||_||_||_||_||_||_||_||_||_||_|====|    |
       | |                                                      |    | main
       | |   ________________________________________________   |    | junction
       |_|  |_||_||_||_||_||_||_||_||_||_||_||_||_||_||_||_|====|____| (x 3-6um)
        ^ moves with num_array        fixed: main junction and right-hand contacts

The main junction (junction_gap, junction_width, finger_width) is mapped in
fluxonium_junction_parameters.png (code/general/plot_junction_parameters.py).
All dimensions are parsed in the design's units (normally mm in Qiskit Metal).
See fluxonium_junction.md for source selection and validation boundaries.
"""
from math import isfinite

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class FluxoniumJunction(QComponent):
    """Fixed-size array cells, variable count, and a shared main-junction gap.

    ``num_array`` counts rectangles PER ARM, not electrical junctions. Changing
    it moves the left return and extends the two array bars; the junction and
    right-hand contacts stay fixed. ``array_length`` is along x, ``array_width``
    along y, and ``array_gap`` is the space between cells on layer_array.
    Layer numbers retain the source drawing's mask meaning, not inferred nets.
    ``junction_width`` is the y width of the right-hand horizontal junction
    bar (source 0.25um); it is resized about the fixed source centre line.
    ``finger_width`` is the x width of the vertical junction finger (source
    0.1um); its right edge stays flush with the top bar end at x=3.484um, so
    it grows leftwards and the bridge to the junction bar is unchanged.
    ``contact_extension`` lengthens both end contacts (layer_array) and their
    flanking recess strips outwards beyond the source ends at y=+/-16um; the
    pins move with them. ``mirror=True`` mirrors the local drawing about x=0
    (before rotation) so the array runs to +x with the junction kept on the
    top arm; it also reverses the finger/bar order along x.
    """

    # Source centre line (um) of the right-hand horizontal junction bar.
    JUNCTION_BAR_YC = 2.571
    # Source right edge (um) of the vertical finger = end of the top bar.
    FINGER_RIGHT_X = 3.484
    # Source finger tip (um); the recess below it keeps junction_gap.
    FINGER_TIP_Y = 2.009

    default_options = Dict(
        num_array='43', array_length='0.998um', array_width='2um',
        array_gap='0.202um', junction_gap='0.075um', junction_width='0.25um',
        finger_width='0.1um', contact_extension='0um', mirror='False',
        layer='1', layer_recess='2', layer_array='3', layer_array_bar='4',
        jj_metal='True',
    )
    component_metadata = Dict(short_name='fluxonium_junction',
                              _qgeometry_table_poly='True')

    @staticmethod
    def reference():
        """Return the immutable, device-local source polygon snapshot (um)."""
        return {'source_file': 'layout3 ftf.gds',
         'source_sha256': '6546a6202ce78dc716b100501d88884f8ed6f5854e17b0e59d962bce10ea2b26',
         'source_cell': 'TOP',
         'unit': 'um',
         'selection_window_um': [[-545, -1017], [-475, -983]],
         'translation_um': [485, 1000],
         'regenerated': '2026-11-18; snapshot file was missing from the working tree and was rebuilt from '
                        'the same source hash',
         'polygons_um': {'1': [[[0.501, 3.0], [0.501, 3.5], [3.484, 3.5], [3.484, 3.0]],
                               [[3.384, 2.009], [3.384, 3.0], [3.484, 3.0], [3.484, 2.009]],
                               [[3.732, 2.446], [3.732, 2.696], [5.73, 2.696], [5.73, 2.446]],
                               [[5.73, -3.0], [5.73, 3.412], [6.481, 3.412], [6.481, -3.0]],
                               [[0.501, -3.6], [0.501, -3.0], [6.481, -3.0], [6.481, -3.6]]],
                         '2': [[[2.713, 1.759], [2.713, 2.95], [3.308, 2.95], [3.308, 1.759]],
                               [[3.308, 1.759], [3.308, 1.934], [3.558, 1.934], [3.558, 1.759]],
                               [[6.535, -4.086], [6.535, 3.912], [7.035, 3.912], [7.035, -4.086]],
                               [[5.181, 2.771], [5.181, 3.812], [5.681, 3.812], [5.681, 2.771]],
                               [[5.181, -2.95], [5.181, 2.371], [5.681, 2.371], [5.681, -2.95]],
                               [[0.501, 3.6], [0.501, 16.0], [1.005, 16.0], [1.005, 3.6]],
                               [[-1.005, 5.1], [-1.005, 16.0], [-0.501, 16.0], [-0.501, 5.1]],
                               [[0.501, -16.0], [0.501, -3.7], [1.005, -3.7], [1.005, -16.0]],
                               [[-1.005, -16.0], [-1.005, -5.1], [-0.501, -5.1], [-0.501, -16.0]],
                               [[3.558, 1.759], [3.558, 2.371], [4.131, 2.371], [4.131, 1.759]],
                               [[3.558, 2.371], [3.558, 2.771], [3.657, 2.771], [3.657, 2.371]],
                               [[3.558, 2.771], [3.558, 3.812], [4.131, 3.812], [4.131, 2.771]],
                               [[-54.797, -5.05], [-54.797, 5.049], [-54.111, 5.049], [-54.111, -5.05]],
                               [[-53.118, -2.95], [-53.118, 2.95], [-52.612, 2.95], [-52.612, -2.95]]],
                         '3': [[[-54.1, 3.0], [-54.1, 5.0], [-52.1, 5.0], [-52.1, 3.0]],
                               [[-54.1, -3.0], [-54.1, 3.0], [-53.108, 3.0], [-53.108, -3.0]],
                               [[-54.1, -5.0], [-54.1, -3.0], [-52.1, -3.0], [-52.1, -5.0]],
                               [[-1.5, -5.0], [-1.5, -3.0], [-0.501, -3.0], [-0.501, -5.0]],
                               [[-2.699, -5.0], [-2.699, -3.0], [-1.701, -3.0], [-1.701, -5.0]],
                               [[-3.9, -5.0], [-3.9, -3.0], [-2.9, -3.0], [-2.9, -5.0]],
                               [[-5.1, -5.0], [-5.1, -3.0], [-4.1, -3.0], [-4.1, -5.0]],
                               [[-6.299, -5.0], [-6.299, -3.0], [-5.3, -3.0], [-5.3, -5.0]],
                               [[-7.5, -5.0], [-7.5, -3.0], [-6.5, -3.0], [-6.5, -5.0]],
                               [[-8.699, -5.0], [-8.699, -3.0], [-7.701, -3.0], [-7.701, -5.0]],
                               [[-9.9, -5.0], [-9.9, -3.0], [-8.9, -3.0], [-8.9, -5.0]],
                               [[-11.1, -5.0], [-11.1, -3.0], [-10.101, -3.0], [-10.101, -5.0]],
                               [[-12.299, -5.0], [-12.299, -3.0], [-11.3, -3.0], [-11.3, -5.0]],
                               [[-13.5, -5.0], [-13.5, -3.0], [-12.501, -3.0], [-12.501, -5.0]],
                               [[-14.699, -5.0], [-14.699, -3.0], [-13.701, -3.0], [-13.701, -5.0]],
                               [[-15.9, -5.0], [-15.9, -3.0], [-14.9, -3.0], [-14.9, -5.0]],
                               [[-17.1, -5.0], [-17.1, -3.0], [-16.1, -3.0], [-16.1, -5.0]],
                               [[-18.299, -5.0], [-18.299, -3.0], [-17.3, -3.0], [-17.3, -5.0]],
                               [[-19.5, -5.0], [-19.5, -3.0], [-18.5, -3.0], [-18.5, -5.0]],
                               [[-20.699, -5.0], [-20.699, -3.0], [-19.701, -3.0], [-19.701, -5.0]],
                               [[-21.9, -5.0], [-21.9, -3.0], [-20.9, -3.0], [-20.9, -5.0]],
                               [[-23.1, -5.0], [-23.1, -3.0], [-22.101, -3.0], [-22.101, -5.0]],
                               [[-24.3, -5.0], [-24.3, -3.0], [-23.3, -3.0], [-23.3, -5.0]],
                               [[-25.499, -5.0], [-25.499, -3.0], [-24.501, -3.0], [-24.501, -5.0]],
                               [[-26.699, -5.0], [-26.699, -3.0], [-25.701, -3.0], [-25.701, -5.0]],
                               [[-27.9, -5.0], [-27.9, -3.0], [-26.9, -3.0], [-26.9, -5.0]],
                               [[-29.1, -5.0], [-29.1, -3.0], [-28.101, -3.0], [-28.101, -5.0]],
                               [[-30.299, -5.0], [-30.299, -3.0], [-29.3, -3.0], [-29.3, -5.0]],
                               [[-31.499, -5.0], [-31.499, -3.0], [-30.501, -3.0], [-30.501, -5.0]],
                               [[-32.701, -5.0], [-32.701, -3.0], [-31.701, -3.0], [-31.701, -5.0]],
                               [[-33.9, -5.0], [-33.9, -3.0], [-32.9, -3.0], [-32.9, -5.0]],
                               [[-35.1, -5.0], [-35.1, -3.0], [-34.101, -3.0], [-34.101, -5.0]],
                               [[-36.299, -5.0], [-36.299, -3.0], [-35.3, -3.0], [-35.3, -5.0]],
                               [[-37.499, -5.0], [-37.499, -3.0], [-36.501, -3.0], [-36.501, -5.0]],
                               [[-38.699, -5.0], [-38.699, -3.0], [-37.701, -3.0], [-37.701, -5.0]],
                               [[-39.9, -5.0], [-39.9, -3.0], [-38.9, -3.0], [-38.9, -5.0]],
                               [[-41.1, -5.0], [-41.1, -3.0], [-40.101, -3.0], [-40.101, -5.0]],
                               [[-42.299, -5.0], [-42.299, -3.0], [-41.3, -3.0], [-41.3, -5.0]],
                               [[-43.499, -5.0], [-43.499, -3.0], [-42.501, -3.0], [-42.501, -5.0]],
                               [[-44.7, -5.0], [-44.7, -3.0], [-43.701, -3.0], [-43.701, -5.0]],
                               [[-45.9, -5.0], [-45.9, -3.0], [-44.901, -3.0], [-44.901, -5.0]],
                               [[-47.1, -5.0], [-47.1, -3.0], [-46.101, -3.0], [-46.101, -5.0]],
                               [[-48.299, -5.0], [-48.299, -3.0], [-47.299, -3.0], [-47.299, -5.0]],
                               [[-49.499, -5.0], [-49.499, -3.0], [-48.501, -3.0], [-48.501, -5.0]],
                               [[-50.7, -5.0], [-50.7, -3.0], [-49.701, -3.0], [-49.701, -5.0]],
                               [[-51.9, -5.0], [-51.9, -3.0], [-50.901, -3.0], [-50.901, -5.0]],
                               [[-1.5, 3.0], [-1.5, 5.0], [-0.501, 5.0], [-0.501, 3.0]],
                               [[-2.699, 3.0], [-2.699, 5.0], [-1.701, 5.0], [-1.701, 3.0]],
                               [[-3.9, 3.0], [-3.9, 5.0], [-2.9, 5.0], [-2.9, 3.0]],
                               [[-5.1, 3.0], [-5.1, 5.0], [-4.1, 5.0], [-4.1, 3.0]],
                               [[-6.299, 3.0], [-6.299, 5.0], [-5.3, 5.0], [-5.3, 3.0]],
                               [[-7.5, 3.0], [-7.5, 5.0], [-6.5, 5.0], [-6.5, 3.0]],
                               [[-8.699, 3.0], [-8.699, 5.0], [-7.701, 5.0], [-7.701, 3.0]],
                               [[-9.9, 3.0], [-9.9, 5.0], [-8.9, 5.0], [-8.9, 3.0]],
                               [[-11.1, 3.0], [-11.1, 5.0], [-10.101, 5.0], [-10.101, 3.0]],
                               [[-12.299, 3.0], [-12.299, 5.0], [-11.3, 5.0], [-11.3, 3.0]],
                               [[-13.5, 3.0], [-13.5, 5.0], [-12.501, 5.0], [-12.501, 3.0]],
                               [[-14.699, 3.0], [-14.699, 5.0], [-13.701, 5.0], [-13.701, 3.0]],
                               [[-15.9, 3.0], [-15.9, 5.0], [-14.9, 5.0], [-14.9, 3.0]],
                               [[-17.1, 3.0], [-17.1, 5.0], [-16.1, 5.0], [-16.1, 3.0]],
                               [[-18.299, 3.0], [-18.299, 5.0], [-17.3, 5.0], [-17.3, 3.0]],
                               [[-19.5, 3.0], [-19.5, 5.0], [-18.5, 5.0], [-18.5, 3.0]],
                               [[-20.699, 3.0], [-20.699, 5.0], [-19.701, 5.0], [-19.701, 3.0]],
                               [[-21.9, 3.0], [-21.9, 5.0], [-20.9, 5.0], [-20.9, 3.0]],
                               [[-23.1, 3.0], [-23.1, 5.0], [-22.101, 5.0], [-22.101, 3.0]],
                               [[-24.3, 3.0], [-24.3, 5.0], [-23.3, 5.0], [-23.3, 3.0]],
                               [[-25.499, 3.0], [-25.499, 5.0], [-24.501, 5.0], [-24.501, 3.0]],
                               [[-26.699, 3.0], [-26.699, 5.0], [-25.701, 5.0], [-25.701, 3.0]],
                               [[-27.9, 3.0], [-27.9, 5.0], [-26.9, 5.0], [-26.9, 3.0]],
                               [[-29.1, 3.0], [-29.1, 5.0], [-28.101, 5.0], [-28.101, 3.0]],
                               [[-30.299, 3.0], [-30.299, 5.0], [-29.3, 5.0], [-29.3, 3.0]],
                               [[-31.499, 3.0], [-31.499, 5.0], [-30.501, 5.0], [-30.501, 3.0]],
                               [[-32.701, 3.0], [-32.701, 5.0], [-31.701, 5.0], [-31.701, 3.0]],
                               [[-33.9, 3.0], [-33.9, 5.0], [-32.9, 5.0], [-32.9, 3.0]],
                               [[-35.1, 3.0], [-35.1, 5.0], [-34.101, 5.0], [-34.101, 3.0]],
                               [[-36.299, 3.0], [-36.299, 5.0], [-35.3, 5.0], [-35.3, 3.0]],
                               [[-37.499, 3.0], [-37.499, 5.0], [-36.501, 5.0], [-36.501, 3.0]],
                               [[-38.699, 3.0], [-38.699, 5.0], [-37.701, 5.0], [-37.701, 3.0]],
                               [[-39.9, 3.0], [-39.9, 5.0], [-38.9, 5.0], [-38.9, 3.0]],
                               [[-41.1, 3.0], [-41.1, 5.0], [-40.101, 5.0], [-40.101, 3.0]],
                               [[-42.299, 3.0], [-42.299, 5.0], [-41.3, 5.0], [-41.3, 3.0]],
                               [[-43.499, 3.0], [-43.499, 5.0], [-42.501, 5.0], [-42.501, 3.0]],
                               [[-44.7, 3.0], [-44.7, 5.0], [-43.701, 5.0], [-43.701, 3.0]],
                               [[-45.9, 3.0], [-45.9, 5.0], [-44.901, 5.0], [-44.901, 3.0]],
                               [[-47.1, 3.0], [-47.1, 5.0], [-46.101, 5.0], [-46.101, 3.0]],
                               [[-48.299, 3.0], [-48.299, 5.0], [-47.299, 5.0], [-47.299, 3.0]],
                               [[-49.499, 3.0], [-49.499, 5.0], [-48.501, 5.0], [-48.501, 3.0]],
                               [[-50.7, 3.0], [-50.7, 5.0], [-49.701, 5.0], [-49.701, 3.0]],
                               [[-51.9, 3.0], [-51.9, 5.0], [-50.901, 5.0], [-50.901, 3.0]],
                               [[-0.501, 3.0], [-0.501, 16.0], [0.501, 16.0], [0.501, 3.0]],
                               [[-0.501, -16.0], [-0.501, -3.0], [0.501, -3.0], [0.501, -16.0]]],
                         '4': [[[-54.1, 2.95], [-54.1, 5.05], [-1.069, 5.05], [-1.069, 2.95]],
                               [[-54.1, -5.05], [-54.1, -2.95], [-1.069, -2.95], [-1.069, -5.05]]]}}

    def make(self):
        p = self.p
        n = float(p.num_array)
        if not isfinite(n) or not n.is_integer() or n < 1:
            raise ValueError('num_array must be a positive integer per arm')
        n = int(n)
        u = self.parse_value('1um')
        length, width, spacing, gap, jj_width, f_width = (
            float(v) / u for v in
            (p.array_length, p.array_width, p.array_gap, p.junction_gap,
             p.junction_width, p.finger_width))
        if not all(isfinite(v) and v > 0 for v in
                   (length, width, spacing, gap, jj_width, f_width)):
            raise ValueError('Array dimensions, array_gap, junction_gap, '
                             'junction_width and finger_width must be positive')
        if width > 2:
            raise ValueError('array_width must be <= 2um for the reference 8um arm spacing')

        def box(x0, y0, x1, y1, what):
            # draw.box silently reorders inverted corners, so check explicitly.
            if not (x1 - x0 > 1e-9 and y1 - y0 > 1e-9):
                raise ValueError(f'{what} has non-positive size; reduce '
                                 'junction_gap, junction_width or finger_width')
            return draw.box(x0, y0, x1, y1)
        def flag(value, name):
            text = str(value).lower()
            if text not in ('true', 'false', '1', '0'):
                raise ValueError(f'{name} must be True or False')
            return text in ('true', '1')
        ext = float(p.contact_extension) / u
        if not isfinite(ext) or ext < 0:
            raise ValueError('contact_extension must be >= 0')
        sx = -1 if flag(p.mirror, 'mirror') else 1
        if not flag(p.jj_metal, 'jj_metal'):
            return

        source = self.reference()['polygons_um']
        groups = {int(k): {f'reference_{i}': draw.Polygon(x)
                           for i, x in enumerate(v)} for k, v in source.items()}
        # One common gap controls every red/purple clearance of the junction:
        # top landing, both sides and the tip of the finger, the three sides
        # of the bridge, and both strips beside the bar at the right contact.
        # The right-hand bar's y edges follow junction_width (source 0.25um).
        y0 = self.JUNCTION_BAR_YC - jj_width/2
        y1 = self.JUNCTION_BAR_YC + jj_width/2
        groups[1]['reference_2'] = box(3.732, y0, 5.73, y1, 'junction bar')
        # The vertical finger's left edge follows finger_width (source 0.1um).
        fx1 = self.FINGER_RIGHT_X
        fx0 = fx1 - f_width
        tip = self.FINGER_TIP_Y
        groups[1]['reference_1'] = box(fx0, tip, fx1, 3.0, 'junction finger')
        recess = groups[2]
        recess['reference_0'] = box(2.713, 1.759, fx0-gap, 3-gap, 'recess 0')
        recess['reference_1'] = box(fx0-gap, 1.759, fx1+gap, tip-gap, 'recess 1')
        recess['reference_9'] = box(3.484+gap, 1.759, 4.131, y0-gap, 'recess 9')
        recess['reference_10'] = box(3.484+gap, y0-gap, 3.732-gap, y1+gap,
                                     'recess 10')
        recess['reference_11'] = box(3.484+gap, y1+gap, 4.131, 3.812, 'recess 11')
        # Recess strips beside the bar near the right contact.
        recess['reference_3'] = box(5.181, y1+gap, 5.681, 3.812, 'recess 3')
        recess['reference_4'] = box(5.181, -2.95, 5.681, y0-gap, 'recess 4')
        # Contact-side recess strips 5-8 end at y=+/-16 in the source.
        for i in (5, 6, 7, 8):
            x0_, ya, x1_, yb = recess[f'reference_{i}'].bounds
            recess[f'reference_{i}'] = box(x0_, ya - (ext if ya < 0 else 0), x1_,
                                           yb + (ext if yb > 0 else 0), f'recess {i}')

        pitch = length + spacing
        # Preserve the 0.202um terminal clearance at the left return. This
        # also lets array_length change without moving the rightmost contact.
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
            x0_, ya, x1_, yb = draw.Polygon(coords).bounds
            array[f'contact_{i}'] = box(x0_, ya - (ext if ya < 0 else 0), x1_,
                                        yb + (ext if yb > 0 else 0), f'contact {i}')
        for arm, y in (('top', 4), ('bottom', -4)):
            for i in range(n):
                right = -0.501 - i*pitch
                array[f'{arm}_array_{i:04d}'] = draw.box(
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
                shape = draw.scale(shape, sx*u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layers[source_layer],
                               chip=p.chip, subtract=False)

        for name, y, direction in (('top', 16+ext, 1), ('bottom', -16-ext, -1)):
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
