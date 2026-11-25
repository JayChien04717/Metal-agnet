r"""Transmon SQUID base metal with recesses for the Dolan junction leads.

Sketch of the SQUID block of MyTmonCouplerDolan (pre-flip frame; the stem / island is on top)::

              island electrode
         ______________________________
        |      |__|       |__|         |   <- island notches (notch_width x notch_depth),
        |______|  |_______|  |_________|      shifted by island_notch_offset in x
              JJ1  ---cut---  JJ2
         _______  _________  __________
        |       ||  hole   ||          |   <- ground-side notches, shifted by
        |       |__|       |__|        |      ground_notch_offset in x
        |__________________________________ ground electrode -> chip ground
"""
from qiskit_metal import Dict, draw
from shapely.geometry import Point

from my_tmon_dolan import MyTmonCouplerDolan


class MyTmonCouplerDolanFab(MyTmonCouplerDolan):
    """Recess the Dolan contacts with signed x offsets on the island and ground sides."""

    default_options = Dict(contact_notch_width='2um', contact_notch_depth='4um',
                           island_notch_offset='0um', ground_notch_offset='2um')

    def _squid_geometry(self):
        island, ground, hole, jj1, jj2 = super()._squid_geometry()
        half_width = self.p.contact_notch_width / 2
        depth = self.p.contact_notch_depth
        island_offset = self.p.island_notch_offset
        ground_offset = self.p.ground_notch_offset
        for junction in (jj1, jj2):
            (x, island_edge), (_, ground_edge) = junction.coords
            left, right = x - half_width, x + half_width
            island = island.difference(draw.box(left + island_offset, island_edge - depth,
                                                right + island_offset, island_edge))
            ground = ground.difference(draw.box(left + ground_offset, ground_edge,
                                                right + ground_offset, ground_edge + depth))
        return island, ground, hole, jj1, jj2

    def contact_hooks(self):
        """Placed contact notches, in design units (mm), for drawing e-beam junctions on them.

        Returns a list of dicts with ``name`` ('island_jj1', 'ground_jj1', 'island_jj2',
        'ground_jj2'), ``mouth`` (notch opening at the junction gap), ``end`` (closed notch end)
        and ``width``. Call after the component is built.
        """
        p = self.p
        _, _, _, jj1, jj2 = MyTmonCouplerDolan._squid_geometry(self)   # JJ lines are not notched
        depth = p.contact_notch_depth
        island_offset, ground_offset = p.island_notch_offset, p.ground_notch_offset
        hooks = []
        for tag, line in (('jj1', jj1), ('jj2', jj2)):
            (x, island_edge), (_, ground_edge) = line.coords
            pts = self._place([Point(x + island_offset, island_edge),
                               Point(x + island_offset, island_edge - depth),
                               Point(x + ground_offset, ground_edge),
                               Point(x + ground_offset, ground_edge + depth)])
            xy = [(q.x, q.y) for q in pts]
            hooks += [dict(name=f'island_{tag}', mouth=xy[0], end=xy[1], width=p.contact_notch_width),
                      dict(name=f'ground_{tag}', mouth=xy[2], end=xy[3], width=p.contact_notch_width)]
        return hooks
