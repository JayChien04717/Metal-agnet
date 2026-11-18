"""Transmon SQUID base metal with recesses for the Dolan junction leads."""
from qiskit_metal import Dict, draw

from my_tmon import MyTmonCoupler


class MyTmonCouplerFab(MyTmonCoupler):
    """Recess the Dolan contacts with a common signed x offset on the ground side."""

    default_options = Dict(contact_notch_width='2um', contact_notch_depth='4um',
                           ground_notch_offset='2um')

    def _squid_geometry(self):
        island, ground, hole, jj1, jj2 = super()._squid_geometry()
        half_width = self.p.contact_notch_width / 2
        depth = self.p.contact_notch_depth
        offset = self.p.ground_notch_offset
        for junction in (jj1, jj2):
            (x, island_edge), (_, ground_edge) = junction.coords
            left, right = x - half_width, x + half_width
            island = island.difference(draw.box(left, island_edge - depth, right, island_edge))
            ground = ground.difference(draw.box(left + offset, ground_edge,
                                                right + offset, ground_edge + depth))
        return island, ground, hole, jj1, jj2
