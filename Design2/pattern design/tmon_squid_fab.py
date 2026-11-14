"""Design2 transmon with the Design3 SQUID landing geometry."""
from qiskit_metal import Dict, draw

from my_tmon import MyTmonCoupler
from my_tmon_fab import MyTmonCouplerFab
from tmon_c import Tmon_c


class TmonSQUID(Tmon_c):
    """Keep the Design2 body and align the Design3 SQUID to its ground edge."""

    default_options = Dict(
        squid_width='40um',
        squid_loop_shape='rectangle',
        squid_loop_width='20um',
        squid_loop_height='24um',
        squid_contact_overlap='3um',
        squid_shift_toward_stem='2um',
        squid_span='40um',
        jj_cut='5um',
        jj_width='7um',
    )

    def _squid_geometry(self):
        p = self.p
        yc = p.pad_head_length + p.jj_gap - p.squid_span / 2
        half_span = p.squid_span / 2 + p.squid_contact_overlap
        half_width = p.squid_width / 2
        half_cut = p.jj_cut / 2
        hole = draw.rectangle(p.squid_loop_width, p.squid_loop_height, 0, yc)
        island = draw.box(-half_width, yc - half_span, half_width, yc - half_cut).difference(hole)
        ground = draw.box(-half_width, yc + half_cut, half_width, yc + half_span).difference(hole)
        jj_x = (p.squid_loop_width + p.squid_width) / 4
        junctions = [draw.LineString([(x, yc - half_cut), (x, yc + half_cut)]) for x in (-jj_x, jj_x)]
        return draw.translate([island, ground, hole, *junctions], 0, -p.squid_shift_toward_stem)

    def make_pocket(self):
        # Tmon_c supplies the original stem, arms and etch with jj_metal=False.
        super().make_pocket()
        p = self.p
        island, ground, hole, jj1, jj2 = self._squid_geometry()
        # The original Design2 pocket is 10 um shorter along the stem.
        # Use the unchanged pocket wall as the ground landing reference.
        parts = [island, ground, hole, jj1, jj2]
        parts = draw.rotate(parts, p.orientation, origin=(0, 0))
        island, ground, hole, jj1, jj2 = draw.translate(parts, p.pos_x, p.pos_y)

        polygons = self.design.qgeometry.tables['poly']
        stem, = polygons.index[(polygons.component == self.id) & (polygons.name == 'tmon')]
        etch, = polygons.index[(polygons.component == self.id) & (polygons.name == 'pad_etch')]
        polygons.at[stem, 'geometry'] = polygons.at[stem, 'geometry'].union(island).difference(hole)
        # Retain the grounded landing as part of the ground plane.
        polygons.at[etch, 'geometry'] = polygons.at[etch, 'geometry'].difference(ground)

        junctions = self.design.qgeometry.tables['junction']
        self.design.qgeometry.tables['junction'] = junctions.loc[junctions.component.ne(self.id)]
        self.add_qgeometry('junction', dict(rect_jj1=jj1, rect_jj2=jj2), width=p.jj_width)


class TmonSQUIDFab(TmonSQUID):
    """Apply Design3's same-direction Dolan recesses to the SQUID landings."""

    default_options = MyTmonCouplerFab.default_options

    def _squid_geometry(self):
        island, ground, hole, jj1, jj2 = super()._squid_geometry()
        p = self.p
        half_width = p.contact_notch_width / 2
        for junction in (jj1, jj2):
            (x, island_edge), (_, ground_edge) = junction.coords
            left, right = x - half_width, x + half_width
            island = island.difference(draw.box(left, island_edge - p.contact_notch_depth,
                                                right, island_edge))
            ground = ground.difference(draw.box(left + p.ground_notch_offset, ground_edge,
                                                right + p.ground_notch_offset,
                                                ground_edge + p.contact_notch_depth))
        return island, ground, hole, jj1, jj2
