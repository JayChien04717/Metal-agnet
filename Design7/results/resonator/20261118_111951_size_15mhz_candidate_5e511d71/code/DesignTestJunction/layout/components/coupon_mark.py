"""Etched F and T identification marks for the test coupons."""

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class CouponMark(QComponent):
    component_metadata = Dict(short_name='couponmark', _qgeometry_table_poly='True')

    def make(self):
        p = self.p
        x, y, s = p.width / 2, p.height / 2, p.stroke
        top = draw.box(-x, y - s, x, y)
        if self.options.letter == 'F':
            letter = top.union(draw.box(-x, -y, -x + s, y))
            letter = letter.union(draw.box(-x, -s / 2, x - s, s / 2))
        elif self.options.letter == 'T':
            letter = top.union(draw.box(-s / 2, -y, s / 2, y))
        else:
            raise ValueError('Coupon mark must be F or T')
        shape = draw.translate(letter, p.pos_x, p.pos_y)
        self.add_qgeometry('poly', {'letter': shape}, subtract=True, layer=p.layer)
