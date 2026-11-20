"""Local flux-bias line geometry shared by transmon and coupler QComponents.

Frame: the my_tmon.py drawing frame (stem / south arm along +y, SQUID ground edge at y = d; the caller flips y,
mirrors x if wanted, rotates and translates). The line runs at y = d + fbl_sep ... d + fbl_sep + w.
The ground strip between the SQUID pocket etch (top edge at y = etch_edge, default d) and the flux-line
gap is d + fbl_sep - cpw_gap - etch_edge. geometry() raises metal_width.MetalWidthError when it is below
the minimum metal width (4 µm, skill/chip_drc.md C12a), so no option set can draw a thinner strip.
B2_FLUX_OPTIONS: fbl_sep 10.5 µm (2026-11-18 ground-strip candidate; the historical 9 µm left a 3 µm strip).

Sketch (pre-flip frame, SQUID above; ## = etched gap, == = line metal, // = ground)::

      ________ SQUID pocket etch (top edge at etch_edge) ________
      /////////////// ground strip, at least 4 um ///////////////
      ##########################################
      ====================================  ##     'L': bar under the SQUID, shorted into
      ##############################  ==  ##          the ground at its left end
                                   ## ==  ##
                                   ## ==  ##     pin at y_end (fbl_height below)

      'T': feed at x = 0, two shorted arms          'U': hairpin, return leg shorted at fbl_drop
          ##########################                ###############################
          ==========================                ==============================
          ##########  ==  ##########                ##  ==  ##################  ==  ##
                   ## == ##                         ##  ==  ## (short)        ##  ==  ##
                   ## == ##  pin                                              ##  ==  ## pin

fbl_type
  'L' (my_tmon.MyTmonCoupler, unchanged): bar x in [-fl/2 - w, fl/2 + w] under the SQUID, shorted into the ground at
      its -x end, leg to the pin at x in [fl/2, fl/2 + w].  Current returns through the ground on one side only.
  'T': the feed comes up at x = 0 (shift with fbl_x so the SQUID sits over ONE arm) and splits into two arms
      x in [-fl/2, fl/2], both ends shorted into the ground.  The return current splits symmetrically: the two
      current loops have opposite moments (quadrupole-like far field -> less crosstalk); only ~half the current
      passes under the SQUID.
Returns (line, gap, port, fake) shapely objects; fake = shorted end segment (LOM only, None for 'T').
"""
from qiskit_metal import draw

from metal_width import min_metal_width, require_strip


B2_FLUX_OPTIONS = dict(fbl_type='L', fbl_sep='10.5um', fbl_height='80um',
                       fbl_length='45um', fbl_x='12um',
                       cpw_width='10um', cpw_gap='6um')


def geometry(d, o, etch_edge=None, unit=1e-3, owner='flux line'):
    """etch_edge: top of the SQUID pocket etch above the line (default d); unit: design units per µm."""
    s, fl, fh, w, g = o['fbl_sep'], o['fbl_length'], o['fbl_height'], o['cpw_width'], o['cpw_gap']
    etch_edge = d if etch_edge is None else etch_edge
    require_strip(owner, 'ground strip between SQUID pocket and flux-line gap (d + fbl_sep - cpw_gap - etch_edge)',
                  d + s - g - etch_edge, min_metal_width(unit))
    y0, y1 = d + s, d + s + w                 # bar bottom / top (pre-flip: away from the SQUID is +y)
    y_end = d + 2 * (fh + s)                  # pin end of the feed
    opt = lambda k, default: o[k] if k in o else default      # qiskit parsed options have no usable .get
    kind = str(opt('fbl_type', 'L')).upper()
    if kind == 'L':
        x0, x1 = -fl / 2 - w, fl / 2 + w
        line = draw.union(draw.box(x0, y0, x1, y1), draw.box(fl / 2, y0, x1, y_end))
        gap = draw.union(draw.box(x0, y0 - g, x1 + g, y0),
                         draw.box(x0, y1, fl / 2, y1 + g),
                         draw.box(fl / 2 - g, y1, fl / 2, y_end),
                         draw.box(x1, y0 - g, x1 + g, y_end))
        port = draw.LineString([(x1, y_end), (fl / 2, y_end)])
        fake = draw.LineString([(x0, y1), (x0, y0)])
    elif kind == 'T':
        xa, xb = -fl / 2, fl / 2
        line = draw.union(draw.box(xa, y0, xb, y1), draw.box(-w / 2, y0, w / 2, y_end))
        gap = draw.union(draw.box(xa, y0 - g, xb, y0),                     # between bar and SQUID ground strip
                         draw.box(xa, y1, -w / 2, y1 + g),                  # below the left arm
                         draw.box(w / 2, y1, xb, y1 + g),                   # below the right arm
                         draw.box(-w / 2 - g, y1, -w / 2, y_end),           # feed CPW gaps
                         draw.box(w / 2, y1, w / 2 + g, y_end))
        port = draw.LineString([(w / 2, y_end), (-w / 2, y_end)])          # same order as 'L' (normal -> pin)
        fake = None
    elif kind == 'U':
        # hairpin: feed leg (right) up to the bar under the SQUID, return leg (left) back down by fbl_drop and shorted
        # into the ground there, far from the loops (Dai et al., PRX Quantum 2, 040313: ground connection away from the
        # control loops). Feed and return run antiparallel fl apart -> small magnetic moment of the line itself.
        drop = opt('fbl_drop', 120e-3)
        xl, xr = -fl / 2 - w, fl / 2 + w            # outer edges of the two legs
        y_short = y1 + drop
        if y_short >= y_end - w:
            raise ValueError('U flux line: fbl_drop too long for fbl_height')
        line = draw.union(draw.box(xl, y0, xr, y1),                       # bar under the SQUID
                          draw.box(fl / 2, y0, xr, y_end),                # feed leg (pin)
                          draw.box(xl, y0, -fl / 2, y_short))             # return leg, shorted at y_short
        gap = draw.union(draw.box(xl - g, y0 - g, xr + g, y0),            # above the bar
                         draw.box(xl - g, y0 - g, xl, y_short),           # outside the return leg
                         draw.box(xr, y0 - g, xr + g, y_end),             # outside the feed leg
                         draw.box(-fl / 2, y1, fl / 2, y1 + g),           # under the bar, between the legs
                         draw.box(-fl / 2, y1, -fl / 2 + g, y_short),     # inside the return leg
                         draw.box(fl / 2 - g, y1, fl / 2, y_end))         # inside the feed leg
        port = draw.LineString([(xr, y_end), (fl / 2, y_end)])
        fake = None
    else:
        raise ValueError(f'fbl_type {kind!r}: use L, T or U')
    items = draw.translate([line, gap, port] + ([fake] if fake is not None else []), o['fbl_x'], 0)
    if str(opt('mirror', 'False')).lower() in ('true', '1'):
        items = draw.scale(items, -1, 1, origin=(0, 0))
        items[2] = draw.LineString(list(items[2].coords)[::-1])
        if fake is not None:
            items[3] = draw.LineString(list(items[3].coords)[::-1])
    line, gap, port = items[:3]
    return line, gap, port, (items[3] if fake is not None else None)
