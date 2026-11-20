"""Design3e with a straight ㄈ bracket around the inner end of each FX bottom pad."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import sys

from shapely.geometry import LineString
from shapely.ops import unary_union

sys.path.insert(0, _I('Design3'))
from run_q3d import Dict, build_variant, layout as baseline
from fluxline import B2_FLUX_OPTIONS


# Upper and lower arms on each side share the same reach past the pad edge, in µm.
ARM_LENGTHS = (6.7, 18.9)


def bracket_arm(options, pad_edge, pad_half_height, gap, length):
    """Straight ㄈ around the inner end of the FX bottom pad (local coordinates, pad at +x).

    A straight vertical segment runs `gap` from the pad's inner edge; the upper and lower arms
    are straight, equal in length, and reach `length` past that edge (rounded tips). A straight
    neck joins the vertical segment to the coupler bar at the pad's centre line (y = 0)."""
    half = options.pad_arm_width / 2
    x_v = pad_edge - gap - half                      # centre line of the vertical segment
    y_a = pad_half_height + gap + half               # centre lines of the two arms
    x_tip = pad_edge + length - half                 # rounded tip ends `length` past the edge
    if x_tip <= x_v or x_v <= options.pad_arm_length / 2:
        raise ValueError('Bracket arm needs room between the bar end and the pad.')
    bracket = LineString([(x_tip, y_a), (x_v, y_a), (x_v, -y_a), (x_tip, -y_a)]).buffer(
        half, cap_style=1, join_style=2, mitre_limit=2.0)
    neck = LineString([(options.pad_arm_length / 2 - .018, 0), (x_v, 0)]).buffer(half, cap_style=2)
    return unary_union([bracket, neck])


# E-beam junction energies for the fab render (GHz per JJ; skill/ebeam_junction.md).
# FX small JJ 4 GHz is the Design3 input. TC 11.528 / 11.528 GHz is the symmetric assignment the user
# chose on 2026-11-20 (Design3's CZ model uses 16.401 / 6.655 GHz; the half-flux EJ is then lost).
EBEAM_JUNCTIONS = dict(fluxonium={'FX1': 4.0, 'FX2': 4.0}, squid={'TC': (11.528, 11.528)})


def build_design(rounded=True, arm_lengths=ARM_LENGTHS, tc_readout=False, ebeam_junctions=None,
                 ebeam_layers=None, **options):
    """Build the full chip (main FTF layout). Dimensions here are in micrometres.

    tc_readout=False (default): Design3 form, no transmon readout pad / resonator.
    tc_readout=True: add the Design3_2 transmon readout pad, tee and resonator.
    ebeam_junctions: fab render only (skill/gds_export.md section 2.1). None (default) draws no
    junction layers. True uses EBEAM_JUNCTIONS; a dict {'fluxonium': {...}, 'squid': {...}}
    overrides it. Requires jj_metal=False and fab_contacts=True; the drawings sit on the fab
    components' contact hooks (code/general/ebeam_junction_placement.py), on layers 11-14 unless
    ebeam_layers is given. The placement record is stored in design.metadata['ebeam_junctions']."""
    # Charge-line tip clearance cl_sep (um): 10 leaves the 4 um minimum ground strip (cl_sep - cpw_gap 6 um).
    # Historical 9 / 8 um (3 / 2 um strips) are below the rule and now raise in MyFluxonium.
    for i, clearance in enumerate((10, 10), 1):
        fx_options = options.setdefault(f'fx{i}_options', {})
        fx_options.setdefault(
            'charge_line_options',
            Dict(baseline.FX_OPTIONS['charge_line_options'], cl_sep=f'{clearance}um'),
        )
    tmon_options = options.setdefault('tmon_options', {})
    tmon_options.setdefault('paddle_gap', '55um')
    tmon_options.setdefault('pad_head_width', '15um')
    tmon_options.setdefault('flux_bias_line_options', Dict(B2_FLUX_OPTIONS, make_fbl=True))
    if not rounded:
        options.setdefault('include_test_junctions', False)
    if ebeam_junctions and (options.get('jj_metal', True) or not options.get('fab_contacts', False)):
        raise ValueError('ebeam_junctions needs jj_metal=False and fab_contacts=True (fab render only)')
    design = build_variant(
        gaps=(4, 4), lengths=arm_lengths, rounded=rounded, tc_readout=tc_readout,
        arm_builder=bracket_arm,   # straight ㄈ, upper arm = lower arm on each side
        **options,
    )
    if ebeam_junctions:
        from ebeam_junction_placement import add_ebeam_junctions
        spec = EBEAM_JUNCTIONS if ebeam_junctions is True else ebeam_junctions
        design.metadata['ebeam_junctions'] = add_ebeam_junctions(design, layers=ebeam_layers, **spec)
    return design
