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

sys.path.insert(0, _I('Design3_3'))
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


def build_design(rounded=True, arm_lengths=ARM_LENGTHS, **options):
    """Build the full chip. Dimensions here are in micrometres."""
    for i, clearance in enumerate((9, 8), 1):
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
    return build_variant(
        gaps=(4, 4), lengths=arm_lengths, rounded=rounded,
        arm_builder=bracket_arm,   # straight ㄈ, upper arm = lower arm on each side
        **options,
    )
