"""250 um spacing, 9 um hammer gap: HFSS eigenmode with and without the hammer ground straps.

Question: is the ~130 MHz splitting of hammer250_gap9 a ground-current coupling (ground between the filter and the
resonator cut by stub_f + hammer pocket + stub_r), so that with straps across the two stubs J drops to the
capacitor-model value (Q3D C_j 2.40 fF -> J ~ 16 MHz)?

Cases (each its own HFSS design in Design4/simulation/aedt/Design4.aedt, same setup and mesh seeds):
  control       no straps, res trim -9 um (repeat of hammer250_gap9/matched_eigen_m009 with this script)
  straps_<t>    straps (purcell_6500.draw_hammer_straps), res trim t um for each --trims value
Eigenmode: 2 modes from --minf, first-order basis, Delta F --maxdf %, 3 to --maxp (8) passes, length mesh 20 um lines /
9 um hammer + hugger, the two feed ports as 50 Ohm sheets (so the modes carry kappa), qiskit-metal auto wire bonds
on every CPW, filter / resonator shorted ends shorted.
Afterwards: fit the complex mode frequencies of the strap scan with [[f_p + i k_p/2, J], [J, f_r(t) + i k_r/2]],
f_r(t) = f_r0 + a t, and compare J with the capacitor model.

Output: results/<date>/hammer250_gap9_straps/<case>/ (eigen_modes.json, eigen_convergence.csv, layout.json)
        results/<date>/hammer250_gap9_straps/summary.json

Usage (qmetal env, AEDT open or not):
  python eigenmode/hfss_hammer_straps.py                         all cases (~20 min each)
  python eigenmode/hfss_hammer_straps.py --cases straps --trims=-50,-9 --straps cpw   (use '=' for negative trims)
  AEDT faults on the 2nd 62-strap design in one session: run '--straps cpw' one trim per fresh AEDT
  python eigenmode/hfss_hammer_straps.py --dry-run               layout and strap checks only, no AEDT
  python eigenmode/hfss_hammer_straps.py --fit-only              refit the saved cases
  python eigenmode/hfss_hammer_straps.py --layout-fig            only Design4/fig/hammer250_gap9_straps_before_after.png
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # Design4/simulation
DESIGN = os.path.dirname(SIM); HFSS = os.path.dirname(DESIGN)
sys.path[:0] = [os.path.join(DESIGN, 'pattern design'), os.path.join(HFSS, 'Design5', 'simulation')]
import purcell_6500 as P

LENGTHS = (4.191257715327539, 4.5059730122065185)       # drawn filter / resonator (mm), as hammer250_gap9
SPACING_MM, GAP_UM, NO_HAMMER = 0.25, 9, False      # set from --spacing / --gap / --no-hammer
MESH_LINE, MESH_CAP = '20um', '9um'
CAP_COMPS = ('cj', 'cp')
PORTS = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]
REF = os.path.join(SIM, 'results', '2026-11-11', 'hammer250_gap9')


def prefix():
    return f'P{SPACING_MM * 1e3:.0f}' + ('nh' if NO_HAMMER else f'g{GAP_UM:g}')


def build(trim_um, autobonds=True):
    if NO_HAMMER:      # hammer electrodes removed, ground metal kept between the two stub pockets
        d = P.build_no_hammer(LENGTHS, spacing_mm=SPACING_MM, res_trim=trim_um / 1000)
        t = d.qgeometry.tables['path']
        if autobonds:
            t.loc[t['subtract'].astype(bool), 'hfss_wire_bonds'] = True
        comps = [n for n, c in d.components.items() if type(c).__name__ not in ('ResEnd', 'LaunchpadWirebond')]
        return d, comps
    d = P.build(LENGTHS, paper_hammer=True, paper_hammer_gap_um=GAP_UM, spacing_mm=SPACING_MM,
                res_trim=trim_um / 1000)
    t = d.qgeometry.tables['path']
    if autobonds:     # qiskit auto bonds (>400 um segments only); the 'cpw' strap set replaces them
        t.loc[t['subtract'].astype(bool), 'hfss_wire_bonds'] = True
    comps = [n for n, c in d.components.items() if type(c).__name__ not in ('ResEnd', 'LaunchpadWirebond')]
    return d, comps


def mesh_ops(renderer, d, comps):
    id2n = {c.id: n for n, c in d.components.items()}
    lines, caps = [], []
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            comp = id2n[int(r.component)]
            if not bool(r.subtract) and comp in comps:
                (caps if comp in CAP_COMPS else lines).append(f'{r.name}_{comp}')
    renderer.modeler.mesh_length('cpw_line_mesh', lines, MaxLength=MESH_LINE)
    renderer.modeler.mesh_length('cap_mesh', caps, MaxLength=MESH_CAP)
    print(f'mesh: {MESH_LINE} on {len(lines)} line objects, {MESH_CAP} on {len(caps)} capacitor objects', flush=True)


def run_case(name, trim_um, straps, args, out):
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    import save_aedt
    d, comps = build(trim_um, autobonds=not (straps and args.straps == 'cpw'))
    info = dict(case=name, lengths_mm=[LENGTHS[0], LENGTHS[1] + trim_um / 1000], res_trim_um=trim_um,
                cpw_spacing_um=SPACING_MM * 1000, hammer_gap_um=None if NO_HAMMER else GAP_UM, no_hammer=NO_HAMMER, hammer_straps=straps,
                strap_set=args.straps if straps else None,
                straps=P.ground_straps(d, args.straps) if straps else [], mesh_seed=[MESH_LINE, MESH_CAP],
                min_freq_GHz=args.minf, n_modes=args.modes, max_delta_f_pct=args.maxdf)
    os.makedirs(out, exist_ok=True)
    json.dump(info, open(os.path.join(out, 'layout.json'), 'w'), indent=1)

    class Eig(EigenmodeSim):
        def _analyze(self):
            if straps:
                names = P.draw_hammer_straps(self.renderer, d, which=args.straps)
                print(f'ground straps ({args.straps}) drawn: {names}', flush=True)
            mesh_ops(self.renderer, d, comps)
            super()._analyze()

    eig = Eig(d, 'hfss')
    for k, v in dict(min_freq_ghz=args.minf, n_modes=args.modes, max_delta_f=args.maxdf, max_passes=args.maxp,
                     min_passes=3, min_converged=2, pct_refinement=20, basis_order=1, vars=Dict()).items():
        eig.setup[k] = v
    t0 = time.time()
    design_name = save_aedt.design_name(f'{prefix()}_{name}')
    eig.run(name=design_name, components=comps, open_terminations=[], port_list=PORTS)
    setup = eig.renderer.pinfo.setup
    conv = setup.get_convergence()[0]; conv.to_csv(os.path.join(out, 'eigen_convergence.csv'))
    f, im = setup.get_solutions().eigenmodes()
    modes = dict(frequencies_GHz=[float(x) for x in f], imag_GHz=[float(x) for x in im],
                 passes=int(len(conv)), elements=int(conv.iloc[-1, 0]), last_delta_f_pct=float(conv.iloc[-1, -1]),
                 minutes=round((time.time() - t0) / 60, 1), aedt_design=design_name + '_hfss')
    json.dump(modes, open(os.path.join(out, 'eigen_modes.json'), 'w'), indent=1)
    print(f'{name}: ' + ', '.join(f'{a:.5f} GHz (kappa {2e3 * b:.1f} MHz)' for a, b in zip(f, im))
          + f'; split {1e3 * (max(f) - min(f)):.1f} MHz; {modes["passes"]} passes, {modes["minutes"]} min', flush=True)
    return modes


def recover_case(name, out):
    """Read the solved passes of an aborted / failed case from the open AEDT project (newest design of that case)."""
    from pyEPR.ansys import HfssApp
    proj = HfssApp().get_app_desktop().get_active_project()
    hits = sorted(x.name for x in proj.get_designs() if x.name.startswith(f'{prefix()}_{name}_'))
    if not hits:
        raise RuntimeError(f'no AEDT design for {name}')
    des = [x for x in proj.get_designs() if x.name == hits[-1]][0]
    setup = des.get_setup('Setup')
    conv = setup.get_convergence()[0]; conv.to_csv(os.path.join(out, 'eigen_convergence.csv'))
    f, im = setup.get_solutions().eigenmodes()
    modes = dict(frequencies_GHz=[float(x) for x in f], imag_GHz=[float(x) for x in im],
                 passes=int(len(conv)), elements=int(conv.iloc[-1, 0]), last_delta_f_pct=float(conv.iloc[-1, -1]),
                 aedt_design=hits[-1], note='recovered from AEDT after the solve was stopped')
    json.dump(modes, open(os.path.join(out, 'eigen_modes.json'), 'w'), indent=1)
    print(f'{name} (recovered, {modes["passes"]} passes): '
          + ', '.join(f'{a:.5f} GHz (kappa {2e3 * b:.1f} MHz)' for a, b in zip(f, im))
          + f'; split {1e3 * (max(f) - min(f)):.1f} MHz', flush=True)
    return modes


def point_J(freqs, imags):
    """Exact 2-mode inversion of ONE eigen solve (resonator intrinsic loss set to 0):
    trace and determinant of [[f_p + i k_p/2, J], [J, f_r]] -> f_p, k_p, f_r, J. Robust when the
    scan model (fixed f_p, linear f_r) does not hold, e.g. weak coupling (J << k_p / 4)."""
    from scipy.optimize import least_squares
    z = sorted(np.array(freqs) + 1j * np.array(imags), key=lambda c: c.real)
    s_, p_ = z[0] + z[1], z[0] * z[1]

    def r(x):
        fp, kp, fr, J = x; a = fp + 0.5j * kp
        e1, e2 = a + fr - s_, a * fr - J ** 2 - p_
        return 1e3 * np.array([e1.real, e1.imag, e2.real, e2.imag])
    wide = int(np.argmax([c.imag for c in z]))
    sol = least_squares(r, [z[wide].real, 2 * (z[0].imag + z[1].imag), z[1 - wide].real, 0.002])
    fp, kp, fr, J = sol.x
    return dict(f_p_GHz=fp, kappa_p_MHz=1e3 * kp, f_r_GHz=fr, detuning_MHz=1e3 * (fr - fp), J_MHz=1e3 * abs(J))


def fit_scan(points):
    """points: [(trim_um, [f1, f2], [im1, im2])] -> 2x2 coupled-mode fit (GHz in, MHz out)."""
    from scipy.optimize import least_squares
    t = np.array([p[0] for p in points], float)
    z = np.array([sorted(np.array(p[1]) + 1j * np.array(p[2]), key=lambda c: c.real) for p in points])

    def eig(x, tt):
        fp, kp, fr0, a, J, kr = x
        out = []
        for ti in tt:
            ev = np.linalg.eigvals([[fp + 0.5j * kp, J], [J, fr0 + a * ti + 0.5j * kr]])
            out.append(sorted(ev, key=lambda c: c.real))
        return np.array(out)

    def res(x):
        r = (eig(x, t) - z).ravel() * 1e3        # MHz
        return np.r_[r.real, r.imag]

    c = z.real.mean()
    best = None
    for J0 in (0.005, 0.02, 0.06):
        x0 = [c, 0.035, c, -0.0013, J0, 0.001]
        s = least_squares(res, x0, bounds=([c - .3, 0, c - .3, -.005, 0, 0], [c + .3, .2, c + .3, 0, .2, .05]))
        if best is None or s.cost < best.cost:
            best = s
    fp, kp, fr0, a, J, kr = best.x
    t0 = (fp - fr0) / a
    return dict(f_p_GHz=fp, kappa_p_MHz=1e3 * kp, f_r0_GHz=fr0, df_r_dtrim_MHz_per_um=1e3 * a, J_MHz=1e3 * J,
                kappa_r_MHz=1e3 * kr, crossing_trim_um=t0, rms_MHz=float(np.sqrt(np.mean(res(best.x) ** 2))),
                n_points=len(points),
                note='J from the real mode frequencies and widths; with one trim only J is not separable from detuning')


def summarize(root, cases):
    got = {}
    for name in cases:
        p = os.path.join(root, name, 'eigen_modes.json')
        if os.path.exists(p):
            got[name] = dict(json.load(open(p)), **json.load(open(os.path.join(root, name, 'layout.json'))))
    s = dict(cases={k: dict(trim_um=v['res_trim_um'], straps=v['hammer_straps'], f_GHz=v['frequencies_GHz'],
                            kappa_MHz=[2e3 * x for x in v['imag_GHz']],
                            split_MHz=1e3 * (max(v['frequencies_GHz']) - min(v['frequencies_GHz'])))
                 for k, v in got.items()})
    scan = [(v['res_trim_um'], v['frequencies_GHz'], v['imag_GHz']) for v in got.values() if v['hammer_straps']]
    if len(scan) >= 3:
        s['straps_fit'] = fit_scan(scan)
    s['per_point'] = {k: dict(trim_um=v['res_trim_um'], **point_J(v['frequencies_GHz'], v['imag_GHz']))
                      for k, v in got.items()}
    cmp_p = os.path.join(REF, 'capacitance_comparison.json')
    if os.path.exists(cmp_p) and SPACING_MM == 0.25 and GAP_UM == 9 and not NO_HAMMER:
        c = json.load(open(cmp_p))
        s['capacitor_model'] = dict(q3d_cj_fF=c['q3d_cj_fF'], J_MHz=c['model_J_MHz'], kappa_p_MHz=c['model_kappa_p_MHz'],
                                    mhz_per_fF=c['mhz_per_fF'], no_strap_hfss_J_MHz=c['hfss_J_MHz'])
        if 'straps_fit' in s:
            s['straps_fit']['equivalent_Cj_fF'] = s['straps_fit']['J_MHz'] / c['mhz_per_fF']
    json.dump(s, open(os.path.join(root, 'summary.json'), 'w'), indent=1, default=float)
    print('\n=== summary ===')
    for k, v in s['cases'].items():
        print(f'{k:>12}: trim {v["trim_um"]:+.0f} um, straps {v["straps"]}, modes '
              + ', '.join(f'{f:.4f}' for f in v['f_GHz']) + f' GHz, split {v["split_MHz"]:.1f} MHz')
    if 'straps_fit' in s:
        fz = s['straps_fit']
        print(f'straps fit: J {fz["J_MHz"]:.1f} MHz, kappa_p {fz["kappa_p_MHz"]:.1f} MHz, f_p {fz["f_p_GHz"]:.4f} GHz, '
              f'crossing at trim {fz["crossing_trim_um"]:+.0f} um, rms {fz["rms_MHz"]:.2f} MHz'
              + (f', equivalent C_j {fz["equivalent_Cj_fF"]:.2f} fF' if 'equivalent_Cj_fF' in fz else ''))
    for k, v in s['per_point'].items():
        print(f'{k:>12}: per-point inversion J {v["J_MHz"]:.2f} MHz, detuning {v["detuning_MHz"]:+.1f} MHz, '
              f'kappa_p {v["kappa_p_MHz"]:.1f} MHz', flush=True)
    if 'capacitor_model' in s:
        m = s['capacitor_model']
        print(f'capacitor model: Q3D C_j {m["q3d_cj_fF"]:.2f} fF -> J {m["J_MHz"]:.1f} MHz, kappa_p {m["kappa_p_MHz"]:.1f} MHz '
              f'(no-strap HFSS J {m["no_strap_hfss_J_MHz"]:.1f} MHz)')
    print('summary:', os.path.join(root, 'summary.json'))
    return s


def plot_layout(path, trim_um=-9, which='all'):
    """Before / after figure: overview of the channel and a zoom on the hammer, ground split vs strapped."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from shapely.geometry import box
    from shapely.ops import unary_union
    d, _ = build(trim_um)
    etch = P._etched_region(d)
    metal = unary_union([r.geometry.buffer(r.width / 2, cap_style=2) for r in d.qgeometry.tables['path'].itertuples()
                         if not r.subtract] + [r.geometry for r in d.qgeometry.tables['poly'].itertuples() if not r.subtract])
    straps = P.ground_straps(d, which)
    y0 = P.TAP_Y; xf, xr = 0.0, SPACING_MM

    def draw_straps(a, ms=7, lw=3):
        for s in straps:
            p, o = np.array(s['pos']), np.array(s['ori']) * s['span'] / 2
            xs, ys = [p[0] - o[0], p[0] + o[0]], [p[1] - o[1], p[1] + o[1]]
            a.plot(xs, ys, color='crimson', lw=lw, solid_capstyle='round', zorder=5)
            a.plot(xs, ys, 'o', color='crimson', ms=ms, zorder=6)

    def polys(g):
        return [g] if g.geom_type == 'Polygon' else [p for p in getattr(g, 'geoms', []) if p.geom_type == 'Polygon']

    def fill(ax, g, **kw):
        for p in polys(g):
            ax.fill(*p.exterior.xy, lw=0, **kw)
            for r in p.interiors:
                ax.fill(*r.xy, lw=0, color='white')

    fig = plt.figure(figsize=(14, 8.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1, 1], height_ratios=[1, 1])
    # overview
    ax = fig.add_subplot(gs[:, 0])
    win = box(-1.9, -0.25, 1.9, 4.0)
    fill(ax, etch.intersection(win), color='0.6'); fill(ax, metal.intersection(win), color='#c08a2b')
    ax.add_patch(plt.Rectangle((-0.06, y0 - 0.1), SPACING_MM + 0.12, 0.2, fill=False, ec='r', lw=1.5))
    ax.annotate('hammer (zoom)', (SPACING_MM + 0.06, y0), (0.7, y0 + 0.5), color='r', fontsize=9,
                arrowprops=dict(arrowstyle='->', color='r'))
    ax.plot([s['pos'][0] for s in straps], [s['pos'][1] for s in straps], 'x', color='crimson', ms=6, mew=2,
            label=f'added ground straps ({len(straps)})')
    ax.legend(loc='upper right', fontsize=8)
    if which in ('all', 'cpw'):   # inset: the feed tee n0 and the stub to the hugger
        ins = ax.inset_axes([0.03, 0.2, 0.36, 0.26])
        nwin = box(-0.12, -0.08, 0.12, 0.17)
        fill(ins, etch.intersection(nwin), color='0.6'); fill(ins, metal.intersection(nwin), color='#c08a2b')
        draw_straps(ins, ms=4, lw=2)
        ins.set_xlim(-0.12, 0.12); ins.set_ylim(-0.08, 0.17); ins.set_aspect('equal')
        ins.set_xticks([]); ins.set_yticks([]); ins.set_title('feed tee n0 (AFTER)', fontsize=8)
    ax.set_xlim(-1.9, 1.9); ax.set_ylim(-0.25, 4.0); ax.set_aspect('equal')
    ax.set_title(f'channel, {SPACING_MM * 1e3:.0f} µm spacing, gap {GAP_UM} µm\n'
                 'grey = etch, gold = metal (CPW centre / hammer)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    # zooms: before / after
    zx0, zx1, zy0, zy1 = xf - 0.04, xr + 0.04, y0 - 0.075, y0 + 0.075
    zwin = box(zx0, zy0, zx1, zy1)
    inner = box(xf + 0.011, zy0, xr - 0.011, zy1).difference(etch)
    upper = unary_union([p for p in polys(inner) if p.centroid.y > y0])
    lower = unary_union([p for p in polys(inner) if p.centroid.y < y0])
    after = dict(hammer='AFTER: hammer ground straps', all='AFTER: hammer + T-junction straps',
                 cpw='AFTER: straps on every CPW')[which]
    for col, (title, strapped) in enumerate([('BEFORE: no straps', False), (after, True)]):
        for row, scale in enumerate((1.0, None)):
            a = fig.add_subplot(gs[row, col + 1])
            if row == 1:   # second row: schematic of the ground islands between the two lines
                fill(a, etch.intersection(zwin), color='0.6')
                if strapped:
                    fill(a, inner, color='#5aa05a')
                else:
                    fill(a, upper, color='#4c78c8'); fill(a, lower, color='#e07b39')
                fill(a, metal.intersection(zwin), color='#c08a2b')
                txt = ('upper and lower ground between the lines: one net (green)' if strapped else
                       'ground between the lines cut in two (blue / orange):\nreturn current detours → extra J ≈ 60 MHz')
                a.set_title(txt, fontsize=8.5)
            else:
                fill(a, etch.intersection(zwin), color='0.6'); fill(a, metal.intersection(zwin), color='#c08a2b')
                a.set_title(title, fontsize=11, fontweight='bold', color='#2a7a2a' if strapped else '#b03030')
            if strapped:
                draw_straps(a)
                if row == 0:
                    a.annotate(f'bond wire, span {straps[0]["span"] * 1e3:.0f} µm\nØ {straps[0]["diameter"] * 1e3:.0f} µm, PEC',
                               (straps[1]['pos'][0], y0 - 0.03), (SPACING_MM / 2, y0 - 0.068), fontsize=8,
                               color='crimson', ha='center', arrowprops=dict(arrowstyle='->', color='crimson'))
            else:
                a.annotate('etch: stub_f + pocket + stub_r', (SPACING_MM / 2, y0 - 0.03), (SPACING_MM / 2, y0 - 0.065),
                           ha='center', fontsize=8, color='#b03030', arrowprops=dict(arrowstyle='->', color='#b03030'))
            a.text(xf, zy1 - 0.006, 'filter', ha='center', va='top', fontsize=8)
            a.text(xr, zy1 - 0.006, 'resonator', ha='center', va='top', fontsize=8)
            a.set_xlim(zx0, zx1); a.set_ylim(zy0, zy1); a.set_aspect('equal')
            a.tick_params(labelsize=7)
    fig.suptitle('Hammer coupler before / after ground straps (Wallraff: air bridges C_AB at the hammer T-junctions)',
                 fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140); plt.close(fig)
    print('layout figure:', path, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', default='control,straps', help='control, straps or both (comma separated)')
    ap.add_argument('--trims', default='-50,-9,30', help='resonator trims (um) for the strap scan')
    ap.add_argument('--control-trim', type=float, default=-9)
    ap.add_argument('--minf', type=float, default=5.95)
    ap.add_argument('--modes', type=int, default=2)
    ap.add_argument('--maxdf', type=float, default=0.1)
    ap.add_argument('--spacing', type=float, default=0.25, help='filter-resonator centreline spacing (mm)')
    ap.add_argument('--gap', type=float, default=9, help='hammer gap (um)')
    ap.add_argument('--no-hammer', action='store_true', help='hammer electrodes removed (parasitic line coupling)')
    ap.add_argument('--straps', choices=('hammer', 'all', 'cpw'), default='cpw',
                    help="'hammer': 2 hammer-stub straps; 'all': + every T-junction arm; 'cpw': + every CPW segment (replaces auto bonds)")
    ap.add_argument('--maxp', type=int, default=8, help='max adaptive passes (Jay 2026-11-11: 8 is enough)')
    ap.add_argument('--out', default=os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'), 'hammer250_gap9_straps'))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--fit-only', action='store_true')
    ap.add_argument('--layout-fig', action='store_true', help='only write the before / after layout figure')
    args = ap.parse_args()
    global SPACING_MM, GAP_UM, NO_HAMMER
    SPACING_MM, GAP_UM, NO_HAMMER = args.spacing, args.gap, args.no_hammer
    fig_path = os.path.join(DESIGN, 'fig', f'hammer{SPACING_MM * 1e3:.0f}_gap{GAP_UM:g}_straps_before_after.png')
    if args.layout_fig:
        plot_layout(fig_path, which=args.straps); return
    if not NO_HAMMER:
        plot_layout(fig_path, which=args.straps)
    try:
        os.makedirs(args.out, exist_ok=True); import shutil; shutil.copy(fig_path, args.out)
    except OSError:
        pass

    want = args.cases.split(',')
    plan = []
    if 'control' in want:
        plan.append(('control', args.control_trim, False))
    if 'straps' in want:
        prefix = dict(hammer='straps', all='allwb', cpw='cpwwb')[args.straps]
        plan += [(f'{prefix}_{float(t):+.0f}um'.replace('+', 'p').replace('-', 'm'), float(t), True)
                 for t in args.trims.split(',')]

    for name, trim, straps in plan:                       # layout checks before any solve (fails fast)
        d, comps = build(trim)
        st = P.ground_straps(d, args.straps) if straps else []
        print(f'{name}: trim {trim:+.0f} um, {len(st)} straps, {len(comps)} components: '
              + ', '.join(s['where'] for s in st), flush=True)
    if args.dry_run:
        return
    if args.fit_only:                                   # pick up cases that stopped early (solved passes in AEDT)
        for name, _, _ in plan:
            out = os.path.join(args.out, name)
            if os.path.exists(os.path.join(out, 'layout.json')) and not os.path.exists(os.path.join(out, 'eigen_modes.json')):
                try:
                    recover_case(name, out)
                except Exception as e:
                    print(f'{name} not recovered: {e!r}', flush=True)
    else:
        import hfss_dirs, save_aedt
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(SIM))
        for name, trim, straps in plan:
            try:
                run_case(name, trim, straps, args, os.path.join(args.out, name))
            except Exception as e:
                print(f'!!! {name} failed: {e!r}', flush=True)
        try:
            save_aedt.snapshot(args.out, 'hammer250_gap9_straps')
        except Exception as e:
            print('AEDT project not snapshotted:', repr(e)[:150])
    summarize(args.out, [p[0] for p in plan])


if __name__ == '__main__':
    main()
