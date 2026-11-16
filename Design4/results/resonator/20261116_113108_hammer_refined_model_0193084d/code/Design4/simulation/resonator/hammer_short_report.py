"""Compare a requested short hammer candidate with the original effective coupling."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from shapely import wkt
import hammer_deembed as H
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('candidate',type=Path);ap.add_argument('baseline_poles',type=Path);ap.add_argument('candidate_poles',type=Path)
    args=ap.parse_args();candidate=args.candidate.resolve();base=args.baseline_poles.resolve();new=args.candidate_poles.resolve()
    params=json.loads((candidate/'run.json').read_text())['parameters']
    oldpole=json.loads((base/'pole_fit.json').read_text());newpole=json.loads((new/'pole_fit.json').read_text())
    full=Path(newpole['source']);baseline=Path(oldpole['source'])
    gap=params['gap_um'];a=oldpole['fits'][-1];b=newpole['fits'][-1]
    fraction=(b['J_MHz']-a['J_MHz'])/a['J_MHz']
    folder=create_run('Design4','resonator','hammer_short_report',entrypoint=__file__,
        inputs=[candidate/'run.json',base/'pole_fit.json',new/'pole_fit.json'],upstream=[str(candidate),str(base),str(new)],
        parameters=dict(EM_rerun=False,target='original effective coupling',screening_relative_tolerance=.05,
            tolerance_source='assistant-selected exploratory comparison criterion; not user or fabrication acceptance'))
    print('RUN_FOLDER='+str(folder),flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    geom=json.loads((candidate/'raw/full_geometry.json').read_text())
    spacing=params['spacing_mm'];tap=H.P.TAP_Y
    ax=axes[0];ax.set_facecolor('#dddddd')
    for subtract in (True,False):
        for kind in ('poly','path'):
            for r in geom[kind]:
                if r['subtract']!=subtract:continue
                shape=wkt.loads(r['wkt'])
                if kind=='path':shape=shape.buffer(float(r['width'])/2,cap_style=2)
                for part in getattr(shape,'geoms',[shape]):
                    if hasattr(part,'exterior'):ax.fill(*part.exterior.xy,facecolor='white' if subtract else '#444444',edgecolor='#555555',lw=.5)
    audit=json.loads((candidate/'geometry_audit.json').read_text())
    for s in audit['ground_straps']:
        p=np.array(s['pos']);v=np.array(s['ori'])*s['span']/2
        ax.plot([p[0]-v[0],p[0]+v[0]],[p[1]-v[1],p[1]+v[1]],color='tab:blue',lw=2)
    for x0,x1,text,y in [(0,.05,'50 um',tap-.06),(spacing-.05,spacing,'50 um',tap-.06),(.075,.075+gap/1000,f'gap {gap:g} um',tap+.05)]:
        ax.annotate('',xy=(x0,y),xytext=(x1,y),arrowprops={'arrowstyle':'<->','color':'tab:red'})
        ax.text((x0+x1)/2,y+.004,text,ha='center',fontsize=9,color='tab:red')
    ax.annotate('',xy=(spacing/2,tap-.015),xytext=(spacing/2,tap+.015),arrowprops={'arrowstyle':'<->','color':'tab:red'})
    ax.text(spacing/2,tap-.035,'head height 30 um',ha='center',fontsize=9,color='tab:red')
    ax.set_xlim(-.04,spacing+.04);ax.set_ylim(tap-.085,tap+.085);ax.set_aspect('equal')
    ax.set_title('Simulation geometry; blue = ground bridges');ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)')
    for source,label in [(baseline,'Original'),(full,f'Arm50 / head30 / gap{gap:g}')]:
        df=pd.read_csv(source/'raw/S_fine.csv',index_col=0)
        axes[1].plot(df.index,20*np.log10(abs(df.S21.map(complex))),label=label)
    axes[1].legend();axes[1].grid(alpha=.25);axes[1].set_xlabel('Frequency (GHz)');axes[1].set_ylabel('HFSS S21 (dB)')
    axes[1].set_title('Main physical line lengths fixed')
    fig.tight_layout();fig.savefig(folder/'candidate.png',dpi=180)
    conv=pd.read_csv(full/'raw/hfss_convergence.csv');after=pd.read_csv(full/'raw/hfss_convergence_after_fine.csv')
    baseline_conv=pd.read_csv(baseline/'raw/hfss_convergence.csv')
    caps=pd.read_csv(candidate/'cpw_sweep.csv').iloc[0].to_dict()
    result=dict(candidate=str(candidate),HFSS_run=str(full),baseline_HFSS_run=str(baseline),
        baseline_pole_run=str(base),candidate_pole_run=str(new),geometry=params,
        baseline_J_MHz=a['J_MHz'],candidate_J_MHz=b['J_MHz'],relative_J_change=fraction,
        within_5pct_exploratory_screen=abs(fraction)<=.05,
        interpretation='Effective J in H/h=[[fp+i*kappa_p/2,J],[J,fr]]; zero intrinsic readout loss assumed.',
        candidate_pole_degree_sensitivity=newpole['J_relative_degree_sensitivity'],
        candidate_S21_pole_fit_rms=b['S21_rms'],capacitance=caps,
        final_delta_s=float(conv.iloc[-1]['Max Mag. Delta S']),delta_s_target=.01,
        baseline_delta_s=float(baseline_conv.iloc[-1]['Max Mag. Delta S']),
        delta_s_target_met=bool(conv.iloc[-1]['Max Mag. Delta S']<=.01),mesh_reused_for_fine=conv.equals(after),
        limitations=['No field-mode assignment or direct-frequency validation; baseline and candidate mesh error remain.',
            'Conditional effective coupling is not an independent direct HFSS scalar.',
            'Negative excess shunts, if present, are signed CPW-reference differences, not physical isolated capacitors.',
            'Requested short arms change main-line spacing below the older 400 um preference; meander pitch is retained.',
            'Two local ground bridges use 12 um simulated diameter to fit the landing region; fabrication signoff pending.'])
    H.write(folder/'comparison.json',result)
    text=f'''# Short hammer candidate: arm50 / head30 / gap{gap:g} um

Actual full-channel HFSS and local Q3D were rerun. Both main line lengths stay 4.535956670430899 / 4.852581441726434 mm. The requested 50 um arms are measured from each head outer edge to its main-line tap. Head width25 um, CPW16/9 um, ground clearance16 um remain unchanged. Centre spacing is {spacing*1000:g} um. Two hammer ground bridges use12 um wire diameter because the17 um landing interval fails the previous15 um wire plus edge-margin check. Other bridges retain their existing definitions. This is simulation geometry, not an adopted fabrication version.

| Quantity | Original | Candidate |
|---|---:|---:|
| Effective J (MHz) | {a['J_MHz']:.6f} | {b['J_MHz']:.6f} |
| Bare Purcell frequency (GHz) | {a['bare_fp_GHz']:.6f} | {b['bare_fp_GHz']:.6f} |
| Bare readout frequency (GHz) | {a['bare_fr_GHz']:.6f} | {b['bare_fr_GHz']:.6f} |

Relative J change: {fraction:+.2%}. These are shared-complex-pole fits of S11/S21 followed by a conditional two-mode inversion, with negligible intrinsic readout loss and real reciprocal coupling. They are not half the visible dip splitting. No physical lengths were changed by fitting. A synthetic known-J recovery check passed. Polynomial background degree sensitivity is {newpole['J_relative_degree_sensitivity']:.3g}; this does not include EM uncertainty.

Candidate final HFSS Delta S: {result['final_delta_s']:.6f}, target0.01. Baseline Delta S: {result['baseline_delta_s']:.6f}. Fine/adaptive history unchanged: {result['mesh_reused_for_fine']}. No direct-frequency or field-mode validation was performed. A5% comparison is only an assistant-selected exploratory screen, not a user-specified acceptance or fabrication criterion.

Local mutual capacitance: {caps['mutual_fF']:.6f} fF. Excess shunts: {caps['excess_a_fF']:.6f} / {caps['excess_b_fF']:.6f} fF. Signed negative excess values can occur after subtracting a more capacitive straight-CPW reference; they are not standalone physical negative capacitors. Coupling preservation is assessed using full-channel spectra, not by enforcing the old mutual capacitance.

[Geometry and actual spectra](candidate.png) · [Exact results and source runs](comparison.json) · [Run provenance](run.json)
'''
    (folder/'README.md').write_text(text,encoding='utf-8')
    finish_run(folder,status='completed',summary=result)


if __name__=='__main__':main()
