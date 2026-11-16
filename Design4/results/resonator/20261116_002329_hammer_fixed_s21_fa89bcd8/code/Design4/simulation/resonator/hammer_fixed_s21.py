"""Fixed-geometry, fixed-length Design4 comparison. One HFSS mesh for all sweeps.

Local Q3D hammer data are upstream. Identical cin/cp coupons are extracted once.
Raw/excess hammer shunts are compared with exactly the same TL network, without
fitting lengths, phase, amplitudes, impedance, velocity or any capacitance.
"""
from pathlib import Path
import argparse
import json
import time
import traceback
import shutil
import numpy as np
import pandas as pd
import hammer_deembed as H
from hammer_deembed import P, ROOT, write, geometry_dump
from run_record import create_run, finish_run


def native_signature(design):
    editor=design._modeler
    data={'units':editor.GetModelUnits(),'objects':{},'boundaries':{}}
    for group in ('Solids','Sheets'):
        for obj in sorted(editor.GetObjectsInGroup(group) or []):
            vertices=sorted(tuple(round(float(x),9) for x in editor.GetVertexPosition(v))
                for v in editor.GetVertexIDsFromObject(obj))
            areas=sorted(round(float(editor.GetFaceArea(f)),9) for f in editor.GetFaceIDs(obj))
            props={k:editor.GetPropertyValue('Geometry3DAttributeTab',obj,k)
                for k in ('Material','Solve Inside')}
            data['objects'][obj]=dict(group=group,vertices=vertices,face_areas=areas,properties=props)
    boundaries=list(design._boundaries.GetBoundaries())
    for name,kind in zip(boundaries[::2],boundaries[1::2]):
        server='BoundarySetup:'+name
        props={k:design._design.GetPropertyValue('HfssTab',server,k)
            for k in design._design.GetProperties('HfssTab',server)}
        data['boundaries'][name]=dict(kind=kind,properties=props)
    return data


def model(f, cap, hammer, lengths, z0, velocity, feed_lengths, spacing=.4, cap_half_span=0.):
    """Nodal TL network; pin planes at enclosing caps, centre planes at hammer."""
    f = np.asarray(f, float)
    omega = 2*np.pi*f*1e9
    # 0,1 external ports; 2,3 cin a/b; 4 tee n0; 5,6 cp a/b;
    # 7,8 filter/resonator taps; 9,10 hammer centre planes.
    Y = np.zeros((len(f),11,11), complex)
    def tl(a,b,length):
        theta = omega/velocity*length*1e-3
        diag = -1j/z0/np.tan(theta); off = 1j/z0/np.sin(theta)
        Y[:,a,a] += diag; Y[:,b,b] += diag
        Y[:,a,b] += off; Y[:,b,a] += off
    def pi(a,b,values):
        cc, ga, gb = np.asarray(values)*1e-15
        Y[:,a,a] += 1j*omega*(cc+ga); Y[:,b,b] += 1j*omega*(cc+gb)
        Y[:,a,b] -= 1j*omega*cc; Y[:,b,a] -= 1j*omega*cc
    tl(0,2,feed_lengths[0]); pi(2,3,cap); tl(3,4,feed_lengths[1]); tl(4,1,feed_lengths[2])
    tl(4,5,feed_lengths[3]); pi(5,6,cap); tl(6,7,P.FILTER_TAP)
    Y[:,7,7] += -1j/z0/np.tan(omega/velocity*(lengths[0]-P.FILTER_TAP)*1e-3)
    Y[:,8,8] += 1j/z0*np.tan(omega/velocity*P.RES_TAP*1e-3)
    Y[:,8,8] += -1j/z0/np.tan(omega/velocity*(lengths[1]-P.RES_TAP)*1e-3)
    tl(7,9,spacing/2-cap_half_span); pi(9,10,hammer); tl(10,8,spacing/2-cap_half_span)
    port = Y[:,:2,:2] - Y[:,:2,2:] @ np.linalg.solve(Y[:,2:,2:], Y[:,2:,:2])
    eye = np.broadcast_to(np.eye(2), port.shape)
    # Lumped HFSS ports are 50 ohm, even when the physical line Z0 differs.
    S = np.linalg.solve(eye+50*port, eye-50*port)
    return S[:,0,0], S[:,1,0]


def compare(folder, upstream, cap, raw_hammer, d):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.special import ellipk
    from scipy.signal import find_peaks
    sweep = pd.read_csv(upstream/'cpw_sweep.csv')
    H.plot_sweep(sweep.to_dict('records'),folder)
    audit = {}
    for case in sorted((upstream/'raw').iterdir()):
        files = sorted(case.glob('pass_*.csv'), key=lambda p:int(p.stem.split('_')[-1]))
        if len(files)<2: continue
        a,b = [pd.read_csv(p,index_col=0) for p in files[-2:]]
        sig = [n for n in b.index if 'ground' not in n.lower()]
        aa,bb = [x.loc[sig,sig].to_numpy(float) for x in (a,b)]
        audit[case.name] = dict(last_pass_relative_change=float(np.max(abs(bb-aa)/np.maximum(abs(bb),1e-30))),
            source=[str(p.relative_to(upstream)) for p in files[-2:]])
    write(folder/'q3d_pass_audit.json',dict(cases=audit,
        correction='Upstream caps.json relative-change denominator was 1e-9 but all-pass matrices are in F. Recomputed here from immutable raw exports; final fF capacitance values unaffected.'))
    sensitivity = []
    for length in sweep.length_um:
        changes = []
        for comp in ('cj','ref'):
            files = sorted((upstream/'raw'/f'{comp}_{length:g}um').glob('pass_*.csv'),key=lambda p:int(p.stem.split('_')[-1]))
            grounds = []
            for file in files[-2:]:
                m = pd.read_csv(file,index_col=0)
                sig = [n for n in m.index if 'ground' not in n.lower()]
                grounds.append(m.loc[sig,sig].sum(axis=1).to_numpy()*1e15)
            changes.append(abs(grounds[1]-grounds[0]))
        sensitivity.append(dict(length_um=float(length),
            excess_last_pass_bound_fF=(changes[0]+changes[1][0]/2).tolist()))
    write(folder/'subtraction_sensitivity.json',sensitivity)
    # Deterministic selection: longest requested extension, with convergence reported separately.
    row = sweep.sort_values('length_um').iloc[-1]
    width, gap = (d.parse_value(P.CPW[k]) for k in ('trace_width','trace_gap'))
    eps_eff = (1+P.d5.SI_EPS)/2
    k = width/(width+2*gap)
    z0 = 30*np.pi/np.sqrt(eps_eff)*ellipk(1-k*k)/ellipk(k*k)
    v = 299792458/np.sqrt(eps_eff)
    cprime_fit = float(np.polyfit(sweep.length_um, sweep.reference_half_fF, 1)[0]*1000)
    pl = lambda n: P.d5.f3.path_length(d,n)
    feeds = [pl('feed_in'),pl('feed_mid')+P.TEE_HALF,pl('feed_out')+P.TEE_HALF,pl('stub_p')+.1]
    full = pd.read_csv(folder/'raw/S_fine.csv',index_col=0)
    f = full.index.to_numpy(float)
    hfss = np.array([complex(str(x).replace('i','j')) for x in full['S21']])
    metrics = dict(line_lengths_mm=P.LENGTHS, feed_lengths_mm=feeds, Z0_ohm=z0, velocity_m_s=v,
        reference_slope_fF_per_mm=cprime_fit, conformal_Cprime_fF_per_mm=1e12/(z0*v),
        line_model_source='quasi-static conformal mapping, eps_eff=(1+eps_r)/2; no fitting',
        selected_length_um=float(row.length_um), cin_cp_fF=cap, models={},
        q3d_pass_audit=audit,
        subtraction_last_pass_sensitivity=sensitivity,
        limitations=['cin/cp local endpoint contributions not deembedded in either model',
            'distributed mutual coupling, tee discontinuities, bridge loading not explicitly modelled',
            'quasi-static CPW approximation; fixed lengths do not guarantee exact phase velocity',
            'raw model uses actual component pin planes; excess model uses centre planes, preserving physical layout'])
    fig, ax = plt.subplots(2,1,figsize=(10,7),sharex=True)
    ax[0].plot(f,20*np.log10(np.maximum(abs(hfss),1e-12)),label='HFSS fixed geometry',color='black')
    exported = pd.DataFrame({'frequency_GHz':f, 'HFSS_S21':hfss})
    half_span = (P.PAPER_HAMMER['gap_um']/2+P.PAPER_HAMMER['width_um']+P.d5.CJ_LEAD)/1000
    for label, values, plane in [('raw', raw_hammer, half_span),
                           ('deembedded',[row.mutual_fF,row.excess_a_fF,row.excess_b_fF],0.)]:
        s11,s21 = model(f,cap,values,P.LENGTHS,z0,v,feeds,cap_half_span=plane)
        exported[label+'_S21'] = s21
        ax[0].plot(f,20*np.log10(np.maximum(abs(s21),1e-12)),label=label)
        ax[1].plot(f,abs(s21-hfss),label=label)
        ix = find_peaks(-abs(s21))[0]
        metrics['models'][label] = dict(hammer_fF=list(map(float,values)),
            tee_to_cap_reference_mm=.2-plane,
            rms_complex_S21=float(np.sqrt(np.mean(abs(s21-hfss)**2))), dips_GHz=f[ix].tolist())
    ax[0].set_ylabel('S21 (dB)'); ax[1].set_ylabel('|model - HFSS|'); ax[1].set_xlabel('Frequency (GHz)')
    for a in ax: a.grid(alpha=.3); a.legend()
    fig.suptitle('Design4 hammer: fixed lengths, no fitted parameters')
    fig.tight_layout(); fig.savefig(folder/'fixed_length_comparison.png',dpi=160)
    exported.to_csv(folder/'model_comparison.csv',index=False)
    write(folder/'comparison.json',metrics)
    return metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('upstream',type=Path)
    args = ap.parse_args()
    upstream = args.upstream.resolve()
    if not (upstream/'cpw_sweep.csv').exists(): raise FileNotFoundError('Upstream CPW sweep not complete')
    parameters = dict(hammer=P.PAPER_HAMMER, spacing_mm=.4, lengths_mm=P.LENGTHS, trim_um=0,
        fit_parameters=[], mesh='one driven adaptive setup; wide/fine sweeps reuse it',
        max_passes=3, mesh_seed_um=10, port_impedance_ohm=50, straps='all')
    parameters['original_hammer_q3d']=dict(max_passes=20,percent_error=.1,
        reason='match the refined hammer extraction accuracy target for a fair comparison')
    folder = create_run('Design4','resonator','hammer_fixed_s21',entrypoint=__file__,
        inputs=[upstream/'cpw_sweep.csv',upstream/'run.json'],upstream=[str(upstream)],parameters=parameters)
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        import save_aedt, hfss_dirs
        import hfss_hammer_straps as E
        from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
        from scipy.signal import find_peaks
        hfss_dirs.use_d_drive(); save_aedt.use_project(save_aedt.project_file(ROOT/'Design4/simulation'))
        E.SPACING_MM=.4; E.GAP_UM=P.PAPER_HAMMER['gap_um']; E.LENGTHS=P.LENGTHS; E.NO_HAMMER=False
        d, comps = E.build(0)
        write(folder/'raw/channel_geometry.json',geometry_dump(d))
        record=json.loads((folder/'run.json').read_text(encoding='utf-8'))
        record['environment']['solver_version']=save_aedt._desktop().get_version()
        record['geometry_revision']=H.hashlib.sha256((folder/'raw/channel_geometry.json').read_bytes()).hexdigest()
        write(folder/'run.json',record)
        write(folder/'parameter_provenance.json',dict(
            geometry=dict(role='input',unit='mm and um as explicitly labelled in run.json',
                source_ref='Design4/layout/purcell_6500.py:P.LENGTHS,P.PAPER_HAMMER,P.CPW',
                status='recorded',applicability='unknown',geometry_revision=record['geometry_revision']),
            hammer_capacitance=dict(role='result',unit='fF',source_ref=str(upstream),status='recorded',
                method='Q3D Maxwell matrix, ground row sum with ground reference excluded'),
            excess_shunts=dict(role='derived',unit='fF',source_ref=str(upstream/'cpw_sweep.csv'),
                method='signal ground capacitance minus one half of matched continuous CPW capacitance',status='recorded'),
            fixed_lengths=dict(role='input',unit='mm',value=P.LENGTHS,source_ref='purcell_6500.py:LENGTHS',
                constraints='No retuning or fitting authorized in this experiment',status='recorded')))
        bare, cname = H.coupon(P.d5.CJ_LEAD)
        bare_result = H.solve_cap(bare,cname,folder/'raw/original_hammer', 'D4BARE_'+folder.name[-8:],max_passes=20,percent_error=.1)
        raw_hammer = [bare_result['mutual_fF'],*bare_result['ground_fF']]
        # cin and cp are identical local structures, so solve only one coupon.
        local = H.designs.DesignPlanar(); local.chips.main.material=P.d5.SI_MATERIAL
        local.chips.main.size.update(size_x='5.2mm',size_y='4.3mm',size_z='-500um',center_x='0mm',center_y='0mm')
        P.enclosing(local,'cp',0,0,0)
        caps = H.solve_cap(local,'cp',folder/'raw/cp_coupon','D4CP_'+folder.name[-8:])
        cap = [caps['mutual_fF'],*caps['ground_fF']]
        reuse={}
        class Driven(ScatteringImpedanceSim):
            def _analyze(self):
                P.draw_hammer_straps(self.renderer,d,which='all')
                source_name='W400g42_s21_p0um_1656_hfss'
                try:
                    project=self.renderer.pinfo.project
                    source=next(x for x in project.get_designs() if x.name==source_name)
                    previous=native_signature(source)
                    current=native_signature(self.renderer.pinfo.design)
                    write(folder/'raw/native_geometry_audit.json',dict(source_design=source_name,
                        equal=previous==current,source=previous,current=current))
                    if previous==current:
                        reuse['source']=source_name
                        print('REUSE validated historical HFSS spectrum:',source_name,flush=True)
                        return
                except Exception as exc:
                    write(folder/'raw/reuse_unavailable.json',{'reason':repr(exc)})
                E.mesh_ops(self.renderer,d,comps)
                super()._analyze()
        sim = Driven(d,'hfss')
        sim.renderer._options.update(P.WB_OPTIONS)
        sim.setup.update(freq_ghz=6.08,max_delta_s=.01,max_passes=3,min_passes=3,min_converged=1,
            pct_refinement=20,basis_order=1,vars=H.Dict())
        sim.setup.sweep_setup.update(name='Wide',start_ghz=5.3,stop_ghz=7.0,count=1701,type='Fast',save_fields=True)
        sim.run(name='D4FIX_'+folder.name[-8:],components=comps,open_terminations=[],port_list=E.PORTS)
        if reuse:
            source_folder=ROOT/'Design4/results/resonator/2026-11-12/cai_final_s21'
            shutil.copy2(source_folder/'S.csv',folder/'raw/S_wide.csv')
            source=next(x for x in sim.renderer.pinfo.project.get_designs() if x.name==reuse['source'])
            setup=source.get_setup('Setup'); before=setup.get_convergence()[0]
            before.to_csv(folder/'raw/hfss_convergence.csv')
            fine=setup.insert_sweep(start_ghz=5.85,stop_ghz=6.3,count=18001,
                name='DE_'+folder.name[-8:],type='Fast',save_fields=False)
            fine.analyze_sweep()
            freq,values=fine.get_network_data(['S11','S21'])
            pd.DataFrame(np.asarray(values).T,index=np.asarray(freq)/1e9,columns=['S11','S21']).to_csv(folder/'raw/S_fine.csv')
            after=setup.get_convergence()[0]; after.to_csv(folder/'raw/hfss_convergence_after_fine.csv')
            write(folder/'raw/hfss.json',dict(design=reuse['source'],source=str(source_folder),
                reused=True,method='native object/vertex/face-area/material/boundary signature equality',
                mesh_policy='new 25 kHz Fast sweep on existing solved setup',adaptive_history_unchanged=before.equals(after),
                raw_sha256=H.hashlib.sha256((source_folder/'S.csv').read_bytes()).hexdigest()))
            record=json.loads((folder/'run.json').read_text(encoding='utf-8'))
            record['upstream'].append(str(source_folder)); write(folder/'run.json',record)
            metrics=compare(folder,upstream,cap,raw_hammer,d)
            finish_run(folder,status='completed',summary=metrics)
            return
        _,_,S = sim.renderer.get_params(['S11','S21'])
        S.to_csv(folder/'raw/S_wide.csv')
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(folder/'raw/hfss_convergence.csv')
        f = S.index.to_numpy(float); s = S['S21'].to_numpy(complex)
        indices,_ = find_peaks(-20*np.log10(np.maximum(abs(s),1e-12)),prominence=1)
        # Include all visible features and preserve broad coverage if no dip resolves.
        lo,hi = (max(5.3,float(f[indices].min())-.1), min(7.,float(f[indices].max())+.25)) if len(indices) else (5.3,7.)
        sim.renderer.add_sweep(setup_name='Setup',name='Fine',start_ghz=lo,stop_ghz=hi,
            count=int(np.ceil((hi-lo)/.000025))+1,type='Fast',save_fields=False)
        sim.renderer.analyze_sweep('Fine','Setup')
        _,_,S = sim.renderer.get_params(['S11','S21']); S.to_csv(folder/'raw/S_fine.csv')
        after = sim.renderer.pinfo.setup.get_convergence()[0]
        after.to_csv(folder/'raw/hfss_convergence_after_fine.csv')
        if not conv.equals(after):
            raise RuntimeError('Adaptive convergence record changed during fine sweep; inspect mesh reuse')
        sim.renderer.pinfo.project._project.Save()
        write(folder/'raw/hfss.json',dict(design=sim.renderer.pinfo.design.name,
            solver_version=save_aedt._desktop().get_version(),wide_GHz=[5.3,7.],fine_GHz=[lo,hi],
            mesh_reused='Fine sweep on the identical solved Setup; no new adaptive setup',
            convergence=conv.to_dict(),nominal_pass_count_is_not_convergence=True))
        metrics = compare(folder,upstream,cap,raw_hammer,d)
        finish_run(folder,status='completed',summary=metrics)
    except Exception:
        (folder/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        finish_run(folder,status='failed',summary={'error':traceback.format_exc()})
        raise


if __name__=='__main__': main()
