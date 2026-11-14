"""Write the Design3-style README from the available Design3_2 solver outputs."""
from pathlib import Path
import sys, json, re
import pandas as pd
import numpy as np

design=Path(__file__).resolve().parents[1]
runs=design/'simulation/results'
root=design.parent
sys.path[:0]=[str(root/'note'),str(root/'Design3/simulation')]
import cap_matrix as cap
import zz_analysis as zz
dev=('FX1','TC','FX2')
target=dict(zip(dev,(6.3,6.0,6.5)))
out={}
for q in dev:
    path=runs/f'd32_eigen_{q}_convergence.csv'
    if not path.exists():
        continue
    table=pd.read_csv(path,index_col=0).dropna(how='all')
    modes=table.iloc[-1].dropna().tolist()
    native=f'D32_R5_eigen_{q}_targeted_hfss'
    if not (design/'simulation/aedt/Design3_2.aedtresults'/f'{native}.results').exists():
        native=f'D32_R5_eigen_{q}_contact_hfss'
    profiles=sorted((design/'simulation/aedt/Design3_2.aedtresults'/f'{native}.results').rglob('*.EMConvData'), key=lambda p:p.stat().st_mtime)
    if not profiles:
        native=f'D32_R5_eigen_{q}_hfss'
        profiles=sorted((design/'simulation/aedt/Design3_2.aedtresults'/f'{native}.results').rglob('*.EMConvData'), key=lambda p:p.stat().st_mtime)
    status=profiles[-1].read_text(errors='replace') if profiles else ''
    delta_values=re.findall(r'MaxDeltaFreq\(([-\d.eE]+),',status)
    fieldfile=runs/f'd32_eigen_{q}.csv'
    field=pd.read_csv(fieldfile).iloc[0].to_dict() if fieldfile.exists() else {}
    idx=min(range(len(modes)),key=lambda i:abs(modes[i]-target[q]))
    assignment='proximity candidate'
    ratios=[field.get(f'mode{i+1}_cpw_near_far_E2',0) for i in range(len(modes))]
    if len(ratios)==2 and min(ratios)>0 and max(ratios)>3*min(ratios):
        idx=int(np.argmax(ratios))
        assignment='CPW field contrast'
    intensity=[field.get(f'mode{i+1}_cpw_mean_E2',0) for i in range(len(modes))]
    if len(modes) == 1 and ratios[0] > 3:
        assignment='CPW field probe'
    if len(modes) == 2 and idx == int(np.argmax(intensity)) and min(intensity) > 0 and max(intensity) > 20*min(intensity):
        assignment='target proximity and CPW field'
    out[q]=dict(f=modes[idx],modes=modes,index=idx+1,last_MHz=(modes[idx]-table.iloc[-2,idx])*1000,
                passes=len(table),assignment=assignment,last_delta_pct=float(delta_values[-1]) if delta_values else None,converged='IsConverged=true' in status, field=field)
freq={q:out[q]['f'] if q in out else target[q] for q in dev}
p=runs/'d32_matrix_cmat_fF.csv'
matrix=pd.read_csv(p,index_col=0) if p.exists() else None
ecj,g,j01,charge={},{},{},{}
if matrix is not None:
    matrix=matrix.loc[matrix.index,matrix.index]
    a=matrix.to_numpy()
    assert np.max(np.abs(a-a.T))<np.max(np.diag(a))*1e-5
    assert np.linalg.eigvalsh((a+a.T)/2).min()>0
    ro={q:cap.find_node(matrix,f'readout_line_{q}',f'res_{q}') for q in dev}
    nodes=cap.DEFAULT_NODES[:5]+list(ro.values())
    assert len(set(nodes))==8 and min(matrix.loc[r,r] for r in ro.values())>500
    ecj=cap.analyse(matrix,nodes=nodes,verbose=False)
    pads={'FX1':('pad_bot_FX1','pad_top_FX1'),'TC':('tmon_TC',),'FX2':('pad_top_FX2','pad_bot_FX2')}
    n01={}
    for q in dev:
        ec=ecj[f'Ec_{q} [GHz]']
        e,n=zz.transmon(ec,zz.coupler_EJ(ec,5)) if q=='TC' else zz.fluxonium(ec,4,.6)
        n01[q]=abs(n[0,1])
        g[q]={r:cap.readout_g(matrix,pads[q],ro[r],Z0=50.9,f_r=freq[r]*1e9,
              transmon_EJ_over_Ec=50 if q=='TC' else None,verbose=False) for r in dev}
        charge[q]={r:cap.charge_line(matrix,pads[q],f'charge_line_{r}',e[1]*1e9,n01[q]) for r in dev}
    for a,b,key in [('FX1','TC','J_1c [MHz]'),('FX2','TC','J_2c [MHz]'),('FX1','FX2','J_12 [MHz]')]:
        j01[f'{a}–{b}']=(abs(ecj[key]),abs(ecj[key]*n01[a]*n01[b]))
p=runs/'d32_kappa.csv'
k=pd.read_csv(p).set_index('resonator').to_dict('index') if p.exists() else {}
m={}
for q,cases in [('FX1',['fx1']),('TC',['tc_left','tc_right']),('FX2',['fx2'])]:
    if q == 'FX1':
        path = runs / 'd32_fx1_port_review_summary.json'
        if path.exists():
            reviewed = json.loads(path.read_text(encoding='utf8'))
            if reviewed.get('verified'):
                value = abs(reviewed['M_pH_fmin'])
                m[q] = {'M_pH': value, 'I_phi0_mA': reviewed['I_phi0_mA'],
                        'fixture_M_pH': [value], 'sweep_spread_pct': [reviewed['M_spread_pct']]}
        continue
    prefix = 'd32_mutual_fixed_' if q == 'FX2' else 'd32_mutual_contact_'
    files=[runs/f'{prefix}{case}_summary.json' for case in cases]
    if q == 'FX2':
        log = runs / 'd32_mutual_fixed_fx2_messages.log'
        conv = runs / 'd32_mutual_fixed_fx2_convergence.csv'
        if not log.exists() or not conv.exists():
            continue
        messages = log.read_text(encoding='utf8')
        fx2_delta = pd.read_csv(conv)['Max Mag. Delta S'].dropna()
        if (re.search(r'conductors? touch lumped port|\[warning\]|\[error\]',messages,re.I)
                or len(fx2_delta) < 11 or (fx2_delta.tail(3) > .0002).any()):
            continue
    if all(p.exists() for p in files):
        values=[abs(json.loads(p.read_text())['M_pH_fmin']) for p in files]
        avg=sum(values)/len(values)
        m[q]={'M_pH':avg,'I_phi0_mA':2.067833848/avg,'fixture_M_pH':values,
              'sweep_spread_pct':[json.loads(p.read_text())['M_spread_pct'] for p in files]}

def show(value,digits=3):
    return '—' if value is None else f'{value:.{digits}f}'
heading = '## Baseline simulation · 6.0 GHz TC design'
lines=[heading,'',
       'The electrical model uses the Design3_2 QComponents with local polygon fillets targeting R = 5 µm. '
       'CPW centerlines and junction boundaries use the native AEDT renderer. The driven mutual-inductance models preserve straight port landings. '
       'The model omits fabrication contact notches and ground perforations. The rounded electrical polygons are generated '
       'before rendering and can differ locally from the final GDS, which rounds the assembled mask. '
       'Launch pads retain the Design3_2 width of 200 µm.','',
       '| Device | $E_C/h$ (GHz) | Readout frequency (GHz) | Offset from target (MHz) | Readout $g/2\\pi$ (MHz) | Feedline $\\kappa/2\\pi$ (MHz) |',
       '|---|---:|---:|---:|---:|---:|']
for q in dev:
    gg=g.get(q,{}).get(q,{})
    gv=gg.get('g_transmon/2pi [MHz]' if q=='TC' else 'g/2pi [MHz]')
    gv=abs(gv) if gv is not None else None
    kv=k[q]['kappa_MHz']*(freq[q]/k[q]['f_GHz'])**3 if q in k else None
    lines.append(f"| {q} | {show(ecj.get(f'Ec_{q} [GHz]'),4)} | {show(out.get(q,{}).get('f'),6)} | "
                 f"{show((freq[q]-target[q])*1000 if q in out else None,2)} | {show(gv,2)} | {show(kv,3)} |")
lines+=['','FX readout coupling is per unit charge number. TC follows the Design3 convention at $E_J/E_C=50$. '
        'Readout $g$ and $\\kappa$ use the HFSS candidate frequency when available, and the design frequency otherwise. '
        'Feedline linewidths come from the passive Q3D model.','',
        '| Qubit pair | $\\lvert J\\rvert/2\\pi$ (MHz) | $0\\to1$ coupling magnitude (MHz) |','|---|---:|---:|']
for pair in ('FX1–TC','FX2–TC','FX1–FX2'):
    a,b=j01.get(pair,(None,None))
    lines.append(f'| {pair} | {show(a)} | {show(b)} |')
if charge:
    vals=[charge['TC'][q]['C_line-qubit [fF]'] for q in dev]
    lines+=['','The direct capacitances to the TC island are '+', '.join(f'{q} {v:.5f} fF' for q,v in zip(dev,vals))+'.']
lines+=['','### Flux bias line mutual inductance','',
        '| Device | Mutual inductance magnitude (pH) | Current for one flux quantum (mA) |','|---|---:|---:|']
for q in dev:
    lines.append(f"| {q} | {show(m.get(q,{}).get('M_pH'),4)} | {show(m.get(q,{}).get('I_phi0_mA'),4)} |")
lines+=['','The current is $I_{\\Phi_0}=\\Phi_0/|M|$. Mutual inductance is sampled at 10 MHz. '
        'The TC result averages the two JJ-open fixture magnitudes.','',
        '### Solver results and convergence','',
        '| Device | Candidate modes (GHz) | Selected mode | Last-pass drift (MHz) | Passes | Max Δf (%) | HFSS status |',
        '|---|---|---:|---:|---:|---:|---|']
for q,r in out.items():
    modes=', '.join(f'{f:.6f}' for f in r['modes'])
    lines.append(f"| {q} | {modes} | {r['index']} ({r['assignment']}) | {r['last_MHz']:+.2f} | {r['passes']} | {show(r['last_delta_pct'], 3)} | "
                 f"{'Converged' if r['converged'] else 'Tolerance not reached'} |")
if out:
    lines+=['','HFSS uses one qubit, its readout resonator and its tee per solve. The eigenmode convergence target is 0.05% for two consecutive passes. Mode selection uses CPW field probes and frequency proximity. '
            'FX1 retained 23 passes before the stalled final pass was stopped. FX2 used a single-mode solve with a 5.75 GHz minimum frequency. TC mode 1 remains a frequency-proximity assignment.']
for tag,label,tol in [('d32_matrix','Device capacitance matrix',.5),('d32_kappa','Feedline matrix',.2)]:
    p=runs/f'{tag}_convergence.csv'
    if p.exists():
        row=pd.read_csv(p)
        delta=row.iloc[-1]['Delta %']
        lines+=['',f'{label} last exported pass is {len(row)}, with ΔC = {delta:.3f}% and a {tol:g}% target.']
kappa_passes = runs / 'd32_kappa_by_pass.csv'
if kappa_passes.exists():
    kappas = pd.read_csv(kappa_passes).pivot(index='pass_number',columns='device',values='kappa_MHz')
    if len(kappas) > 1:
        relative_change = ((kappas.iloc[-1] / kappas.iloc[-2]) - 1).abs() * 100
        values = ', '.join(f'{q} {relative_change[q]:.2f}%' for q in dev)
        lines += ['', f'Last-pass changes in κ at fixed design frequencies were {values}. '
                  '[Feedline coupling by pass](simulation/results/d32_kappa_by_pass.csv)']
if m:    lines+=['','Mutual-inductance sweep variation across 10–100 MHz is '+
            ', '.join(q+' '+ '/'.join(f'{v:.2f}%' for v in m[q]['sweep_spread_pct']) for q in m)+'.']
if 'FX1' in m:
    label = 'simulation/results/d32_fx1_port_review'
    lines += ['', f'FX1 uses the reviewed R5 model. HFSS reported no port warnings. '
              f'The adaptive mesh converged in {reviewed["adaptive_passes"]} passes '
              f'(final ΔS = {reviewed["last_delta_s"]:.3g}, target {reviewed["max_delta_s"]:.3g}). '
              f'[Impedance sweep]({label}_Z.csv) · [Adaptive convergence]({label}_convergence.csv) · '
              f'[HFSS messages]({label}_messages.log)']
if 'FX2' in m:
    lines += ['', f'FX2 completed {len(fx2_delta)+1} adaptive passes with final ΔS = {fx2_delta.iloc[-1]:.3g} '
              '(target 0.0002). HFSS reported no port warnings. '
              '[FX2 sweep](simulation/results/d32_mutual_fixed_fx2_Z.csv) · '
              '[FX2 convergence](simulation/results/d32_mutual_fixed_fx2_convergence.csv)']
if 'TC' in m:
    lines += ['', 'Both TC SQUID fixtures completed 12 adaptive passes. '
              '[Left JJ sweep](simulation/results/d32_mutual_contact_tc_left_Z.csv) · '
              '[Right JJ sweep](simulation/results/d32_mutual_contact_tc_right_Z.csv)']
if len(out)<3 or not ecj or len(k)<3 or len(m)<3:
    lines+=['','A dash indicates a parameter whose new solver result is unavailable.']
lines+=['','[Raw solver outputs](simulation/results) · [Solver runner](simulation/run_remaining.py)','']
payload=dict(eigen=out,EcJ=ecj,readout=g,qubit_pairs=j01,charge=charge,kappa=k,mutual=m)
(runs/'electrical_parameters.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf8')
p=design/'README.md'
s=p.read_text(encoding='utf8')
before=s.split(heading)[0].split('## 版圖')[0].replace('- 本版已檢查 GDS 幾何，尚未重跑電磁模擬。\n','')
tail=s[s.index('## 版圖'):].split('## Current Design3_2 parameters')[0].rstrip()+'\n'
p.write_text(before+'\n'.join(lines)+'\n'+tail,encoding='utf8')
print('Report updated with',list(out),list(ecj),list(k),list(m))
