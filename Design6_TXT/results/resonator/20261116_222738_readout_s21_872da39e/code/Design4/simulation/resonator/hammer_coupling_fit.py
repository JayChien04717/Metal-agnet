"""Extract an effective two-mode J from recorded complex HFSS S21, not dip splitting."""
from pathlib import Path
import argparse
import json
import traceback
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import hammer_fixed_s21 as F
import modes5 as M
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    args=ap.parse_args();source=args.source.resolve()
    folder=create_run('Design4','resonator','hammer_coupling_fit',entrypoint=__file__,
        inputs=[source/'raw/S_fine.csv',source/'run.json'],upstream=[str(source)],
        parameters=dict(method='modes5.direct_fit complex two-mode input/output model',
            geometry_changed=False,main_lengths_changed=False,units_J='MHz (Hamiltonian divided by h)',
            model_assumptions='one lossy Purcell mode, negligible intrinsic readout loss; smooth feed background'))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        raw=pd.read_csv(source/'raw/S_fine.csv',index_col=0)
        f=raw.index.to_numpy(float);z=raw.S21.map(complex).to_numpy()
        peaks,prop=find_peaks(-20*np.log10(abs(z)),prominence=3)
        assert len(peaks)>=2, 'Two resolved dips required'
        peaks=peaks[np.argsort(prop['prominences'])[-2:]]
        mask=np.zeros(len(f),bool);mask[::10]=True
        for i in peaks:mask|=abs(f-f[i])<.002
        fit,x=M.direct_fit(f[mask],z[mask])
        prediction=np.conj(M.io_model(f,x))
        result=dict(source=str(source),fit=fit,model_parameters=x.tolist(),
            full_grid_rms=float(np.sqrt(np.mean(abs(prediction-z)**2))),
            fit_sample_count=int(mask.sum()),raw_sample_count=len(f),
            sampling='every tenth background point plus every point within 2 MHz of both dips',
            coupling_interpretation='effective two-mode fit J; not independent EM coupling or half dip splitting',
            validation='exploratory: HFSS mesh/direct-frequency and model-identifiability limitations retained')
        # Window sensitivity, with all narrow-feature samples still present.
        fit2,x2=M.direct_fit(f[mask],z[mask],window=.15)
        result['narrow_window_fit']=fit2
        result['window_relative_J_difference']=abs(fit2['J']-fit['J'])/max(abs(fit['J']),1e-12)
        F.write(folder/'coupling.json',result)
        pd.DataFrame({'frequency_GHz':f,'HFSS_S21':z,'fitted_S21':prediction}).to_csv(folder/'fit.csv',index=False)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)
        ax[0].plot(f,20*np.log10(abs(z)),label='HFSS')
        ax[0].plot(f,20*np.log10(abs(prediction)),label=f"effective two-mode fit: J={fit['J']:.3f} MHz")
        ax[1].plot(f,abs(prediction-z));ax[0].set_ylabel('S21 (dB)');ax[1].set_ylabel('complex residual')
        ax[1].set_xlabel('Frequency (GHz)');ax[0].legend()
        for a in ax:a.grid(alpha=.25)
        fig.tight_layout();fig.savefig(folder/'coupling_fit.png',dpi=150)
        finish_run(folder,status='completed',summary=result)
        print(json.dumps({'J_MHz':fit['J'],'rms_rel':fit['rms_rel'],'window_relative_J_difference':result['window_relative_J_difference']}),flush=True)
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error)
        finish_run(folder,status='failed',summary={'error':error});raise


if __name__=='__main__':main()
