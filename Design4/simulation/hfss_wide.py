"""Wide Interpolating S21 of the combined chain: Purcell pass-band (6.4-7.0 GHz)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hfss_combined as hc  # noqa: E402

hc.SETUP.update(freq_ghz=6.65, max_delta_s=0.01, max_passes=10)
hc.SWEEP.update(start_ghz=6.4, stop_ghz=7.0, count=601, type='Interpolating')
hc.run('wide3')
