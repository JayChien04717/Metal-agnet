"""Design3 flux-line mutual inductance: converged 10-100 MHz run + broadband run, then the figure."""
from qiskit_metal import Dict
import hfss_mutual_fx as m
import plot_mutual_z

m.run('mutual_d3')
m.run('mutual_d3_broadband', setup=dict(sweep_setup=Dict(name='Sweep', start_ghz=0.01, stop_ghz=6.0, count=300,
                                                         step_ghz=None, type='Interpolating', save_fields=False)))
plot_mutual_z.main('mutual_d3_broadband', 'mutual_d3')
