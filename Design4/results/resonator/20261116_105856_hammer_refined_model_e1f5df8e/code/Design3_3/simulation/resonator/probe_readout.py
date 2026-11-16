"""Identify solved readout eigenmodes by CPW-slot electric field."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import os

import numpy as np
import pandas as pd
from shapely.geometry import Point

from q3d import base
from layout import build_design
from recheck_readout_flux import PROJECT, RESULTS, READOUTS


def sample_locations(trace):
    rows = []
    for frac in np.linspace(.08, .92, 25):
        s = frac*trace.length
        mid = np.array(trace.interpolate(s).coords[0])
        d = np.array(trace.interpolate(s+.003).coords[0])-np.array(trace.interpolate(s-.003).coords[0])
        n = np.array([-d[1],d[0]])/np.linalg.norm(d)
        for side in (-1, 1):
            near = mid+n*side*.008
            far = mid+n*side*.032
            if trace.distance(Point(far)) < .023:
                continue
            rows += [dict(s_fraction=frac, side=side, region=region, x_mm=p[0], y_mm=p[1], z_mm=.002)
                     for region,p in [('near',near),('far',far)]]
    return pd.DataFrame(rows)


def probe(project, layout, device):
    design = project.get_design(f'readout_{device}_hfss')
    setup = design.get_setup('Setup')
    sol = setup.get_solutions()
    native = design._design
    fr = native.GetModule('FieldsReporter')
    table = layout.qgeometry.tables['path']
    trace = table.loc[table.component.eq(layout.components[f'res_{device}'].id) & table.name.eq('trace'),
                      'geometry'].item()
    points = sample_locations(trace)
    points_file = PROJECT.parent/f'{device}.pts'
    coordinates_mm = points[['x_mm','y_mm','z_mm']].to_numpy()
    with points_file.open('w', encoding='ascii') as stream:
        stream.write('Unit=mm\n')
        np.savetxt(stream, coordinates_mm)
    frequencies, _ = sol.eigenmodes()
    rows = []
    for mode, frequency in enumerate(frequencies,1):
        sol.set_mode(mode)
        fr.CalcStack('Clear')
        fr.EnterQty('E')
        fr.CalcOp('CmplxMag')
        fr.CalcOp('Mag')
        output = PROJECT.parent/f'fields_{device}_{mode}.fld'
        fr.ExportToFile(str(output), str(points_file), setup.solution_name, ['Phase:=','0deg'], True)
        fields = np.loadtxt(output,skiprows=1)
        if (fields.shape != (len(points),4) or not np.isfinite(fields).all()
                or not np.allclose(fields[:,:3], coordinates_mm*1e-3, atol=1e-12)):
            raise ValueError(f'Unexpected field sample shape for {device} mode {mode}')
        row = points.assign(mode=mode,frequency_GHz=frequency,E=fields[:,-1])
        rows.append(row)
    frame=pd.concat(rows,ignore_index=True)
    frame.to_csv(_relocate(RESULTS/f'readout_{device}_slots.csv'),index=False)
    print(frame.assign(E2=frame.E**2).groupby(['mode','frequency_GHz','region']).E2.mean())


def main():
    os.chdir(RESULTS)
    project=base.save_aedt._desktop().get_active_project()
    if project.name != PROJECT.stem:
        base.save_aedt.use_project(PROJECT)
        project=base.save_aedt._desktop().get_active_project()
    design=build_design(for_process=True)
    for device in READOUTS:
        probe(project,design,device)


if __name__ == '__main__':
    main()
