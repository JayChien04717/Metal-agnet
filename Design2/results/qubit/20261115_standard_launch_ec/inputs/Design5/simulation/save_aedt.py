"""AEDT project handling for the HFSS / Q3D run scripts.

Rule (Jay 2026-11-09): every run of a Design goes into ONE AEDT session and ONE project file,
<Design>/simulation/aedt/<Design>.aedt - nothing is closed or reopened between runs; each run adds its own
design (name with the run time, so reruns of the same tag do not collide).

  use_project(path)            open the project if it is not open yet (create + save it if the file does not exist)
                               and make it the active project, so qiskit-metal renders into it. Call before sim.run.
  snapshot(folder, name)       Save the project in place and copy the .aedt (geometry, setups, mesh operations;
                               not the large .aedtresults) to <folder>/<name>.aedt -> kept with the run folder.
  design_name(tag)             '<tag>_<HHMM>' (qiskit adds '_hfss' / '_q3d')
  save_active(folder, name)    old behaviour (SaveAs a copy with results, optionally close); kept for older runs.
As a script:  python save_aedt.py <run folder> [name]   (save_active without closing)
"""
import os, sys, time, shutil


def _desktop():
    from pyEPR.ansys import HfssApp
    return HfssApp().get_app_desktop()


def project_file(sim_dir):
    """<Design>/simulation/aedt/<Design>.aedt for a simulation folder (env AEDT_PROJECT overrides it, e.g. to keep
    working in a project that already holds a design to link to)."""
    if os.environ.get('AEDT_PROJECT'):
        return os.environ['AEDT_PROJECT']
    design = os.path.basename(os.path.dirname(os.path.abspath(sim_dir)))
    return os.path.join(os.path.abspath(sim_dir), 'aedt', f'{design}.aedt')


def design_name(tag):
    return f'{tag}_{time.strftime("%H%M")}'


def use_project(path):
    path = os.path.abspath(path).replace('\\', '/'); name = os.path.splitext(os.path.basename(path))[0]
    desk = _desktop(); D = desk._desktop
    open_names = [p.name for p in desk.get_projects()]
    if name not in open_names:
        if os.path.exists(path):
            D.OpenProject(path)
        else:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            pr = D.NewProject(); pr.SaveAs(path, True)
    D.SetActiveProject(name)
    print(f'AEDT project: {path} (open: {[p.name for p in desk.get_projects()]})', flush=True)
    ensure_si_material(desk.get_active_project()._project)
    return name


def ensure_si_material(proj):
    """Substrate material for the layouts (design5_layout.SI_MATERIAL, env SI_EPS, default 11.45): a copy of the
    library 'silicon' with only the permittivity changed (the loss / conductivity of the library material are kept).
    Made once per project; the layouts name it 'silicon_eps11p45'."""
    eps = float(os.environ.get('SI_EPS', 11.45))
    if abs(eps - 11.9) < 1e-9:
        return 'silicon'
    mat = 'silicon_eps' + f'{eps:g}'.replace('.', 'p')
    dm = proj.GetDefinitionManager()
    if not dm.DoesMaterialExist(mat):
        dm.CloneMaterial('silicon', mat)
        dm.EditMaterial(mat, [f'NAME:{mat}', 'permittivity:=', f'{eps:g}'])
    print(f'substrate material {mat}: permittivity {eps:g}', flush=True)
    return mat


def snapshot(folder, name):
    proj = _desktop().get_active_project()
    proj._project.Save()
    src = os.path.join(proj._project.GetPath(), proj._project.GetName() + '.aedt')
    os.makedirs(folder, exist_ok=True)
    dst = os.path.join(folder, f'{name}.aedt'); shutil.copy(src, dst)
    print('project saved in place:', src, '| snapshot:', dst, flush=True)
    return dst


def save_active(folder, name, close=True):
    desk = _desktop()
    for p in desk.get_projects():
        print('open project:', p.name, '| designs:', [d.name for d in p.get_designs()])
    proj = desk.get_active_project()
    os.makedirs(folder, exist_ok=True)
    target = os.path.join(os.path.abspath(folder), f'{name}.aedt').replace('\\', '/')
    proj._project.SaveAs(target, True)
    print('saved', target)
    if close:
        try:
            desk._desktop.CloseProject(name)
        except Exception as e:
            print('project left open:', repr(e)[:120])
    return target


if __name__ == '__main__':
    run_dir = sys.argv[1]
    save_active(run_dir, sys.argv[2] if len(sys.argv) > 2 else os.path.basename(run_dir), close=False)
