"""Save the project open in the running AEDT (the last HFSS / Q3D run), with its solved results (.aedtresults).

As a function (called by the HFSS scripts right after the sweep):  save_active(folder, name, close=True)
  -> <folder>/<name>.aedt + <name>.aedtresults; close=True closes the project afterwards so the files can be moved.
As a script:  python save_aedt.py <run folder> [name]      e.g. results/2026-11-09/1134_chip6_v1 chip6_v1
"""
import os, sys


def save_active(folder, name, close=True):
    from pyEPR.ansys import HfssApp
    desk = HfssApp().get_app_desktop()
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
