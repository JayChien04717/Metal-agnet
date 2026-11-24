"""Print the AEDT message window (errors / warnings) of the open Design6_TXT project: python aedt_messages.py [design]"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import sys
from pyEPR.ansys import HfssApp

desk = HfssApp().get_app_desktop()._desktop
design = sys.argv[1] if len(sys.argv) > 1 else ''
for m in desk.GetMessages('Design6_TXT', design, 0)[-40:]:
    print(m)
