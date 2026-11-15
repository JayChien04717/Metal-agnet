"""Print the AEDT message window (errors / warnings) of the open Design6 project: python aedt_messages.py [design]"""
import sys
from pyEPR.ansys import HfssApp

desk = HfssApp().get_app_desktop()._desktop
design = sys.argv[1] if len(sys.argv) > 1 else ''
for m in desk.GetMessages('Design6', design, 0)[-40:]:
    print(m)
