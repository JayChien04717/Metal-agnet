"""Put the AEDT temp directory (where HFSS writes the large hfss_*.pjt solver scratch, tens of GB per run) on D:
instead of C:\\Users\\...\\AppData\\Local\\Temp. Called by hfss_chan5.py / hfss_chip6.py before every run.

    use_d_drive()             -> sets AEDT temp dir to D:/HFSS_temp (env HFSS_TEMP overrides), prints before / after
    python hfss_dirs.py       -> read-only: prints the current AEDT temp and project directories

Note: SetTempDirectory changes the AEDT application option (it stays set for later sessions too). It does not move a
solve that is already running; it applies to the next Analyze.
"""
import os

DEFAULT_TEMP = 'D:/HFSS_temp'


def _desktop():
    from pyEPR.ansys import HfssApp
    return HfssApp().get_app_desktop()            # connects to the running AEDT (starts one if none)


def use_d_drive(temp=None):
    temp = (temp or os.environ.get('HFSS_TEMP') or DEFAULT_TEMP).replace('\\', '/')
    os.makedirs(temp, exist_ok=True)
    desk = _desktop()
    before = desk.temp_directory
    desk.temp_directory = temp
    after = desk.temp_directory
    ok = os.path.normcase(os.path.normpath(after)) == os.path.normcase(os.path.normpath(temp))
    print(f'AEDT temp directory: {before} -> {after}  ({"OK" if ok else "NOT CHANGED"})', flush=True)
    if not ok:
        raise RuntimeError(f'AEDT temp directory is still {after}, wanted {temp}')
    return dict(before=before, after=after)


if __name__ == '__main__':
    d = _desktop()
    print('AEDT version  :', d.get_version())
    print('temp directory:', d.temp_directory)
    print('project dir   :', d.project_directory)
