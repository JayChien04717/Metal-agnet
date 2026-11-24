"""Close the AEDT session a run started, when the run ends.

Every qiskit-metal / pyEPR connection to AEDT goes through ``pyEPR.ansys.HfssApp`` (COM ``Dispatch``). If no usable
AEDT is running, that call launches a new ``ansysedt.exe``, which keeps running after Python exits. Agents that run
scripts from their own terminals (or kill Python on a timeout) leave one idle AEDT behind per run.

``install()`` (called when ``project_paths`` is imported, so every simulation script gets it) wraps
``HfssApp.__init__``. On the first connection of a process it compares the ``ansysedt.exe`` processes before and
after the connection:

* attached to an AEDT that was already open (e.g. the one you opened by hand) -> left open, nothing else happens;
* this run launched a new AEDT -> at exit (normal end or exception) all open projects are saved and closed and AEDT
  is asked to quit. A detached watchdog process waits for Python to exit (AEDT, a COM server, only quits after its
  clients release it) and force-kills that AEDT with its solver children if it is still running WATCH_GRACE_S
  later - which also covers a Python process that was hard-killed (timeout) without running its exit handlers.

Only the AEDT process this run started is ever closed. Set ``AEDT_KEEP_OPEN=1`` to keep a newly started AEDT open.

As a script:
    python aedt_session.py --list                    list running AEDT processes (PID, start time, CPU s, memory)
    python aedt_session.py --watch <ppid> <pid:t>... watchdog mode (started automatically; not for manual use)
Watchdog actions are logged to %TEMP%/aedt_session_watchdog.log.
"""

import atexit
import os
import subprocess
import sys
import tempfile
import time

KEEP_ENV = 'AEDT_KEEP_OPEN'
IMAGE = 'ansysedt.exe'
MODULE = 'pyEPR.ansys'
WATCH_GRACE_S = 90        # watchdog: wait for AEDT to quit after the parent exits, then force-kill
QUIT_WAIT_S = 60          # no watchdog: wait in-process after QuitApplication, then force-kill
LOG = os.path.join(tempfile.gettempdir(), 'aedt_session_watchdog.log')

_STATE = {'baseline': None, 'owned': {}, 'app': None, 'closed': False, 'watchdog': False}


# ---- process helpers (psutil when available, tasklist otherwise) --------------------------------------------------
def aedt_processes():
    """{pid: create_time or None} of running ansysedt.exe processes."""
    try:
        import psutil
        out = {}
        for p in psutil.process_iter(['name', 'create_time']):
            if (p.info['name'] or '').lower() == IMAGE:
                out[p.pid] = p.info['create_time']
        return out
    except ImportError:
        pass
    try:
        txt = subprocess.run(['tasklist', '/FI', f'IMAGENAME eq {IMAGE}', '/FO', 'CSV', '/NH'],
                             capture_output=True, text=True, timeout=30).stdout
    except Exception:
        return {}
    out = {}
    for line in txt.splitlines():
        cols = [c.strip('"') for c in line.split('","')]
        if len(cols) > 1 and cols[0].lower() == IMAGE and cols[1].isdigit():
            out[int(cols[1])] = None
    return out


def _same_process(pid, create_time):
    """True if pid is still the ansysedt.exe process recorded at claim time (guards against PID reuse)."""
    now = aedt_processes()
    if pid not in now:
        return False
    return create_time is None or now[pid] is None or abs(now[pid] - create_time) < 2


def _kill_tree(pid):
    subprocess.run(['taskkill', '/PID', str(pid), '/T', '/F'], capture_output=True, timeout=60)


def _wait_gone(owned, seconds):
    t_end = time.time() + seconds
    while time.time() < t_end:
        if not any(_same_process(p, t) for p, t in owned.items()):
            return True
        time.sleep(2)
    return False


def _log(msg):
    try:
        with open(LOG, 'a', encoding='utf-8') as f:
            f.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")} [{os.getpid()}] {msg}\n')
    except OSError:
        pass


# ---- claim / close -------------------------------------------------------------------------------------------------
def _claim(app):
    before = _STATE['baseline']
    after = aedt_processes()
    new = {p: t for p, t in after.items() if p not in before}
    attached = None
    try:
        attached = int(app._app.GetAppDesktop().GetProcessID())
    except Exception:
        pass
    if attached is not None:
        owned = {attached: after.get(attached)} if attached in new else {}
    else:
        owned = new
    if not owned:
        print(f'AEDT session: attached to an AEDT that was already open'
              f'{f" (PID {attached})" if attached else ""} -> left open at exit', flush=True)
        return
    _STATE['owned'] = owned
    _STATE['app'] = app
    pids = ', '.join(map(str, owned))
    if os.environ.get(KEEP_ENV, '').strip() not in ('', '0'):
        print(f'AEDT session: this run started AEDT (PID {pids}); {KEEP_ENV} is set -> left open at exit',
              flush=True)
        return
    print(f'AEDT session: this run started AEDT (PID {pids}) -> projects saved and AEDT closed at exit '
          f'(set {KEEP_ENV}=1 to keep it)', flush=True)
    atexit.register(close)
    _start_watchdog(owned)


def close():
    """Save and close all projects, quit the AEDT this run started, force-kill it if it does not exit."""
    owned = _STATE['owned']
    if _STATE['closed'] or not owned:
        return
    _STATE['closed'] = True
    try:
        desk = _STATE['app']._app.GetAppDesktop()
        for proj in list(desk.GetProjects()):
            name = proj.GetName()
            try:
                proj.Save()
            except Exception as e:
                print(f'AEDT session: could not save project {name}: {e!r}'[:200], flush=True)
            try:
                desk.CloseProject(name)
            except Exception:
                pass
        desk.QuitApplication()
        print('AEDT session: projects saved and closed, AEDT asked to quit', flush=True)
    except Exception as e:
        print(f'AEDT session: graceful quit failed ({e!r})'[:200], flush=True)
    if _STATE['watchdog']:
        # AEDT (a COM server) exits only after this process drops its COM references, i.e. after Python exits;
        # the watchdog waits for that and force-kills it if it is still running after WATCH_GRACE_S.
        return
    if _wait_gone(owned, QUIT_WAIT_S):
        print(f'AEDT session: AEDT (PID {", ".join(map(str, owned))}) closed', flush=True)
        return
    for pid, t in owned.items():
        if _same_process(pid, t):
            _kill_tree(pid)
            print(f'AEDT session: AEDT PID {pid} did not quit -> killed', flush=True)


def _start_watchdog(owned):
    args = [sys.executable, '-I', os.path.abspath(__file__), '--watch', str(os.getpid())]
    args += [f'{p}:{t if t is not None else ""}' for p, t in owned.items()]
    base = 0x00000008 | 0x00000200          # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    kw = dict(stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, close_fds=True)
    for flags in (base | 0x01000000, base):   # try CREATE_BREAKAWAY_FROM_JOB first (outlive a kill-on-close job)
        try:
            subprocess.Popen(args, creationflags=flags, **kw)
            _STATE['watchdog'] = True
            return
        except OSError:
            continue
    print('AEDT session: could not start the watchdog; AEDT is closed only on a normal exit', flush=True)


def _watch(ppid, owned):
    import ctypes
    k32 = ctypes.windll.kernel32
    k32.OpenProcess.restype = ctypes.c_void_p
    h = k32.OpenProcess(0x00100000, False, ppid)          # SYNCHRONIZE
    _log(f'watching parent {ppid} for AEDT {sorted(owned)}')
    if h:
        k32.WaitForSingleObject(ctypes.c_void_p(h), 0xFFFFFFFF)
        k32.CloseHandle(ctypes.c_void_p(h))
    if _wait_gone(owned, WATCH_GRACE_S):
        _log(f'parent {ppid} exited; AEDT {sorted(owned)} quit')
        return
    for pid, t in owned.items():
        if _same_process(pid, t):
            _kill_tree(pid)
            _log(f'parent {ppid} exited; AEDT PID {pid} still running after {WATCH_GRACE_S} s -> killed '
                 '(projects may keep a .lock file)')


# ---- install -------------------------------------------------------------------------------------------------------
def _patch(module):
    cls = module.HfssApp
    if getattr(cls.__init__, '_aedt_session', False):
        return
    original = cls.__init__

    def __init__(self, *args, **kwargs):
        first = _STATE['baseline'] is None
        if first:
            _STATE['baseline'] = aedt_processes()
        original(self, *args, **kwargs)
        if first:
            try:
                _claim(self)
            except Exception as e:
                print(f'AEDT session: tracking failed ({e!r}); AEDT will not be closed automatically'[:200],
                      flush=True)
    __init__._aedt_session = True
    cls.__init__ = __init__


class _Hook:
    """Meta-path finder that patches pyEPR.ansys right after the module is executed."""

    def find_spec(self, name, path=None, target=None):
        if name != MODULE:
            return None
        import importlib.machinery
        spec = importlib.machinery.PathFinder.find_spec(name, path)
        if spec is None or spec.loader is None:
            return spec
        exec_module = spec.loader.exec_module

        def exec_and_patch(module):
            exec_module(module)
            _patch(module)
        spec.loader.exec_module = exec_and_patch
        return spec


def install():
    if sys.platform != 'win32':
        return
    if MODULE in sys.modules:
        _patch(sys.modules[MODULE])
    elif not any(isinstance(f, _Hook) for f in sys.meta_path):
        sys.meta_path.insert(0, _Hook())


def _list():
    procs = aedt_processes()
    if not procs:
        print('no AEDT (ansysedt.exe) running')
        return
    try:
        import psutil
    except ImportError:
        psutil = None
    for pid, t in sorted(procs.items()):
        line = f'PID {pid}'
        if t:
            line += f'  started {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t))}'
        if psutil:
            try:
                p = psutil.Process(pid)
                line += f'  CPU {sum(p.cpu_times()[:2]):.0f} s  mem {p.memory_info().rss / 2**20:.0f} MB'
            except psutil.Error:
                pass
        print(line)


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--watch':
        owned = {}
        for item in sys.argv[3:]:
            pid, _, t = item.partition(':')
            owned[int(pid)] = float(t) if t else None
        _watch(int(sys.argv[2]), owned)
    else:
        _list()
