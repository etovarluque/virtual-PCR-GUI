"""Entry point for `python -m vpcr` and for Run vPCR.bat.

Run vPCR.bat keeps its console open behind the window, so diagnostics go there.
Nothing here assumes a console exists, though: `pythonw -m vpcr` works too, and
then sys.stdout is None and messages fall back to a native dialog.
"""
import io
import subprocess
import sys
import traceback
from datetime import datetime
from pathlib import Path

from vpcr import APP_ROOT

REQUIRED = ("PyQt6", "pandas", "openpyxl")

_MB_ICONERROR = 0x10


def _ensure_streams() -> bool:
    """Give the process usable streams. Returns whether a real console exists."""
    console = sys.__stdout__ is not None
    for name in ("stdout", "stderr"):
        if getattr(sys, name, None) is None:
            setattr(sys, name, io.StringIO())
    return console


def _message(text: str, title: str = "virtualPCR Studio") -> None:
    """Report to the console when there is one; otherwise to a native dialog."""
    if sys.__stdout__ is not None:
        print(text, file=sys.stderr)
        return
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(None, text, title, _MB_ICONERROR)
    except Exception:
        pass


def _missing() -> list[str]:
    import importlib.util
    return [m for m in REQUIRED if importlib.util.find_spec(m) is None]


def _install_dependencies(missing: list[str]) -> bool:
    req = APP_ROOT / "requirements.txt"
    print(f"Missing: {', '.join(missing)}. Installing from {req.name}...")
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
             "-r", str(req)])
    except OSError as e:
        _message(f"Could not run pip: {e}")
        return False
    if proc.returncode != 0:
        _message(f"Installing the dependencies failed.\n"
                 f"Try manually:  {sys.executable} -m pip install -r {req}")
        return False
    return not _missing()


def _report_crash(exc_type, exc, tb) -> None:
    text = "".join(traceback.format_exception(exc_type, exc, tb))
    log = APP_ROOT / "crash.log"
    try:
        with open(log, "a", encoding="utf-8") as fh:
            fh.write(f"\n{'=' * 60}\n{datetime.now():%Y-%m-%d %H:%M:%S}\n{text}")
    except OSError:
        log = None

    # The console already shows the traceback for us; the dialog matters when
    # the user is looking at the window, not at the console behind it.
    print(text, file=sys.stderr)
    try:
        from PyQt6.QtWidgets import QApplication, QMessageBox
        if QApplication.instance() is not None:
            where = f"\n\nWritten to {log}" if log else ""
            QMessageBox.critical(None, "virtualPCR Studio",
                                 f"{exc_type.__name__}: {exc}{where}")
    except Exception:
        pass


def main() -> int:
    _ensure_streams()
    sys.excepthook = _report_crash

    if (missing := _missing()) and not _install_dependencies(missing):
        return 1

    try:
        from vpcr.gui.main_window import main as gui_main
    except Exception as e:
        _report_crash(type(e), e, e.__traceback__)
        return 1
    return gui_main()


if __name__ == "__main__":
    raise SystemExit(main())
