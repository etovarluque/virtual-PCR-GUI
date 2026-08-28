"""Loading a primer file with the user in the loop.

Both the Primers panel and the Manage primers dialog load files, and both must
ask the same questions on the way: which worksheet, and which columns. Keeping
that flow here means the two cannot drift apart.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import pandas as pd
from PyQt6.QtWidgets import (QDialog, QFileDialog, QInputDialog, QMessageBox,
                             QWidget)

from vpcr.core.model import PrimerPair
from vpcr.core.primers import (PrimerConflicts, PrimerFileError,
                               find_primer_conflicts, is_fasta, load_pairs,
                               primer_set_from_pairs, read_table)
from vpcr.gui.column_dialog import ColumnMapDialog

PRIMER_FILTER = "Primer files (*.fasta *.fa *.txt *.tsv *.csv *.xlsx);;All files (*)"


def _format_conflicts(c: PrimerConflicts, limit: int = 12) -> str:
    """A human-readable summary of the duplicate/clashing primers found."""
    lines = ["This file has primers that look like duplicates or clashes:"]
    if c.same_seq_diff_name:
        lines.append("\nSame sequence, different names:")
        for seq, names in c.same_seq_diff_name[:limit]:
            lines.append(f"  • {seq}  →  {', '.join(names)}")
        if len(c.same_seq_diff_name) > limit:
            lines.append(f"  … and {len(c.same_seq_diff_name) - limit} more")
    if c.same_name_diff_seq:
        lines.append("\nSame name, different sequences:")
        for name, seqs in c.same_name_diff_seq[:limit]:
            lines.append(f"  • {name}  →  {', '.join(seqs)}")
        if len(c.same_name_diff_seq) > limit:
            lines.append(f"  … and {len(c.same_name_diff_seq) - limit} more")
    lines.append("\nThe primers were loaded as-is. Check the file if this is unexpected.")
    return "\n".join(lines)


def ask_for_path(parent: QWidget) -> str:
    path, _ = QFileDialog.getOpenFileName(parent, "Select primer file", "",
                                          PRIMER_FILTER)
    return path


def _ask_for_sheet(parent: QWidget, path: str) -> str | int | None:
    """The worksheet to read; None when the user cancelled."""
    if not Path(path).suffix.lower().startswith(".xls"):
        return 0
    try:
        names = pd.ExcelFile(path).sheet_names
    except Exception as e:
        QMessageBox.critical(parent, "Primer file", str(e))
        return None
    if len(names) <= 1:
        return 0
    name, ok = QInputDialog.getItem(parent, "Worksheet", "Sheet:", names, 0, False)
    return name if ok else None


def load_pairs_interactive(
    parent: QWidget, path: str = "",
    log: Callable[[str], None] | None = None,
) -> list[PrimerPair] | None:
    """Pairs from `path`, asking for the file when none is given.

    Returns None when the user cancelled or the file could not be read — in the
    latter case the error has already been shown. When `log` is given, any
    duplicate-primer warning is also written there, so it survives the dialog.
    """
    path = path or ask_for_path(parent)
    if not path:
        return None

    sheet = _ask_for_sheet(parent, path)
    if sheet is None:
        return None

    try:
        if is_fasta(path):
            pairs = load_pairs(path)        # records taken two at a time
        else:
            # Primer sheets are assay-specific, so the mapping is always shown
            # with its guess pre-filled rather than applied silently.
            headers, rows = read_table(path, sheet)
            dlg = ColumnMapDialog(headers, rows, parent)
            if dlg.exec() != QDialog.DialogCode.Accepted:
                return None
            pairs = dlg.pairs
    except PrimerFileError as e:
        QMessageBox.critical(parent, "Primer file", str(e))
        return None
    except Exception as e:
        QMessageBox.critical(parent, "Primer file", f"{e}")
        return None

    if not pairs:
        QMessageBox.warning(parent, "Primer file", "No primer pairs were loaded.")
        return None

    # Same sequence under two names, or one name over two sequences, is almost
    # always a file mistake. Load anyway, but tell the user what to check.
    pset = primer_set_from_pairs(pairs)
    conflicts = find_primer_conflicts([*pset.forwards, *pset.reverses])
    if conflicts.any:
        message = _format_conflicts(conflicts)
        QMessageBox.warning(parent, "Primer file", message)
        if log is not None:
            log(f"WARNING — {Path(path).name}: {message}")
    return pairs
