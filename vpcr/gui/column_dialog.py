"""Manual column mapping for primer tables.

The automatic detection covers the common layouts, but primer sheets are
assay-specific and there is no general format. This dialog shows the detected
mapping, lets the user correct it, and previews the resulting pairs live.
"""
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QFormLayout,
                             QHBoxLayout, QLabel, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

from vpcr.core.primers import PrimerFileError, detect_columns, pairs_from_rows

NONE_LABEL = "— none —"

# key -> label shown in the form. `reverse_extra` is deliberately absent: the
# core still builds an extra pair from a second reverse column (the CLI uses
# it), but detect_columns guesses it from any third DNA-looking column, and a
# wrong guess silently doubled the pairs. Here it is opt-out by omission.
FIELDS = (
    ("name", "Pair name"),
    ("forward_name", "Forward primer name"),
    ("forward", "Forward sequence"),
    ("reverse_name", "Reverse primer name"),
    ("reverse", "Reverse sequence"),
)
REQUIRED = ("forward", "reverse")


class ColumnMapDialog(QDialog):
    def __init__(self, headers: list[str], rows: list[list[str]], parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle("Map primer columns")
        self.resize(900, 620)
        self.headers, self.rows = headers, rows
        self.pairs: list = []

        try:
            proposal = detect_columns(headers, rows)
        except PrimerFileError:
            proposal = {}

        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(
            "Columns detected automatically. Correct them if this assay uses a "
            "different layout; the preview updates as you change them."))

        body = QHBoxLayout()

        form_host = QWidget()
        form = QFormLayout(form_host)
        self.combos: dict[str, QComboBox] = {}
        labels = [f"{i}: {h}" if str(h).strip() else f"{i}: (blank)"
                  for i, h in enumerate(headers)]
        for key, label in FIELDS:
            cb = QComboBox()
            cb.addItem(NONE_LABEL, -1)
            for i, text in enumerate(labels):
                cb.addItem(text, i)
            idx = proposal.get(key, -1)
            cb.setCurrentIndex(idx + 1 if idx >= 0 else 0)
            cb.currentIndexChanged.connect(self._refresh)
            self.combos[key] = cb
            form.addRow(label, cb)
        body.addWidget(form_host, 1)

        right = QVBoxLayout()
        right.addWidget(QLabel("Preview"))
        self.preview = QTableWidget(0, 3)
        self.preview.setHorizontalHeaderLabels(["Pair", "Forward", "Reverse"])
        self.preview.horizontalHeader().setStretchLastSection(True)
        right.addWidget(self.preview, 1)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        right.addWidget(self.status)
        body.addLayout(right, 2)
        lay.addLayout(body)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        lay.addWidget(self.buttons)

        self._refresh()

    def column_map(self) -> dict[str, int]:
        out = {}
        for key, cb in self.combos.items():
            v = cb.currentData()
            if v is not None and v >= 0:
                out[key] = v
        return out

    def _refresh(self) -> None:
        cols = self.column_map()
        ok_btn = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        missing = [k for k in REQUIRED if k not in cols]
        if missing:
            self.pairs = []
            self.preview.setRowCount(0)
            self.status.setText("Select the forward and reverse sequence columns.")
            ok_btn.setEnabled(False)
            return
        try:
            self.pairs = pairs_from_rows(self.headers, self.rows, cols)
        except PrimerFileError as e:
            self.pairs = []
            self.preview.setRowCount(0)
            self.status.setText(f"<span style='color:#b00'>{e}</span>")
            ok_btn.setEnabled(False)
            return

        self.preview.setRowCount(len(self.pairs))
        for r, p in enumerate(self.pairs):
            for c, val in enumerate((p.name, p.forward, p.reverse)):
                item = QTableWidgetItem(val)
                item.setFlags(Qt.ItemFlag.ItemIsEnabled)
                self.preview.setItem(r, c, item)
        self.preview.resizeColumnsToContents()
        self.status.setText(f"{len(self.pairs)} primer pair(s).")
        ok_btn.setEnabled(True)
