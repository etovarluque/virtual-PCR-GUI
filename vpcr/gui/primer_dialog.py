"""The Manage primers window.

Selecting primers is a task in its own right: three tables, a cartesian
product, and sequences 20-30 nt wide. In the left panel it competed for height
with the targets and the parameters and lost, showing two rows of six. Here it
has a resizable window to itself, and the panel keeps only a summary.

The dialog edits a copy: Cancel leaves the run configuration untouched.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from PyQt6.QtCore import QSettings, Qt
from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QHBoxLayout,
                             QHeaderView, QLabel, QMessageBox, QPushButton,
                             QSplitter, QTableWidget, QTableWidgetItem,
                             QVBoxLayout, QWidget)

from vpcr import APP_NAME
from vpcr.core.model import Primer, PrimerPair
from vpcr.core.primers import combine, primer_set_from_pairs
from vpcr.gui.primer_io import load_pairs_interactive

# Every pair is a separate JAR run over every target, so a large cartesian
# product is worth confirming before it starts.
COMBINE_WARN_AT = 12


@dataclass
class PrimerState:
    """The primers loaded, and which of their pairs the run should use.

    Held by the main window rather than by a widget: the tables that display it
    now live in a dialog that is destroyed each time it closes.
    """
    forwards: list[Primer] = field(default_factory=list)
    reverses: list[Primer] = field(default_factory=list)
    pairs: list[PrimerPair] = field(default_factory=list)
    used: list[bool] = field(default_factory=list)

    @classmethod
    def from_pairs(cls, pairs: list[PrimerPair]) -> "PrimerState":
        """A freshly loaded file: every pair it suggests is ticked."""
        pset = primer_set_from_pairs(pairs)
        return cls(pset.forwards, pset.reverses, list(pairs), [True] * len(pairs))

    @property
    def selected(self) -> list[PrimerPair]:
        return [p for p, use in zip(self.pairs, self.used) if use]

    @property
    def is_empty(self) -> bool:
        return not (self.forwards or self.reverses or self.pairs)

    def summary(self) -> str:
        """The one-liner the folded section header shows."""
        f, r = len(self.forwards), len(self.reverses)
        if not self.pairs:
            return "none" if not (f or r) else f"{f}F / {r}R, no pairs"
        return f"{f}F / {r}R • {len(self.selected)} of {len(self.pairs)} pairs"


def _check_item(checked: bool) -> QTableWidgetItem:
    it = QTableWidgetItem()
    it.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
    it.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
    return it


def _is_checked(table: QTableWidget, row: int) -> bool:
    it = table.item(row, 0)
    return bool(it) and it.checkState() == Qt.CheckState.Checked


class PrimerDialog(QDialog):
    """Forward and reverse oligos above, the pairs they build below."""

    def __init__(self, state: PrimerState, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("View/Select primers")
        self._settings = QSettings("IAvH", APP_NAME)
        self.state = state
        # Each entry is (table, button, tag) for a select/deselect-all toggle.
        self._toggles: list[tuple[QTableWidget, QPushButton, str]] = []

        lay = QVBoxLayout(self)

        # No drop zone here: the panel behind this window already has one, and
        # this dialog is opened once the file is loaded, to work on its primers.
        self.splitter = QSplitter(Qt.Orientation.Vertical)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.addWidget(self._build_primers())
        self.splitter.addWidget(self._build_pairs())
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 1)
        lay.addWidget(self.splitter, 1)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("OK")
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        lay.addWidget(self.buttons)

        self._load_state(state)
        self._restore_geometry()

    # -- construction ------------------------------------------------------

    def _primer_table(self, role: str) -> QTableWidget:
        t = QTableWidget(0, 3)
        t.setHorizontalHeaderLabels(["Use", "Name", "Sequence"])
        t.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        t.verticalHeader().setVisible(False)
        t.setToolTip(f"Tick the {role} primers, then press Combine selected")
        return t

    def _build_primers(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)

        self.fwd_table = self._primer_table("forward")
        self.rev_table = self._primer_table("reverse")
        self.fwd_toggle = QPushButton()
        self.rev_toggle = QPushButton()
        cols = QHBoxLayout()
        for title, table, toggle, tag in (
                ("Forward", self.fwd_table, self.fwd_toggle, "FORWARD"),
                ("Reverse", self.rev_table, self.rev_toggle, "REVERSE")):
            col = QVBoxLayout()
            col.addWidget(QLabel(f"<b>{title}</b>"))
            col.addWidget(table)
            toggle.setToolTip(f"Select or deselect every {title.lower()} primer")
            self._register_toggle(table, toggle, tag)
            col.addWidget(toggle)
            cols.addLayout(col)
        lay.addLayout(cols, 1)

        row = QHBoxLayout()
        for text, tip, slot in (
                ("Load file...", "Load primers from FASTA, TSV, CSV or XLSX",
                 lambda: self.load_file()),
                ("Combine selected",
                 "Replace the pairs with every ticked forward x ticked reverse",
                 self.combine_selected),
                ("Add selected pair",
                 "Append the one ticked forward + ticked reverse to the pairs list",
                 self.add_selected_pair)):
            btn = QPushButton(text)
            btn.setToolTip(tip)
            btn.clicked.connect(slot)
            row.addWidget(btn)
        lay.addLayout(row)
        return w

    def _build_pairs(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(QLabel("<b>Pairs to run</b>"))

        self.pair_table = QTableWidget(0, 4)
        self.pair_table.setHorizontalHeaderLabels(["Use", "Name", "Forward", "Reverse"])
        hh = self.pair_table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        # It numbers rows the user never refers to by number, and costs width.
        self.pair_table.verticalHeader().setVisible(False)
        lay.addWidget(self.pair_table, 1)

        self.pair_toggle = QPushButton()
        self.pair_toggle.setToolTip("Select or deselect every pair")
        self._register_toggle(self.pair_table, self.pair_toggle, "")

        row = QHBoxLayout()
        row.addWidget(self.pair_toggle)
        for text, slot in (("Remove", self.remove_pair_rows),
                           ("Clear pairs", self.clear_pairs)):
            btn = QPushButton(text)
            btn.clicked.connect(slot)
            row.addWidget(btn)
        lay.addLayout(row)
        return w

    # -- state <-> widgets -------------------------------------------------

    @staticmethod
    def _fill_primer_table(table: QTableWidget, primers: list[Primer]) -> None:
        # itemChanged drives the toggle labels; muted here so a bulk fill does
        # not fire it once per cell (O(rows^2) with the label recompute).
        table.blockSignals(True)
        try:
            table.setRowCount(0)
            for p in primers:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, _check_item(False))
                for c, val in enumerate((p.name, p.sequence), start=1):
                    item = QTableWidgetItem(val)
                    item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    table.setItem(r, c, item)
        finally:
            table.blockSignals(False)
        table.resizeColumnToContents(0)
        table.resizeColumnToContents(1)

    def _load_state(self, state: PrimerState) -> None:
        self._fill_primer_table(self.fwd_table, state.forwards)
        self._fill_primer_table(self.rev_table, state.reverses)
        self.pair_table.blockSignals(True)
        try:
            self.pair_table.setRowCount(0)
            for pair, used in zip(state.pairs, state.used):
                self.add_pair_row(pair.name, pair.forward, pair.reverse, used)
        finally:
            self.pair_table.blockSignals(False)
        # Pair names join both primers ('MiFish-U-F|MiFish-U-R') and are the
        # point of the table; sized once so the column stays draggable.
        self.pair_table.resizeColumnToContents(1)
        self._refresh_all_toggles()

    def current_state(self) -> PrimerState:
        """Read the tables back. Rows missing a sequence are dropped."""
        pairs: list[PrimerPair] = []
        used: list[bool] = []
        for r in range(self.pair_table.rowCount()):
            name = self._cell(r, 1)
            fwd, rev = self._cell(r, 2).upper(), self._cell(r, 3).upper()
            if not (fwd and rev):
                continue
            pairs.append(PrimerPair(name=name or f"pair_{r + 1}",
                                    forward=fwd, reverse=rev))
            used.append(_is_checked(self.pair_table, r))
        return PrimerState(list(self.state.forwards), list(self.state.reverses),
                           pairs, used)

    def _cell(self, row: int, col: int) -> str:
        item = self.pair_table.item(row, col)
        return item.text().strip() if item else ""

    # -- actions -----------------------------------------------------------

    def load_file(self, path: str = "") -> None:
        # Route the duplicate warning into the main window's session log too.
        pairs = load_pairs_interactive(self, path, log=getattr(self.parent(), "_log", None))
        if pairs is None:
            return
        self.state = PrimerState.from_pairs(pairs)
        self._load_state(self.state)

    def _checked_primers(self, table: QTableWidget, role: str) -> list[Primer]:
        return [Primer(table.item(r, 1).text(), table.item(r, 2).text(), role)
                for r in range(table.rowCount()) if _is_checked(table, r)]

    def _register_toggle(self, table: QTableWidget, btn: QPushButton, tag: str) -> None:
        """Wire a single button to check/uncheck every row of `table`.

        The one button does both jobs: its label tracks the table so it reads
        'Select all' until everything is ticked, then 'Deselect all'.
        """
        self._toggles.append((table, btn, tag))
        btn.clicked.connect(lambda _=False: self._toggle_all(table, btn, tag))
        table.itemChanged.connect(lambda *_: self._refresh_toggle(table, btn, tag))
        self._refresh_toggle(table, btn, tag)

    def _all_checked(self, table: QTableWidget) -> bool:
        n = table.rowCount()
        return n > 0 and all(_is_checked(table, r) for r in range(n))

    def _refresh_toggle(self, table: QTableWidget, btn: QPushButton, tag: str) -> None:
        verb = "Deselect" if self._all_checked(table) else "Select"
        btn.setText(f"{verb} all {tag}".rstrip())

    def _refresh_all_toggles(self) -> None:
        for table, btn, tag in self._toggles:
            self._refresh_toggle(table, btn, tag)

    def _toggle_all(self, table: QTableWidget, btn: QPushButton, tag: str) -> None:
        """Check every row, or uncheck them all when they are already checked."""
        st = (Qt.CheckState.Unchecked if self._all_checked(table)
              else Qt.CheckState.Checked)
        table.blockSignals(True)
        try:
            for r in range(table.rowCount()):
                if (it := table.item(r, 0)) is not None:
                    it.setCheckState(st)
        finally:
            table.blockSignals(False)
        self._refresh_toggle(table, btn, tag)

    def _gather_combination(self) -> list[PrimerPair] | None:
        """Ticked forwards x ticked reverses, validated and size-warned.

        Returns the pairs, or None when nothing is ticked or the user backs out
        of a large product. Shared by Combine (replace) and Add (append).
        """
        fwd = self._checked_primers(self.fwd_table, "forward")
        rev = self._checked_primers(self.rev_table, "reverse")
        if not fwd or not rev:
            QMessageBox.information(
                self, APP_NAME, "Tick at least one forward and one reverse primer.")
            return None
        n = len(fwd) * len(rev)
        if n > COMBINE_WARN_AT:
            ok = QMessageBox.question(
                self, APP_NAME,
                f"{len(fwd)} forward x {len(rev)} reverse = {n} pairs.\n\n"
                f"Each pair is a separate JAR run over every target sequence. "
                f"Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if ok != QMessageBox.StandardButton.Yes:
                return None
        return list(combine(fwd, rev))

    def combine_selected(self) -> None:
        """Replace the pairs with the ticked forwards x ticked reverses."""
        pairs = self._gather_combination()
        if pairs is None:
            return
        self.pair_table.blockSignals(True)
        try:
            self.pair_table.setRowCount(0)
            for p in pairs:
                self.add_pair_row(p.name, p.forward, p.reverse)
        finally:
            self.pair_table.blockSignals(False)
        self.pair_table.resizeColumnToContents(1)
        self._refresh_all_toggles()

    def add_selected_pair(self) -> None:
        """Append the single ticked forward + reverse to the pairs list.

        Deliberately one pair only: with several ticked on either side this
        would just repeat Combine selected, so it asks for exactly one each.
        """
        fwd = self._checked_primers(self.fwd_table, "forward")
        rev = self._checked_primers(self.rev_table, "reverse")
        if len(fwd) != 1 or len(rev) != 1:
            QMessageBox.information(
                self, APP_NAME,
                "Tick exactly one forward and one reverse primer to add a single "
                "pair.\nUse Combine selected to build many at once.")
            return
        pair = combine(fwd, rev)[0]
        key = (pair.forward.upper(), pair.reverse.upper())
        existing = {(self._cell(r, 2).upper(), self._cell(r, 3).upper())
                    for r in range(self.pair_table.rowCount())}
        if key in existing:
            QMessageBox.information(self, APP_NAME, "That pair is already in the list.")
            return
        self.add_pair_row(pair.name, pair.forward, pair.reverse)
        self.pair_table.resizeColumnToContents(1)
        self._refresh_all_toggles()

    def add_pair_row(self, name: str = "", fwd: str = "", rev: str = "",
                     used: bool = True) -> None:
        r = self.pair_table.rowCount()
        self.pair_table.insertRow(r)
        self.pair_table.setItem(r, 0, _check_item(used))
        for c, val in enumerate((name, fwd, rev), start=1):
            self.pair_table.setItem(r, c, QTableWidgetItem(val))

    def remove_pair_rows(self) -> None:
        rows = sorted({i.row() for i in self.pair_table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.pair_table.removeRow(r)
        self._refresh_all_toggles()

    def clear_pairs(self) -> None:
        self.pair_table.setRowCount(0)
        self._refresh_all_toggles()

    # -- geometry ----------------------------------------------------------

    def _restore_geometry(self) -> None:
        s = self._settings
        if (geom := s.value("primer_dialog/geometry")) is not None:
            self.restoreGeometry(geom)
        else:
            self.resize(980, 640)
        if (split := s.value("primer_dialog/splitter")) is not None:
            self.splitter.restoreState(split)
        else:
            self.splitter.setSizes([330, 260])

    def done(self, result: int) -> None:  # noqa: N802
        # Saved however the dialog closes, so a Cancel still remembers the size.
        s = self._settings
        s.setValue("primer_dialog/geometry", self.saveGeometry())
        s.setValue("primer_dialog/splitter", self.splitter.saveState())
        s.sync()
        if result == QDialog.DialogCode.Accepted:
            self.state = self.current_state()
        super().done(result)
