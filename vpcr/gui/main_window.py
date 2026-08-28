"""PyQt6 front end for in silico PCR.

Layout: targets + primer summary + parameters on the left, run log and result
tables on the right. The JAR runs in a worker thread; the batch itself is
parallel, so the GUI thread only marshals progress.

Choosing primers happens in its own window (gui/primer_dialog.py), so the left
panel holds a summary and a drop zone rather than three cramped tables.
"""
from __future__ import annotations

import sys
import traceback
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Callable, NamedTuple

import pandas as pd
from PyQt6 import sip
from PyQt6.QtCore import (QAbstractTableModel, QModelIndex, QSettings, Qt,
                          QThread, QTimer, pyqtSignal)
from PyQt6.QtGui import QAction, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox,
                             QComboBox, QDialog, QFileDialog, QFormLayout,
                             QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                             QListWidget, QMainWindow, QMessageBox,
                             QPlainTextEdit, QProgressBar, QPushButton,
                             QScrollArea, QSplitter, QTabWidget,
                             QTableView, QVBoxLayout, QWidget)

from vpcr import APP_NAME, APP_ROOT, JAR_DIR
from vpcr.gui.alignment_view import AlignmentTab
from vpcr.gui.collapsible import CollapsibleBox
from vpcr.gui.primer_map import PrimerMapTab
from vpcr.gui.dnd import (FASTA_SUFFIXES, JAR_SUFFIXES, DropDirLineEdit,
                          DropLineEdit, DropListWidget, DropZone)
from vpcr.gui.primer_dialog import PrimerDialog, PrimerState
from vpcr.gui.primer_io import load_pairs_interactive
from vpcr.gui.widgets import (NoWheelComboBox, NoWheelSpinBox,
                             apply_button_style)
from vpcr.core.export import (EXCEL_ROW_LIMIT, detail_table, export_run,
                              run_folder_name, summary_table, unique_dir,
                              write_tables)
from vpcr.core.model import PrimerPair, RunConfig
from vpcr.core.primers import MIN_PRIMER_LEN
from vpcr.core.runner import (JarError, discover_javas, jar_class_version,
                              run_batch, select_java)

# Four rows of 30 px. The list used to take whatever height the panel had left
# over — 192 px to show two files — because it carried the only stretch factor.
TARGET_LIST_HEIGHT = 115

# The single source of truth for the Parameters section: the widgets are built
# empty and filled from here, so "Reset to defaults" cannot drift from what the
# app starts with. Tuned for mitochondrial barcoding rather than the JAR's own
# defaults (30-3000 bp, linear).
PARAM_DEFAULTS: dict[str, object] = {
    "min_len": 90,
    "max_len": 2200,
    "err3": 1,
    "workers": 4,
    "heap": "",
    "circular": True,
    "probe": False,
    "seqext": True,
    "keep": False,
}

_PARAM_LABELS = {
    "min_len": "Min amplicon", "max_len": "Max amplicon", "err3": "3' mismatches",
    "workers": "Parallel jobs", "heap": "JVM heap", "circular": "Circular template",
    "probe": "Probe search", "seqext": "Extract amplicon sequences",
    "keep": "Keep raw .out reports",
}


def _fmt_default(value: object) -> str:
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value) if value != "" else "auto"


def _defaults_tooltip() -> str:
    lines = "\n".join(f"  {_PARAM_LABELS[k]}: {_fmt_default(v)}"
                      for k, v in PARAM_DEFAULTS.items())
    return f"Restore every parameter below to:\n{lines}"


def _col(df: pd.DataFrame, name: str) -> pd.Series:
    """A column as a NaN-safe numeric/boolean series, absent columns as False."""
    if name not in df:
        return pd.Series(False, index=df.index)
    return df[name].fillna(0)


# Detail rows the table shows at once. Past this the view holds only a prefix:
# QTableView bogs down on very large models and no one reads 50k rows by eye. The
# full table always reaches the CSV/Excel export, and the count line says so.
DETAIL_VIEW_CAP = 50_000


class DetailFilter(NamedTuple):
    label: str
    apply: Callable[[pd.DataFrame], pd.DataFrame]


# 'Multi-amplicon' = more than one predicted product. 'Multiband' = products of
# more than one distinct size, i.e. what a gel would actually resolve.
DETAIL_FILTERS: tuple[DetailFilter, ...] = (
    DetailFilter("All", lambda d: d),
    DetailFilter("Amplified", lambda d: d[_col(d, "amplified") == True]),  # noqa: E712
    DetailFilter("Not amplified", lambda d: d[_col(d, "amplified") != True]),  # noqa: E712
    DetailFilter("Multi-amplicon (>1 product)", lambda d: d[_col(d, "n_amplicons") > 1]),
    DetailFilter("Multiband (>1 size)", lambda d: d[_col(d, "n_unique_sizes") > 1]),
)


# ---------------------------------------------------------------------------
# table model
# ---------------------------------------------------------------------------

class DataFrameModel(QAbstractTableModel):
    """Read-only view over a DataFrame; cheap for large result sets."""

    def __init__(self, df: pd.DataFrame | None = None) -> None:
        super().__init__()
        self._df = df if df is not None else pd.DataFrame()

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        self._df = df
        self.endResetModel()

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._df)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._df.columns)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        value = self._df.iat[index.row(), index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            if pd.isna(value):
                return ""
            if isinstance(value, float):
                return f"{value:g}"
            return str(value)
        if role == Qt.ItemDataRole.BackgroundRole:
            col = self._df.columns[index.column()]
            if col == "amplified":
                return QColor(219, 245, 224) if bool(value) else QColor(250, 226, 226)
        return None

    def headerData(self, section: int, orientation: Qt.Orientation,
                   role: int = Qt.ItemDataRole.DisplayRole):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return str(self._df.columns[section])
        return str(section + 1)


# ---------------------------------------------------------------------------
# worker
# ---------------------------------------------------------------------------

class RunWorker(QThread):
    progress = pyqtSignal(int, int, str)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, pairs, targets, jar, config, workers, heap, reports_dir):
        super().__init__()
        self._args = (pairs, targets, jar, config, workers, heap, reports_dir)
        self._cancel = False

    def cancel(self) -> None:
        self._cancel = True

    def run(self) -> None:
        pairs, targets, jar, config, workers, heap, reports_dir = self._args
        try:
            results = run_batch(
                pairs, targets, jar, config,
                workers=workers, heap=heap or None, reports_dir=reports_dir,
                progress=lambda d, t, l: self.progress.emit(d, t, l),
                should_cancel=lambda: self._cancel,
            )
            self.finished_ok.emit(results)
        except Exception as e:
            self.failed.emit(f"{e}\n\n{traceback.format_exc()}")


# ---------------------------------------------------------------------------
# main window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.results: list = []
        self.worker: RunWorker | None = None
        self._detail_df: pd.DataFrame | None = None
        # Folder chosen once per run so the .out reports and the export share it.
        self._run_outdir: Path | None = None
        # The primer tables live in a dialog that is rebuilt on every open, so
        # the window, not a widget, owns what was loaded.
        self.primers = PrimerState()
        self._settings = QSettings("IAvH", APP_NAME)
        # Set before any widget exists: resizeEvent fires during construction.
        self._restoring = True
        self._pending_splitter = None
        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(400)
        self._save_timer.timeout.connect(self._save_ui_state)

        self.splitter = QSplitter(Qt.Orientation.Horizontal)
        self.splitter.addWidget(self._build_left())
        self.splitter.addWidget(self._build_right())
        self.splitter.setStretchFactor(0, 0)
        self.splitter.setStretchFactor(1, 1)
        self.splitter.setSizes([470, 810])
        self.splitter.splitterMoved.connect(lambda *_: self._schedule_save())
        self.setCentralWidget(self.splitter)

        self._build_menu()
        self._autodetect_jar()
        self._wire_summaries()
        self._restore_ui_state()   # clears _restoring on the way out
        self.statusBar().showMessage("Ready")

    # -- parameters --------------------------------------------------------

    def parameters_are_default(self) -> bool:
        return all(self._param_value(k) == v for k, v in PARAM_DEFAULTS.items())

    def _param_value(self, key: str) -> object:
        w = getattr(self, key)
        if isinstance(w, QCheckBox):
            return w.isChecked()
        if isinstance(w, QComboBox):
            return w.currentText()
        return w.value()

    def reset_parameters(self, confirm: bool = False) -> None:
        """Restore every Parameters widget to PARAM_DEFAULTS."""
        if confirm and not self.parameters_are_default():
            answer = QMessageBox.question(
                self, APP_NAME, "Restore all parameters to their default values?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes)
            if answer != QMessageBox.StandardButton.Yes:
                return
        for key, value in PARAM_DEFAULTS.items():
            w = getattr(self, key)
            if isinstance(w, QCheckBox):
                w.setChecked(bool(value))
            elif isinstance(w, QComboBox):
                w.setCurrentText(str(value))
            else:
                w.setValue(int(value))

    def _on_section_toggled(self, _expanded: bool) -> None:
        self._refresh_tail_stretch()
        self._save_ui_state()

    def _refresh_tail_stretch(self) -> None:
        """Give the trailing spacer the spare height only when nothing else wants it."""
        growing = any(b.is_expanded() and b.stretch for b in self.sections.values())
        self._left_layout.setStretch(self._tail_index, 0 if growing else 1)

    # -- section summaries -------------------------------------------------

    def _wire_summaries(self) -> None:
        """Keep each header's right-hand hint current, folded or not."""
        refresh = self._refresh_summaries
        self.jar_edit.textChanged.connect(refresh)
        self.run_label.textChanged.connect(refresh)
        self.out_edit.textChanged.connect(refresh)
        self.auto_export.toggled.connect(refresh)
        for w in (self.min_len, self.max_len, self.err3):
            w.valueChanged.connect(refresh)
        for w in (self.circular, self.probe, self.seqext, self.keep):
            w.toggled.connect(refresh)
        m = self.target_list.model()
        m.rowsInserted.connect(refresh)
        m.rowsRemoved.connect(refresh)
        m.modelReset.connect(refresh)
        refresh()

    def _dim_primer_names(self) -> None:
        """The pair list is a hint, not a control: dimmed, from the palette."""
        self.primer_names.setStyleSheet(
            f"color: {self.palette().placeholderText().color().name()};")

    def changeEvent(self, event) -> None:  # noqa: N802
        from PyQt6.QtCore import QEvent
        # Qt can deliver this before _build_left has run.
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "primer_names"):
            self._dim_primer_names()
            # Rebuild the accent hover/press style for the new theme.
            if (app := QApplication.instance()) is not None:
                apply_button_style(app)
        super().changeEvent(event)

    def _refresh_primer_labels(self) -> None:
        """The panel's stand-in for the tables: how many pairs, and which."""
        st = self.primers
        if st.is_empty:
            self.primer_count.setText("No primers loaded")
            self.primer_names.setText("")
            self.primer_names.setToolTip("")
            return

        chosen = st.selected
        noun = "pair" if len(st.pairs) == 1 else "pairs"
        self.primer_count.setText(
            f"<b>{len(chosen)}</b> of {len(st.pairs)} {noun} to run, from "
            f"{len(st.forwards)} forward and {len(st.reverses)} reverse primers")
        if not chosen:
            self.primer_names.setText("no pairs ticked")
            self.primer_names.setToolTip("")
            return
        shown = ", ".join(p.name for p in chosen[:2])
        if len(chosen) > 2:
            shown += f", +{len(chosen) - 2} more"
        self.primer_names.setText(shown)
        self.primer_names.setToolTip("\n".join(p.name for p in chosen))

    def _refresh_summaries(self, *_) -> None:
        # Reached from target_list's model signals (rowsInserted/Removed/
        # modelReset). On window teardown the model is cleared and emits
        # modelReset after the C++ widgets are gone; the surviving Python
        # wrapper would then touch dead objects and raise. Bail out, as
        # DropListWidget._place_caption does for the same signal.
        if sip.isdeleted(self):
            return
        jar = self.jar_edit.text().strip()
        self.sections["jar"].set_summary(Path(jar).name if jar else "not set")

        n = self.target_list.count()
        self.sections["targets"].set_summary(
            "none" if not n else f"{n} file{'s' if n > 1 else ''}")

        self.sections["primers"].set_summary(self.primers.summary())
        self._refresh_primer_labels()

        bits = [f"{self.min_len.value()}–{self.max_len.value()} bp",
                f"{self.err3.value()} mm",
                "circular" if self.circular.isChecked() else "linear"]
        if self.probe.isChecked():
            bits.append("probe")
        at_default = self.parameters_are_default()
        if not at_default:
            bits.append("modified")
        self.sections["parameters"].set_summary(" • ".join(bits))
        self.reset_btn.setEnabled(not at_default)

        label = self.run_label.text().strip()
        out = [label or Path(self.out_edit.text().strip() or ".").name]
        out.append("auto-export" if self.auto_export.isChecked() else "manual export")
        self.sections["output"].set_summary(" • ".join(out))

    # -- persisted layout --------------------------------------------------

    def _save_ui_state(self) -> None:
        # set_expanded() emits toggled, which lands here: without this guard a
        # restore would overwrite later sections with their default state.
        if self._restoring:
            return
        s = self._settings
        # saveGeometry already carries the maximised/fullscreen flag.
        s.setValue("geometry", self.saveGeometry())
        s.setValue("splitter", self.splitter.saveState())
        for key, box in self.sections.items():
            s.setValue(f"section/{key}", box.is_expanded())
        # Closing the console that launched the app kills the process without a
        # closeEvent, so nothing may flush the settings later. Do it now.
        s.sync()

    def _schedule_save(self) -> None:
        """Coalesce the flood of events a drag produces into one write."""
        if not self._restoring:
            self._save_timer.start()

    def _restore_ui_state(self) -> None:
        s = self._settings
        self._restoring = True
        try:
            if (geom := s.value("geometry")) is not None:
                self.restoreGeometry(geom)
            else:
                self.resize(1280, 760)
            # Applied on first show: a splitter restored before the window has
            # its final width redistributes its panes when that width arrives.
            self._pending_splitter = s.value("splitter")

            for key, box in self.sections.items():
                val = s.value(f"section/{key}")
                if val is not None:
                    # QSettings round-trips booleans as strings on some backends.
                    box.set_expanded(val in (True, "true", "True", 1, "1"))
        finally:
            self._restoring = False
        self._refresh_tail_stretch()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if (state := self._pending_splitter) is not None:
            self._pending_splitter = None
            self._restoring = True
            try:
                self.splitter.restoreState(state)
            finally:
                self._restoring = False

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._schedule_save()

    def moveEvent(self, event) -> None:  # noqa: N802
        super().moveEvent(event)
        self._schedule_save()

    def closeEvent(self, event) -> None:  # noqa: N802
        self._save_timer.stop()
        self._save_ui_state()
        super().closeEvent(event)

    # -- construction ------------------------------------------------------

    def _build_menu(self) -> None:
        m = self.menuBar().addMenu("&File")
        for text, slot in (("Add FASTA targets...", self.add_targets),
                           ("Load primer pairs...", lambda: self.load_primers()),
                           ("Manage primers...", self.manage_primers),
                           ("Export results...", self.export_results)):
            a = QAction(text, self)
            a.triggered.connect(slot)
            m.addAction(a)
        m.addSeparator()
        quit_a = QAction("Quit", self)
        quit_a.triggered.connect(self.close)
        m.addAction(quit_a)

        h = self.menuBar().addMenu("&Help")
        a = QAction("Java / JAR info", self)
        a.triggered.connect(self.show_env)
        h.addAction(a)

    def _build_left(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        self.sections: dict[str, CollapsibleBox] = {}

        def section(key: str, title: str, layout, stretch: int = 0,
                    expanded: bool = True) -> None:
            box = CollapsibleBox(title, expanded)
            box.set_content_layout(layout)
            box.toggled.connect(self._on_section_toggled)
            self.sections[key] = box
            lay.addWidget(box, stretch)
            box.set_stretch(stretch)   # after addWidget: it adjusts the layout

        # JAR: autodetected, so it starts folded away.
        jl = QHBoxLayout()
        self.jar_edit = DropLineEdit(JAR_SUFFIXES)
        self.jar_edit.setMinimumWidth(120)
        self.jar_edit.setPlaceholderText("virtualPCR.jar")
        self.jar_edit.setToolTip("Path to virtualPCR.jar — or drop the file here")
        self.jar_edit.fileDropped.connect(lambda p: self._log(f"JAR: {p}"))
        b = QPushButton("Browse")
        b.clicked.connect(self.pick_jar)
        jl.addWidget(self.jar_edit)
        jl.addWidget(b)
        section("jar", "virtualPCR JAR", jl, expanded=False)

        # targets
        tl = QVBoxLayout()
        self.target_list = DropListWidget(FASTA_SUFFIXES)
        self.target_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.target_list.setToolTip("Drop FASTA files or a folder of them here")
        # A cap, not a floor: QAbstractScrollArea hands out a 256x192 sizeHint,
        # which a smaller minimumHeight does nothing about.
        self.target_list.setMinimumHeight(70)
        self.target_list.setMaximumHeight(TARGET_LIST_HEIGHT)
        self.target_list.filesDropped.connect(
            lambda ps: self._log(f"Added {len(ps)} target file(s)"))
        tl.addWidget(self.target_list)
        row = QHBoxLayout()
        for text, slot in (("Add", self.add_targets), ("Remove", self.remove_targets),
                           ("Clear", lambda: self.target_list.clear())):
            btn = QPushButton(text)
            btn.clicked.connect(slot)
            row.addWidget(btn)
        tl.addLayout(row)
        section("targets", "Target FASTA files", tl)

        # primers: a summary of what is loaded; the tables live in a dialog
        pl = QVBoxLayout()
        self.primer_count = QLabel()
        self.primer_count.setWordWrap(True)
        self.primer_names = QLabel()
        self.primer_names.setObjectName("primerNames")
        self._dim_primer_names()
        pl.addWidget(self.primer_count)
        pl.addWidget(self.primer_names)

        self.primer_drop = DropZone()
        self.primer_drop.fileDropped.connect(self.load_primers)
        self.primer_drop.clicked.connect(lambda: self.load_primers())
        pl.addWidget(self.primer_drop)

        # Slots with optional arguments are wrapped: clicked(bool) would
        # otherwise arrive as the first parameter.
        row = QHBoxLayout()
        for text, tip, slot in (
                ("Load file...", "Load primers from FASTA, TSV, CSV or XLSX",
                 lambda: self.load_primers()),
                ("View/Select primers...", "Choose primers, recombine them, edit the pairs",
                 self.manage_primers),
                ("Clear", "Discard the loaded primers, so another file can be dropped",
                 self.clear_primers)):
            btn = QPushButton(text)
            btn.setToolTip(tip)
            btn.clicked.connect(slot)
            row.addWidget(btn)
        pl.addLayout(row)
        section("primers", "Primers", pl)

        # parameters. The spin boxes and the heap combo ignore the wheel: they
        # sit in a scroll area, and scrolling past one used to edit it.
        form = QFormLayout()
        self.min_len = NoWheelSpinBox(); self.min_len.setRange(20, 50000)
        self.max_len = NoWheelSpinBox(); self.max_len.setRange(20, 50000)
        self.err3 = NoWheelSpinBox(); self.err3.setRange(0, 10)
        self.workers = NoWheelSpinBox(); self.workers.setRange(1, 32)
        self.heap = NoWheelComboBox(); self.heap.addItems(["", "2g", "4g", "8g", "16g"])
        self.circular = QCheckBox("Circular template")
        self.probe = QCheckBox("Probe search")
        self.seqext = QCheckBox("Extract amplicon sequences")
        self.keep = QCheckBox("Keep raw .out reports")
        form.addRow("Min amplicon (bp)", self.min_len)
        form.addRow("Max amplicon (bp)", self.max_len)
        form.addRow("3' mismatches", self.err3)
        form.addRow("Parallel jobs", self.workers)
        form.addRow("JVM heap", self.heap)
        form.addRow(self.circular)
        form.addRow(self.probe)
        form.addRow(self.seqext)
        form.addRow(self.keep)

        reset_row = QHBoxLayout()
        reset_row.addStretch(1)
        self.reset_btn = QPushButton("Reset to defaults")
        self.reset_btn.setToolTip(_defaults_tooltip())
        self.reset_btn.clicked.connect(lambda: self.reset_parameters(confirm=True))
        reset_row.addWidget(self.reset_btn)
        form.addRow(reset_row)

        self.reset_parameters()   # the widgets have no values until now
        section("parameters", "Parameters", form)

        # output
        ol = QVBoxLayout()
        orow = QHBoxLayout()
        self.out_edit = DropDirLineEdit()
        self.out_edit.setText(str(APP_ROOT / "output"))
        self.out_edit.setMinimumWidth(120)
        self.out_edit.setToolTip(
            "Each run is written to its own <primers>_<timestamp> subfolder here.\n"
            "Drop a folder to change it.")
        ob = QPushButton("Browse")
        ob.clicked.connect(self.pick_output_dir)
        orow.addWidget(self.out_edit)
        orow.addWidget(ob)
        ol.addLayout(orow)
        lrow = QHBoxLayout()
        lrow.addWidget(QLabel("Run label"))
        self.run_label = QLineEdit()
        self.run_label.setMinimumWidth(120)
        self.run_label.setPlaceholderText("optional")
        self.run_label.setToolTip("Names the output subfolder instead of the primer pairs")
        lrow.addWidget(self.run_label)
        ol.addLayout(lrow)
        self.auto_export = QCheckBox("Auto-export when the run finishes")
        self.auto_export.setToolTip(
            "Write summary.csv, detail.csv, vpcr_results.xlsx and run_log.txt\n"
            "into a new subfolder as soon as the run completes")
        self.auto_export.setChecked(True)
        ol.addWidget(self.auto_export)
        section("output", "Output", ol)

        # Soaks up the leftover height when no expanded section wants it;
        # without it QVBoxLayout spreads that height between the headers.
        lay.addStretch(0)
        self._left_layout = lay
        self._tail_index = lay.count() - 1

        # The sections scroll: with every one expanded they are taller than a
        # 768 px screen, and folding them should be a choice, not a duty.
        scroll = QScrollArea()
        scroll.setWidget(w)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        # Labels are kept short enough that the panel fits without it, but never
        # forbid the bar: hiding a control is worse than showing a scrollbar.
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        # Run sits outside the scroll area: it is the point of the panel, and
        # inside it the sections could push it past the bottom of the screen.
        host = QWidget()
        host_lay = QVBoxLayout(host)
        host_lay.setContentsMargins(0, 0, 0, 9)
        host_lay.addWidget(scroll, 1)

        run_row = QHBoxLayout()
        run_row.setContentsMargins(9, 0, 9, 0)
        self.run_btn = QPushButton("Run in silico PCR")
        self.run_btn.clicked.connect(self.start_run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.cancel_run)
        run_row.addWidget(self.run_btn, 2)
        run_row.addWidget(self.cancel_btn, 1)
        host_lay.addLayout(run_row)

        self.progress = QProgressBar()
        host_lay.addWidget(self.progress)
        host.setMinimumWidth(360)
        return host

    def _build_right(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        self.tabs = QTabWidget()

        self.summary_model = DataFrameModel()
        self.detail_model = DataFrameModel()
        for model, title in ((self.summary_model, "Summary"),
                             (self.detail_model, "Detail")):
            view = QTableView()
            view.setModel(model)
            view.setSortingEnabled(False)
            view.setAlternatingRowColors(True)
            view.horizontalHeader().setSectionResizeMode(
                QHeaderView.ResizeMode.ResizeToContents)
            self.tabs.addTab(view, title)

        # Primers drawn onto the amplified region, for a chosen pair/sequence.
        self.alignment_tab = AlignmentTab()
        self.tabs.addTab(self.alignment_tab, "Alignment")

        # Every primer's footprint on a chosen target, to see what overlaps.
        self.primer_map_tab = PrimerMapTab()
        self.tabs.addTab(self.primer_map_tab, "Primer map")

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.tabs.addTab(self.log, "Log")
        lay.addWidget(self.tabs, 1)

        row = QHBoxLayout()
        row.addWidget(QLabel("Detail rows:"))
        self.detail_filter = QComboBox()
        self.detail_filter.addItems([f.label for f in DETAIL_FILTERS])
        # Its longest item would otherwise set the window's minimum width.
        self.detail_filter.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.detail_filter.setMinimumContentsLength(14)
        self.detail_filter.currentIndexChanged.connect(self._apply_detail_filter)
        self.detail_filter.setEnabled(False)
        row.addWidget(self.detail_filter)
        self.filter_count = QLabel("")
        row.addWidget(self.filter_count)

        self.export_btn = QPushButton("Export results (Excel + CSV)")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self.export_results)
        row.addStretch(1)
        row.addWidget(self.export_btn)
        lay.addLayout(row)
        return w

    def _apply_detail_filter(self) -> None:
        if self._detail_df is None:
            return
        try:
            with self._busy("Filtering detail rows…"):
                f = DETAIL_FILTERS[self.detail_filter.currentIndex()]
                shown = f.apply(self._detail_df)
                matched = len(shown)
                if matched > DETAIL_VIEW_CAP:
                    self.detail_model.set_dataframe(shown.iloc[:DETAIL_VIEW_CAP])
                    self.filter_count.setText(
                        f"showing first {DETAIL_VIEW_CAP:,} of {matched:,} rows "
                        f"— export for the full table")
                else:
                    self.detail_model.set_dataframe(shown)
                    self.filter_count.setText(f"{matched:,} of {len(self._detail_df):,} rows")
            self.tabs.setCurrentIndex(1)
        except MemoryError:
            self._log("OUT OF MEMORY applying the detail filter; showing nothing.")
            self.detail_model.set_dataframe(pd.DataFrame())
            self.filter_count.setText("out of memory — use a narrower filter or export")

    # -- helpers -----------------------------------------------------------

    def _log(self, msg: str) -> None:
        self.log.appendPlainText(msg)

    @contextmanager
    def _busy(self, message: str):
        """Show a wait cursor + status while a blocking step runs on the GUI thread.

        Building the result tables and writing the export are CPU/IO bound and
        cannot be interrupted; without this the window paints as 'Not responding'
        and Windows offers to kill it. processEvents once lets the cursor and
        message paint before the blocking call starts.
        """
        self.statusBar().showMessage(message)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()
        try:
            yield
        finally:
            QApplication.restoreOverrideCursor()

    def _autodetect_jar(self) -> None:
        """Pick the newest JAR that some installed JRE can actually run."""
        best_feature = max((j.feature for j in discover_javas()), default=0)
        candidates: list[tuple[int, Path]] = []
        for cand in self._jar_candidates():
            try:
                need = jar_class_version(cand)
            except Exception:
                continue
            if need <= best_feature:
                candidates.append((need, cand))
            else:
                self._log(f"Ignoring {cand.name} in {cand.parent.name}/ "
                          f"(needs Java {need}, newest found is {best_feature})")
        if not candidates:
            self._log("No runnable virtualPCR.jar found; select one manually.")
            return
        need, jar = max(candidates)
        self.jar_edit.setText(str(jar))
        self._log(f"Found JAR: {jar} (Java {need}+)")

    @staticmethod
    def _jar_candidates() -> list[Path]:
        """Look in jar/, beside the app, and one level up; never a drive walk."""
        seen: dict[Path, None] = {}
        for base in (JAR_DIR, APP_ROOT, APP_ROOT.parent):
            if not base.is_dir():
                continue
            for pattern in ("virtualPCR.jar", "*/virtualPCR.jar", "*/*/virtualPCR.jar"):
                for p in base.glob(pattern):
                    seen.setdefault(p.resolve(), None)
        return list(seen)

    def pick_jar(self) -> None:
        p, _ = QFileDialog.getOpenFileName(self, "Select virtualPCR.jar", "",
                                           "JAR files (*.jar)")
        if p:
            self.jar_edit.setText(p)

    def add_targets(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA targets", "",
            "FASTA (*.fasta *.fa *.fas *.fna *.txt);;All files (*)")
        self.target_list.add_paths(files)

    def remove_targets(self) -> None:
        for item in self.target_list.selectedItems():
            self.target_list.takeItem(self.target_list.row(item))

    def _set_primers(self, state: PrimerState) -> None:
        self.primers = state
        self._refresh_summaries()

    def load_primers(self, path: str = "") -> None:
        """Load pairs from `path`, or ask for a file when called with none.

        The file's own pairing seeds the run, with every pair ticked; Manage
        primers is where it gets recombined.
        """
        pairs = load_pairs_interactive(self, path, log=self._log)
        if pairs is None:
            return
        state = PrimerState.from_pairs(pairs)
        self._set_primers(state)
        self._log(f"Loaded {len(state.forwards)} forward and {len(state.reverses)} "
                  f"reverse primer(s) from {path or 'file'}, paired as in the file "
                  f"({len(pairs)} pair(s))")

    def clear_primers(self) -> None:
        """Forget the loaded file, so a different one can be dropped."""
        if self.primers.is_empty:
            return
        self._set_primers(PrimerState())
        self._log("Cleared the loaded primers")

    def manage_primers(self) -> None:
        dlg = PrimerDialog(self.primers, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        before = len(self.primers.pairs)
        self._set_primers(dlg.state)
        after = len(self.primers.pairs)
        if after != before:
            self._log(f"Primer pairs: {before} -> {after}")

    def collect_pairs(self) -> list[PrimerPair]:
        return self.primers.selected

    def show_env(self) -> None:
        lines = ["Java runtimes found:"]
        javas = discover_javas()
        lines += [f"  Java {j.feature:<3} {j.executable}" for j in javas] or ["  (none)"]

        jar = self.jar_edit.text().strip()
        if jar and Path(jar).exists():
            try:
                need = jar_class_version(jar)
                lines.append(f"\nJAR: {jar}\nRequires Java {need}+")
                chosen = select_java(jar)
                lines.append(f"Will run with: Java {chosen.feature} ({chosen.executable})")
            except JarError as e:
                lines.append(f"\n{e}")
            except Exception as e:
                lines.append(f"\nJAR: {e}")
        QMessageBox.information(self, "Environment", "\n".join(lines))

    # -- run ---------------------------------------------------------------

    def start_run(self) -> None:
        jar = self.jar_edit.text().strip()
        targets = self.target_list.paths()
        pairs = self.collect_pairs()

        if not jar or not Path(jar).exists():
            QMessageBox.warning(self, APP_NAME, "Select a valid virtualPCR.jar."); return
        if not targets:
            QMessageBox.warning(self, APP_NAME, "Add at least one FASTA target."); return
        if not pairs:
            QMessageBox.warning(self, APP_NAME, "Add at least one primer pair."); return
        if self.min_len.value() > self.max_len.value():
            QMessageBox.warning(self, APP_NAME, "Min amplicon length exceeds max."); return

        # The JAR silently discards primers shorter than this and then fails
        # with 'Failed to open primer's file', which explains nothing.
        short = [f"{p.name}: {p.forward} / {p.reverse}" for p in pairs
                 if min(len(p.forward), len(p.reverse)) < MIN_PRIMER_LEN]
        if short:
            QMessageBox.warning(
                self, APP_NAME,
                f"These primers are shorter than {MIN_PRIMER_LEN} nt and cannot "
                f"be searched:\n\n" + "\n".join(short))
            return

        cfg = RunConfig(
            min_len=self.min_len.value(), max_len=self.max_len.value(),
            n_3prime_errors=self.err3.value(), circular=self.circular.isChecked(),
            probe=self.probe.isChecked(), sequence_extract=self.seqext.isChecked(),
        )
        # One folder per run, named once here so its timestamp is stable between
        # the worker (which writes the .out reports as it goes) and the export at
        # the end. With both options on, 'Keep raw .out' and auto-export land
        # side by side; the .out files go in a 'reports' subfolder of it.
        self._run_outdir = None
        if self.keep.isChecked() or self.auto_export.isChecked():
            self._run_outdir = unique_dir(
                self._output_base(),
                run_folder_name(pairs, datetime.now(), self.run_label.text()))
        reports = (self._run_outdir / "reports") if self.keep.isChecked() else None
        if reports is not None:
            self._log(f"Raw .out reports: {reports}")

        self.progress.setRange(0, len(pairs) * len(targets))
        self.progress.setValue(0)
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.export_btn.setEnabled(False)
        self.alignment_tab.set_results([])
        self.primer_map_tab.set_results([])
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.log))

        # The log is archived next to the results, so it records the run.
        self._log(f"\n{'=' * 62}\nvirtualPCR Studio run — {datetime.now():%Y-%m-%d %H:%M:%S}")
        self._log(f"JAR       : {jar}")
        self._log(f"Amplicon  : {cfg.min_len}-{cfg.max_len} bp, "
                  f"3' mismatches {cfg.n_3prime_errors}")
        self._log(f"Template  : {'circular' if cfg.circular else 'linear'}"
                  f"{', probe search' if cfg.probe else ''}"
                  f"{'' if cfg.sequence_extract else ', no amplicon sequences'}")
        # Cap the echo: thousands of per-item lines stall the log widget at
        # start. The full lists reach summary.csv / the export regardless.
        if len(targets) <= 50:
            for t in targets:
                self._log(f"Target    : {t}")
        else:
            self._log(f"Targets   : {len(targets)} files")
        if len(pairs) <= 50:
            for p in pairs:
                self._log(f"Primer    : {p.name}  {p.forward} / {p.reverse}")
        else:
            self._log(f"Primers   : {len(pairs)} pairs (see summary.csv for the full list)")
        self._log(f"Running {len(pairs)} pair(s) x {len(targets)} file(s)...")

        self.worker = RunWorker(pairs, targets, jar, cfg,
                                self.workers.value(), self.heap.currentText(), reports)
        self.worker.progress.connect(self.on_progress)
        self.worker.finished_ok.connect(self.on_finished)
        self.worker.failed.connect(self.on_failed)
        self.worker.start()

    def cancel_run(self) -> None:
        if self.worker:
            self.worker.cancel()
            self._log("Cancelling after current jobs finish...")

    def on_progress(self, done: int, total: int, label: str) -> None:
        self.progress.setValue(done)
        self._log(f"[{done}/{total}] {label}")
        self.statusBar().showMessage(label)

    def on_finished(self, results: list) -> None:
        self.results = results
        self.cancel_btn.setEnabled(False)
        try:
            if not results:
                self._log("No results.")
                return
            # Building the detail table materialises one row per amplicon; on a
            # huge result set this is where the app would freeze or run out of
            # memory. Do it under a wait cursor and survive a MemoryError.
            with self._busy("Building result tables…"):
                summary = summary_table(results)
                detail = detail_table(results, include_sequence=self.seqext.isChecked())
                self.summary_model.set_dataframe(summary)
                # Sequences are long; hide them on screen but keep them for export.
                self._detail_df = detail.drop(columns=["amplicon_seq"], errors="ignore")
            self.alignment_tab.set_results(results)
            self.primer_map_tab.set_results(results)
            self.detail_filter.setEnabled(True)
            self._apply_detail_filter()
            self.export_btn.setEnabled(True)
            self.tabs.setCurrentIndex(0)
            multiband = int((summary["n_multiband"]).sum())
            multi = int((summary["n_multi_amplicon"]).sum())
            self._log(f"Done. {len(summary)} pair/file combination(s), {len(detail)} detail rows. "
                      f"{multi} sequence(s) with >1 product, {multiband} with >1 distinct size.")
            self.statusBar().showMessage("Done")

            if self.auto_export.isChecked():
                with self._busy("Exporting results…"):
                    outdir = self._export(quiet=True, outdir=self._run_outdir)
                if outdir:
                    self.statusBar().showMessage(f"Exported to {outdir}")
        except MemoryError:
            self._detail_df = None
            self.export_btn.setEnabled(bool(results))
            self._log("OUT OF MEMORY while building the result tables. The run itself "
                      "finished; try fewer pairs/targets per run or turn off "
                      "'Extract amplicon sequences'.")
            QMessageBox.critical(
                self, APP_NAME,
                "Ran out of memory building the result tables.\n\n"
                "The run finished, but the result set is too large to display here. "
                "Split the targets across smaller runs, or turn off "
                "'Extract amplicon sequences', then run again.")
        except Exception as e:
            self._log(f"ERROR building results: {e}\n{traceback.format_exc()}")
            QMessageBox.critical(self, "Results failed", str(e).splitlines()[0] if str(e) else repr(e))
        finally:
            # Whatever happened, the window must be usable again.
            self.run_btn.setEnabled(True)

    def on_failed(self, msg: str) -> None:
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self._log(f"FAILED: {msg}")
        QMessageBox.critical(self, "Run failed", msg.splitlines()[0])

    def pick_output_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Output folder", self.out_edit.text())
        if d:
            self.out_edit.setText(d)

    def _output_base(self) -> Path:
        return Path(self.out_edit.text().strip() or (APP_ROOT / "output"))

    def _export(self, quiet: bool = False, outdir: Path | None = None) -> Path | None:
        """Write the current results into the run's folder, or a fresh one.

        `outdir` is the folder chosen at run start (shared with the .out
        reports). When None — a manual export of a run that kept neither — a
        fresh <label>_<timestamp> folder is created instead.
        """
        base = self._output_base()
        try:
            base.mkdir(parents=True, exist_ok=True)
            if outdir is None:
                outdir, paths = export_run(
                    base, self.results, log_text=self.log.toPlainText(),
                    label=self.run_label.text(),
                    include_sequence=self.seqext.isChecked())
            else:
                paths = write_tables(
                    outdir, self.results, self.log.toPlainText(),
                    include_sequence=self.seqext.isChecked())
        except PermissionError as e:
            # The usual cause is the workbook still open in Excel.
            self._log(f"EXPORT FAILED: {e}")
            QMessageBox.critical(self, "Export failed",
                                 f"{e}\n\nIs the workbook open in Excel?")
            return None
        except Exception as e:
            self._log(f"EXPORT FAILED: {e}")
            QMessageBox.critical(self, "Export failed", str(e))
            return None

        self._log(f"Exported to {outdir}")
        for k, v in paths.items():
            self._log(f"  {k}: {v.name}")
        skipped_xlsx = "xlsx" not in paths
        if skipped_xlsx:
            self._log(f"  xlsx: skipped (> {EXCEL_ROW_LIMIT:,} detail rows; "
                      f"the full table is in detail.csv)")
        # The log file was written before these lines existed; refresh it.
        if "log" in paths:
            paths["log"].write_text(self.log.toPlainText(), encoding="utf-8")
        if not quiet:
            note = ("" if not skipped_xlsx else
                    f"\n\nExcel skipped: over {EXCEL_ROW_LIMIT:,} rows. "
                    f"The full table is in detail.csv.")
            QMessageBox.information(self, APP_NAME, f"Exported to:\n{outdir}{note}")
        return outdir

    def export_results(self) -> None:
        if not self.results:
            QMessageBox.information(self, APP_NAME, "Nothing to export yet."); return
        # Reuse the run's folder when it has one, so a manual export lands with
        # the .out reports / auto-export instead of spawning a parallel folder.
        self._export(outdir=self._run_outdir)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    apply_button_style(app)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
