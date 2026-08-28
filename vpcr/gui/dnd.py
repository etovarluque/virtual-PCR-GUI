"""Drag-and-drop file widgets.

Each widget accepts only the suffixes it can use, so dropping a FASTA on the
JAR field is refused by the cursor rather than by an error dialog later.
Dropping a folder yields the matching files inside it (non-recursively).

A drag carrying nothing usable is still *accepted as an event* (with
IgnoreAction, so the cursor keeps saying no). Qt stops delivering dragMove and
dragLeave to a widget that ignored the dragEnter, and without those a widget
cannot paint itself as refusing, nor undo that paint when the drag departs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, NamedTuple, Sequence

from PyQt6 import sip
from PyQt6.QtCore import QEvent, Qt, pyqtSignal
from PyQt6.QtGui import (QColor, QDragEnterEvent, QDragLeaveEvent,
                         QDragMoveEvent, QDropEvent, QPalette)
from PyQt6.QtWidgets import (QFrame, QLabel, QLineEdit, QListWidget,
                             QListWidgetItem, QTableWidget, QVBoxLayout,
                             QWidget)

JAR_SUFFIXES = (".jar",)
FASTA_SUFFIXES = (".fasta", ".fa", ".fas", ".fna", ".ffn", ".frn", ".seq", ".txt")
PRIMER_SUFFIXES = (".fasta", ".fa", ".txt", ".tsv", ".csv", ".xlsx", ".xls", ".xlsm")

# Dropping a folder must not sweep up primer files or notes that happen to sit
# in it, so '.txt' and '.seq' are accepted only when dropped explicitly.
FASTA_DIR_SUFFIXES = (".fasta", ".fa", ".fas", ".fna", ".ffn", ".frn")

# The caption's second line. DropZone styles it in CSS and the target list
# paints it; both read these, so the two zones cannot drift in size.
HINT_PX = 11
HINT_GAP = 2   # px between the title and the hint


def _matches(path: Path, suffixes: Sequence[str]) -> bool:
    return path.suffix.lower() in suffixes


def _blend(base: QColor, tint: QColor, t: float) -> QColor:
    """`t` of `tint` over `base`. Tinting beats lightening: a lightened accent
    stays saturated, and the label sitting on it loses its contrast."""
    return QColor(round(base.red() * (1 - t) + tint.red() * t),
                  round(base.green() * (1 - t) + tint.green() * t),
                  round(base.blue() * (1 - t) + tint.blue() * t))


def local_files(event: QDropEvent | QDragEnterEvent) -> list[Path]:
    """Every local path in the drag, whatever its suffix.

    Tells 'this drag carries nothing we could ever use' (no local files at all)
    apart from 'this drag carries the wrong kind of file' — only the second is
    worth showing the user as a refusal.
    """
    md = event.mimeData()
    if not md.hasUrls():
        return []
    return [Path(u.toLocalFile()) for u in md.urls() if u.isLocalFile()]


def extract_paths(event: QDropEvent | QDragEnterEvent,
                  suffixes: Sequence[str],
                  dir_suffixes: Sequence[str] | None = None) -> list[Path]:
    """Local files in the drop that carry an accepted suffix.

    Directories contribute the matching files they contain, which is how a
    folder of per-chromosome FASTA files gets added in one gesture.
    """
    md = event.mimeData()
    if not md.hasUrls():
        return []
    inside = dir_suffixes if dir_suffixes is not None else suffixes
    out: list[Path] = []
    for url in md.urls():
        if not url.isLocalFile():
            continue
        p = Path(url.toLocalFile())
        if p.is_dir():
            out += sorted(c for c in p.iterdir() if c.is_file() and _matches(c, inside))
        elif p.is_file() and _matches(p, suffixes):
            out.append(p)
    return out


class ZoneColors(NamedTuple):
    border: QColor
    background: QColor
    ink: QColor
    muted: QColor


def zone_colors(pal: QPalette, state: str) -> ZoneColors:
    """The dashed border, fill and text a drop target wears in `state`.

    One table for every drop target, so the primer zone and the target list are
    the same object with different words in it. Taken from the palette so a dark
    theme keeps its contrast; the reject red is the colour no palette supplies.

    Idle is the base colour a list already has, with the caption dimmed: a drop
    target at rest should read as an empty field, not as a filled panel.
    """
    window = pal.color(QPalette.ColorRole.Window)
    dark = window.lightness() < 128
    # A real grey, halfway from the text colour to the background. Not
    # PlaceholderText: some styles leave it equal to the text colour, and the
    # caption then renders solid black instead of the soft grey it wants.
    muted = _blend(pal.color(QPalette.ColorRole.WindowText), window, 0.5)

    if state == "accept":
        border = pal.color(QPalette.ColorRole.Highlight)
        return ZoneColors(border, _blend(window, border, 0.16),
                          border.lighter(135) if dark else border.darker(135), muted)
    if state == "reject":
        border = QColor("#e06c6c") if dark else QColor("#c0392b")
        return ZoneColors(border, _blend(window, border, 0.14),
                          border.lighter(120) if dark else border, muted)
    return ZoneColors(window.lighter(150) if dark else window.darker(126),
                      pal.color(QPalette.ColorRole.Base), muted, muted)


class ZoneCaption(QWidget):
    """The two-line 'Drop … here / drag & drop …' label a drop target shows.

    Both zones use this same widget, so both render their text through the one
    path. That matters on Windows: a QLabel draws with ClearType (subpixel)
    while QPainter.drawText falls back to grayscale antialiasing, so a
    hand-painted caption never matches a QLabel one on screen — a difference
    .grab() hides, because an offscreen pixmap has no subpixels either.
    """

    def __init__(self, title: str, hint: str, parent=None) -> None:
        super().__init__(parent)
        self._title = QLabel(title)
        self._title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._hint = QLabel(hint)
        self._hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        # Transparent to the mouse: on the list, clicks must reach the viewport.
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(HINT_GAP)
        lay.addWidget(self._title)
        lay.addWidget(self._hint)

    def set_title(self, text: str) -> None:
        self._title.setText(text)

    def restyle(self, c: ZoneColors) -> None:
        self._title.setStyleSheet(f"color: {c.ink.name()};")
        self._hint.setStyleSheet(f"color: {c.muted.name()}; font-size: {HINT_PX}px;")


class _ZoneStyleMixin:
    """Repaints a drop target from the palette as the drag state changes.

    setStyleSheet re-polishes the widget, which arrives back here as a
    PaletteChange. Without both guards that is an unbounded recursion: the flag
    stops the immediate re-entry, and comparing the sheet stops the round trip a
    parent's own setStyleSheet sets off.
    """

    _state = "none"
    _applying = False

    def drag_state(self) -> str:
        return self._state

    def set_drag_state(self, state: str) -> None:
        if state == self._state:
            return
        self._state = state
        self._on_drag_state(state)
        self._apply_palette()

    def _on_drag_state(self, state: str) -> None:
        """Hook for whatever else the state changes, e.g. a caption."""

    def _apply_palette(self) -> None:
        if self._applying:
            return
        self._applying = True
        try:
            css = self._css(zone_colors(self.palette(), self._state))
            if css != self.styleSheet():
                self.setStyleSheet(css)
        finally:
            self._applying = False

    def _css(self, c: ZoneColors) -> str:
        raise NotImplementedError

    def changeEvent(self, event) -> None:  # noqa: N802
        if event.type() == QEvent.Type.PaletteChange:
            self._apply_palette()
        super().changeEvent(event)


class _DropMixin:
    """Accepts a drag only when it carries at least one usable file."""

    _suffixes: Sequence[str] = ()
    _dir_suffixes: Sequence[str] | None = None

    def _paths(self, event) -> list[Path]:
        return extract_paths(event, self._suffixes, self._dir_suffixes)

    def _drag_state(self, event) -> str:
        """'accept' | 'reject' | 'none' — 'none' means: not our business."""
        if self._paths(event):
            return "accept"
        return "reject" if local_files(event) else "none"

    def set_drag_state(self, state: str) -> None:
        """Paint the widget as accepting or refusing. Overridden by DropZone."""

    def _accept(self, event: QDragEnterEvent | QDragMoveEvent) -> None:
        state = self._drag_state(event)
        self.set_drag_state(state)
        if state == "accept":
            event.setDropAction(event.proposedAction())
            event.acceptProposedAction()
        elif state == "reject":
            # Accept the *event* to keep receiving dragMove/dragLeave, but
            # refuse the *action*, so the cursor still says no.
            event.setDropAction(Qt.DropAction.IgnoreAction)
            event.accept()
        else:
            event.ignore()

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        self._accept(event)

    def dragMoveEvent(self, event: QDragMoveEvent) -> None:  # noqa: N802
        self._accept(event)

    def dragLeaveEvent(self, event: QDragLeaveEvent) -> None:  # noqa: N802
        self.set_drag_state("none")
        super().dragLeaveEvent(event)


class DropLineEdit(_DropMixin, QLineEdit):
    """Single-file field: the last accepted path wins."""

    fileDropped = pyqtSignal(str)

    def __init__(self, suffixes: Sequence[str] = JAR_SUFFIXES, parent=None) -> None:
        super().__init__(parent)
        self._suffixes = suffixes
        self.setAcceptDrops(True)

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        self.set_drag_state("none")
        paths = self._paths(event)
        if not paths:
            event.ignore()
            return
        event.acceptProposedAction()
        self.setText(str(paths[-1]))
        self.fileDropped.emit(str(paths[-1]))


class DropDirLineEdit(QLineEdit):
    """Folder field: accepts a folder, or the folder holding a dropped file."""

    dirDropped = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setAcceptDrops(True)

    @staticmethod
    def _dir_from(event) -> Path | None:
        md = event.mimeData()
        if not md.hasUrls():
            return None
        for url in md.urls():
            if not url.isLocalFile():
                continue
            p = Path(url.toLocalFile())
            if p.is_dir():
                return p
            if p.is_file():
                return p.parent
        return None

    def _accept(self, event) -> None:
        if self._dir_from(event) is not None:
            event.setDropAction(event.proposedAction())
            event.acceptProposedAction()
        else:
            event.ignore()

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        self._accept(event)

    def dragMoveEvent(self, event: QDragMoveEvent) -> None:  # noqa: N802
        self._accept(event)

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        d = self._dir_from(event)
        if d is None:
            event.ignore()
            return
        event.acceptProposedAction()
        self.setText(str(d))
        self.dirDropped.emit(str(d))


class DropListWidget(_ZoneStyleMixin, _DropMixin, QListWidget):
    """Multi-file list that is also a drop zone: drops append, skipping dupes.

    Empty, it says so in the middle of itself rather than looking like a widget
    that failed to load. It answers a drag exactly as DropZone does.
    """

    filesDropped = pyqtSignal(list)

    def __init__(self, suffixes: Sequence[str] = FASTA_SUFFIXES,
                 dir_suffixes: Sequence[str] | None = FASTA_DIR_SUFFIXES,
                 title: str = "Drop FASTA files or a folder here",
                 hint: str = "drag & drop  ·  .fasta .fa .fas .fna", parent=None) -> None:
        super().__init__(parent)
        self._suffixes = suffixes
        self._dir_suffixes = dir_suffixes
        self._title = title
        self.setAcceptDrops(True)
        self.setDragDropMode(QListWidget.DragDropMode.DropOnly)
        self.setObjectName("dropList")
        # A real QLabel caption, centred over the viewport and shown only while
        # the list is empty: the same widget DropZone uses, so the same text
        # rendering. A hand-painted caption would differ under ClearType.
        self._caption = ZoneCaption(title, hint, self.viewport())
        # Show or hide the caption whenever the row count changes, however it
        # changed — add_paths, takeItem, clear(), or a drop.
        self.model().rowsInserted.connect(lambda *_: self._place_caption())
        self.model().rowsRemoved.connect(lambda *_: self._place_caption())
        self.model().modelReset.connect(lambda *_: self._place_caption())
        self._apply_palette()

    def _css(self, c: ZoneColors) -> str:
        self._caption.restyle(c)
        return f"""
            QListWidget#dropList {{
                border: 1px dashed {c.border.name()};
                border-radius: 5px;
                background: {c.background.name()};
            }}
        """

    def _on_drag_state(self, state: str) -> None:
        self._caption.set_title({"accept": "Release to add",
                                 "reject": "Not a FASTA file"}.get(state, self._title))

    def _place_caption(self) -> None:
        # Closing the window clears the model, which emits modelReset after the
        # C++ widget is gone; the connected lambda still fires on the surviving
        # Python wrapper. Touching the dead object then raises. Bail out.
        if sip.isdeleted(self):
            return
        vp = self.viewport().rect()
        self._caption.setVisible(self.count() == 0)
        h = self._caption.sizeHint().height()
        self._caption.setGeometry(0, vp.center().y() - h // 2, vp.width(), h)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._place_caption()

    def paths(self) -> list[str]:
        """Full paths, in row order. Rows display the file name only."""
        return [self.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self.count())]

    def add_paths(self, paths: Iterable[Path | str]) -> list[str]:
        """Append paths not already listed. Long paths would make the sidebar
        scroll sideways, so the row shows the name and the tooltip the path."""
        existing = set(self.paths())
        added = []
        for p in paths:
            s = str(Path(p))
            if s in existing:
                continue
            item = QListWidgetItem(Path(s).name)
            item.setData(Qt.ItemDataRole.UserRole, s)
            item.setToolTip(s)
            self.addItem(item)
            existing.add(s)
            added.append(s)
        return added

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        self.set_drag_state("none")
        paths = self._paths(event)
        if not paths:
            event.ignore()
            return
        event.acceptProposedAction()
        added = self.add_paths(paths)
        if added:
            self.filesDropped.emit(added)


class DropTableWidget(_DropMixin, QTableWidget):
    """Primer table: dropping a primer file reloads the pairs from it."""

    fileDropped = pyqtSignal(str)

    def __init__(self, rows: int, cols: int,
                 suffixes: Sequence[str] = PRIMER_SUFFIXES, parent=None) -> None:
        super().__init__(rows, cols, parent)
        self._suffixes = suffixes
        self.setAcceptDrops(True)
        self.setDragDropMode(QTableWidget.DragDropMode.DropOnly)

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        self.set_drag_state("none")
        paths = self._paths(event)
        if not paths:
            event.ignore()
            return
        event.acceptProposedAction()
        self.fileDropped.emit(str(paths[-1]))


class DropZone(_ZoneStyleMixin, _DropMixin, QFrame):
    """A dashed panel that says what it takes, and shows what it thinks of a drag.

    Idle it is quiet and dimmed. A drag carrying a usable file turns it to the
    palette's highlight colour; one carrying the wrong kind of file turns it
    red. Clicking it opens the file dialog, since a target that reacts to the
    pointer should do something.
    """

    fileDropped = pyqtSignal(str)
    clicked = pyqtSignal()

    def __init__(self, suffixes: Sequence[str] = PRIMER_SUFFIXES,
                 title: str = "Drop a primer file here",
                 hint: str = "drag & drop  ·  FASTA, TSV, CSV or XLSX",
                 dir_suffixes: Sequence[str] | None = None, parent=None) -> None:
        super().__init__(parent)
        self._suffixes = suffixes
        self._dir_suffixes = dir_suffixes
        self._title = title
        self.setAcceptDrops(True)
        self.setObjectName("dropZone")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setMinimumHeight(56)

        # The same caption widget the list uses, so both render text identically.
        self._caption = ZoneCaption(title, hint, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self._caption)
        self._apply_palette()

    # -- appearance --------------------------------------------------------

    def _on_drag_state(self, state: str) -> None:
        self._caption.set_title({"accept": "Release to load",
                                 "reject": "Not a primer file"}.get(state, self._title))

    def _css(self, c: ZoneColors) -> str:
        self._caption.restyle(c)
        return f"""
            QFrame#dropZone {{
                background: {c.background.name()};
                border: 1px dashed {c.border.name()};
                border-radius: 5px;
            }}
        """

    # -- interaction -------------------------------------------------------

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        self.set_drag_state("none")
        paths = self._paths(event)
        if not paths:
            event.ignore()
            return
        event.acceptProposedAction()
        self.fileDropped.emit(str(paths[-1]))
