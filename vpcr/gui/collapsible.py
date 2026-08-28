"""A group box that folds away.

Qt has no collapsible container. QGroupBox.setCheckable() looks like one but its
checkbox reads as "enable this section", not "show it", so a clickable header is
used instead: a chevron, a bold title, and a dimmed summary on the right that
tells you what a folded section holds without unfolding it.

Colours come from the palette rather than being hard-coded, so the header keeps
its contrast under a dark theme.
"""
from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QPalette
from PyQt6.QtWidgets import (QFrame, QHBoxLayout, QLabel, QLayout, QSizePolicy,
                             QToolButton, QVBoxLayout, QWidget)


class _Header(QFrame):
    clicked = pyqtSignal()

    def __init__(self, title: str, parent=None) -> None:
        super().__init__(parent)
        self.setObjectName("sectionHeader")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedHeight(26)

        # Drawn by the current style, so it exists whatever the font ships with.
        self.chevron = QToolButton()
        self.chevron.setObjectName("sectionChevron")
        self.chevron.setArrowType(Qt.ArrowType.DownArrow)
        self.chevron.setAutoRaise(True)
        self.chevron.setFixedSize(16, 16)
        self.chevron.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.chevron.clicked.connect(self.clicked)
        self.title = QLabel(title)
        f = self.title.font()
        f.setBold(True)
        self.title.setFont(f)
        self.summary = QLabel("")
        self.summary.setObjectName("sectionSummary")
        self.summary.setAlignment(Qt.AlignmentFlag.AlignRight
                                  | Qt.AlignmentFlag.AlignVCenter)
        # A hint must never widen the panel: let it shrink, and elide instead.
        self.summary.setSizePolicy(QSizePolicy.Policy.Ignored,
                                   QSizePolicy.Policy.Preferred)
        self._summary_text = ""

        lay = QHBoxLayout(self)
        lay.setContentsMargins(6, 0, 8, 0)
        lay.setSpacing(4)
        lay.addWidget(self.chevron)
        lay.addWidget(self.title)
        lay.addStretch(1)
        lay.addWidget(self.summary, 1)

    def set_summary(self, text: str) -> None:
        self._summary_text = text
        self.summary.setToolTip(text)
        self._elide()

    def _elide(self) -> None:
        avail = max(0, self.summary.width())
        fm = self.summary.fontMetrics()
        self.summary.setText(
            fm.elidedText(self._summary_text, Qt.TextElideMode.ElideRight, avail))

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._elide()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)


class CollapsibleBox(QWidget):
    """Titled section whose body can be hidden to reclaim vertical space."""

    toggled = pyqtSignal(bool)

    def __init__(self, title: str, expanded: bool = True, parent=None) -> None:
        super().__init__(parent)
        self._stretch = 0
        self._expanded = expanded

        self.header = _Header(title)
        self.header.clicked.connect(lambda: self.set_expanded(not self._expanded))

        self.body = QFrame()
        self.body.setObjectName("sectionBody")
        self.body.setVisible(expanded)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 6)
        lay.setSpacing(0)
        lay.addWidget(self.header)
        lay.addWidget(self.body)

        self._apply_palette()
        self._sync()

    # -- content -----------------------------------------------------------

    def set_content_layout(self, layout: QLayout) -> None:
        layout.setContentsMargins(8, 6, 4, 6)
        self.body.setLayout(layout)

    def set_summary(self, text: str) -> None:
        """Right-hand hint, e.g. '3 files'. Shown whether folded or not."""
        self.header.set_summary(text)

    # -- folding -----------------------------------------------------------

    def is_expanded(self) -> bool:
        return self._expanded

    def set_expanded(self, expanded: bool) -> None:
        self._expanded = expanded
        self.body.setVisible(expanded)
        self._sync()
        self.toggled.emit(expanded)

    @property
    def stretch(self) -> int:
        """How much spare height this section takes when expanded; 0 = none."""
        return self._stretch

    def set_stretch(self, stretch: int) -> None:
        """Mark the section as one that should grow when expanded."""
        self._stretch = stretch
        self._sync()

    def _sync(self) -> None:
        self.header.chevron.setArrowType(
            Qt.ArrowType.DownArrow if self._expanded else Qt.ArrowType.RightArrow)
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Expanding if self._expanded and self._stretch
            else QSizePolicy.Policy.Maximum)
        # The size policy alone is not enough: a parent QVBoxLayout keeps handing
        # this widget its share of the spare space, leaving a gap where the body
        # used to be. Zero the stretch factor while collapsed.
        parent = self.parentWidget()
        lay = parent.layout() if parent else None
        if isinstance(lay, QVBoxLayout) and (i := lay.indexOf(self)) >= 0:
            lay.setStretch(i, self._stretch if self._expanded else 0)

    # -- appearance --------------------------------------------------------

    def _apply_palette(self) -> None:
        pal = self.palette()
        window = pal.color(QPalette.ColorRole.Window)
        dark = window.lightness() < 128
        header_bg = window.lighter(118) if dark else window.darker(107)
        hover_bg = window.lighter(132) if dark else window.darker(114)
        border = window.lighter(145) if dark else window.darker(122)
        muted = pal.color(QPalette.ColorRole.PlaceholderText)

        self.setStyleSheet(f"""
            QFrame#sectionHeader {{
                background: {header_bg.name()};
                border: 1px solid {border.name()};
                border-radius: 4px;
            }}
            QFrame#sectionHeader:hover {{ background: {hover_bg.name()}; }}
            QToolButton#sectionChevron {{ border: none; background: transparent; }}
            QLabel#sectionSummary {{ color: {muted.name()}; }}
            QFrame#sectionBody {{
                border: 1px solid {border.name()};
                border-top: none;
                border-bottom-left-radius: 4px;
                border-bottom-right-radius: 4px;
                margin: 0px 1px 0px 1px;
            }}
        """)

    def changeEvent(self, event) -> None:  # noqa: N802
        from PyQt6.QtCore import QEvent
        if event.type() == QEvent.Type.PaletteChange:
            self._apply_palette()
        super().changeEvent(event)
