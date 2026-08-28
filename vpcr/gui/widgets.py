"""Input widgets that do not change value when the wheel rolls over them.

A QSpinBox or QComboBox inside a scroll area is a trap: scrolling the panel
past one silently edits it, and a run then uses a parameter nobody chose. Qt
wires the wheel to the value by default and offers no flag to unwire it.

The wheel is ignored unconditionally, not merely when the widget lacks focus.
Focus is no evidence of intent here: Qt hands it to the first spin box as soon
as the panel is shown, and clicking a field to type in it leaves it focused, so
a later scroll over that same field would edit it again — precisely the accident
this exists to prevent. Set these with the keyboard or the arrows.

An ignored wheel event propagates to the parent, so the panel scrolls as if the
widget were not there. WheelFocus is dropped too, or a passing wheel could give
the widget focus on its way past.
"""
from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPalette, QWheelEvent
from PyQt6.QtWidgets import (QApplication, QComboBox, QDoubleSpinBox, QSpinBox)


# The accent blue reused for button feedback; matches the primer-map palette so
# the whole app reads as one system. Lighter under a dark theme for contrast.
_ACCENT_LIGHT = QColor("#4C78A8")
_ACCENT_DARK = QColor("#5B8FC9")


def _mix(base: QColor, tint: QColor, t: float) -> QColor:
    """`t` of `tint` blended over `base` (0 = base, 1 = tint)."""
    return QColor(round(base.red() * (1 - t) + tint.red() * t),
                  round(base.green() * (1 - t) + tint.green() * t),
                  round(base.blue() * (1 - t) + tint.blue() * t))


def button_stylesheet(pal: QPalette) -> str:
    """A QPushButton style with a light-blue hover and a stronger-blue press.

    Colours are derived from the palette so the resting button still matches the
    theme; only the hover and pressed states carry the accent. A light-blue wash
    on hover, then the solid accent (white text) while held, so the two states
    read as clearly distinct feedback.
    """
    window = pal.color(QPalette.ColorRole.Window)
    dark = window.lightness() < 128
    accent = _ACCENT_DARK if dark else _ACCENT_LIGHT
    button = pal.color(QPalette.ColorRole.Button)
    text = pal.color(QPalette.ColorRole.ButtonText)
    border = _mix(button, text, 0.35)
    hover = _mix(button, accent, 0.30 if dark else 0.22)   # light blue wash
    # A real grey for disabled text: PlaceholderText equals the text colour under
    # some styles, which would render a disabled button as solid black.
    muted = _mix(text, window, 0.55)
    return f"""
        QPushButton {{
            background: {button.name()};
            color: {text.name()};
            border: 1px solid {border.name()};
            border-radius: 4px;
            padding: 4px 12px;
        }}
        QPushButton:hover {{
            background: {hover.name()};
            border-color: {accent.name()};
        }}
        QPushButton:pressed {{
            background: {accent.name()};
            color: #ffffff;
            border-color: {accent.darker(120).name()};
        }}
        QPushButton:disabled {{
            color: {muted.name()};
            border-color: {_mix(button, window, 0.5).name()};
        }}
    """


def apply_button_style(app: QApplication) -> None:
    """Install the button feedback style app-wide, covering dialogs too."""
    app.setStyleSheet(button_stylesheet(app.palette()))


class _NoWheelMixin:
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        event.ignore()


class NoWheelSpinBox(_NoWheelMixin, QSpinBox):
    pass


class NoWheelDoubleSpinBox(_NoWheelMixin, QDoubleSpinBox):
    pass


class NoWheelComboBox(_NoWheelMixin, QComboBox):
    pass
