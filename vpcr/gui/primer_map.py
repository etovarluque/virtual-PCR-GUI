"""Primer map tab: where every primer binds on a chosen target, and what overlaps.

Each JAR run is one pair, so a given forward primer's binding site is reported in
every pair that uses it. Aggregating the hits of all pairs over one target
sequence therefore recovers the full set of forward and reverse footprints, from
which overlap is read directly off the coordinates -- nothing is recomputed.

The view draws a position axis with forward primers above it and reverse below,
each in a lane; primers whose intervals overlap are forced onto different lanes,
so a stack of bars over one region is exactly a group of primers competing for
the same site. A summary clusters the footprints into sites and flags overlaps.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QPolygonF
from PyQt6.QtWidgets import (QComboBox, QGraphicsScene, QGraphicsSimpleTextItem,
                             QGraphicsView, QHBoxLayout, QLabel, QPushButton,
                             QVBoxLayout, QWidget)

from vpcr.core.model import PairResult, PrimerHit

# A run can hold far more target sequences than a combo should carry; the most
# informative (most distinct primers) are kept and the rest noted.
MAP_SEQ_CAP = 4000

# Geometry, in device pixels (vertical never scales; only the position axis
# zooms, horizontally, via the view's transform).
_LANE_H = 15        # bar height
_LANE_STEP = 20     # vertical distance between lanes
_AXIS_GAP = 12      # axis to first lane
_ARROW = 7          # arrowhead length (scene x; stretches slightly under zoom)

# Brand-neutral categorical palette; assigned per primer name, stable per view.
_PALETTE = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2",
    "#9D755D", "#EECA3B", "#FF9DA6", "#6BA3D6", "#8CB369", "#D68A58",
]


@dataclass
class _Site:
    start: int
    end: int
    names: set[str]


def _sites(hits: list[PrimerHit]) -> list[_Site]:
    """Cluster hits into maximal overlapping runs (connected by coordinate)."""
    sites: list[_Site] = []
    for h in sorted(hits, key=lambda x: (x.start, x.end)):
        if sites and h.start <= sites[-1].end:
            sites[-1].end = max(sites[-1].end, h.end)
            sites[-1].names.add(h.primer_name)
        else:
            sites.append(_Site(h.start, h.end, {h.primer_name}))
    return sites


def _lanes(hits: list[PrimerHit]) -> list[tuple[PrimerHit, int]]:
    """Assign each hit the first lane whose last occupant ends before it starts."""
    lane_end: list[int] = []
    placed: list[tuple[PrimerHit, int]] = []
    for h in sorted(hits, key=lambda x: (x.start, x.end)):
        lane = next((i for i, e in enumerate(lane_end) if e < h.start), None)
        if lane is None:
            lane = len(lane_end)
            lane_end.append(h.end)
        else:
            lane_end[lane] = h.end
        placed.append((h, lane))
    return placed


def _nice_step(span: int, target_ticks: int = 7) -> int:
    """A 1/2/5 x 10^k step that puts roughly `target_ticks` ticks over `span`."""
    if span <= 0:
        return 1
    raw = span / target_ticks
    mag = 10 ** int(math.floor(math.log10(raw))) if raw >= 1 else 1
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return int(m * mag)
    return int(10 * mag)


@dataclass
class _Group:
    """One target sequence and the deduplicated primer hits that landed on it."""
    target_file: str
    header: str
    forward: list[PrimerHit]
    reverse: list[PrimerHit]

    @property
    def n_primers(self) -> int:
        return len({h.primer_name for h in self.forward} |
                   {h.primer_name for h in self.reverse})


class _TrackView(QGraphicsView):
    """A graphics view that zooms horizontally under the cursor with the wheel."""

    def wheelEvent(self, event) -> None:  # noqa: N802
        factor = 1.25 if event.angleDelta().y() > 0 else 1 / 1.25
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.scale(factor, 1.0)   # x only: the lanes must keep their height


class PrimerMapTab(QWidget):
    """Pick a target; see every primer's footprint on it and what overlaps."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._groups: list[_Group] = []
        self._colors: dict[str, QColor] = {}

        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Target:"))
        self.target_combo = QComboBox()
        self.target_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.target_combo.setMinimumContentsLength(30)
        self.target_combo.currentIndexChanged.connect(self._render_current)
        row.addWidget(self.target_combo, 1)
        self.fit_btn = QPushButton("Fit")
        self.fit_btn.setToolTip("Reset the horizontal zoom to fit the whole region")
        self.fit_btn.clicked.connect(self._render_current)
        row.addWidget(self.fit_btn)
        lay.addLayout(row)

        self.info = QLabel("")
        self.info.setWordWrap(True)
        self.info.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(self.info)

        self.scene = QGraphicsScene(self)
        self.view = _TrackView(self.scene)
        self.view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.view.setBackgroundBrush(QBrush(QColor("#fbfbfb")))
        lay.addWidget(self.view, 1)
        self._placeholder("Run in silico PCR, then choose a target to map its primers.")

    # -- population --------------------------------------------------------

    def set_results(self, results: list[PairResult]) -> None:
        # Deduplicate hits per (file, header) across all pairs; a forward reused
        # by several pairs binds at the same place, so it is one footprint.
        groups: dict[tuple[str, str], dict[str, dict]] = {}
        capped = False
        for res in results or []:
            for seq in res.sequences:
                if not seq.hits:
                    continue
                key = (res.target_file, seq.header)
                if key not in groups:
                    if len(groups) >= MAP_SEQ_CAP:
                        capped = True
                        continue
                    groups[key] = {"forward": {}, "reverse": {}}
                bucket = groups[key]
                for h in seq.hits:
                    bucket[h.orientation][(h.primer_name, h.start, h.end)] = h

        self._groups = [
            _Group(f, hdr, list(b["forward"].values()), list(b["reverse"].values()))
            for (f, hdr), b in groups.items()
        ]
        # Most distinct primers first: that sequence shows the overlap structure best.
        self._groups.sort(key=lambda g: g.n_primers, reverse=True)
        self._assign_colors()

        self.target_combo.blockSignals(True)
        self.target_combo.clear()
        for g in self._groups:
            nf = len({h.primer_name for h in g.forward})
            nr = len({h.primer_name for h in g.reverse})
            self.target_combo.addItem(
                f"{Path(g.target_file).name} ▸ {g.header}  ({nf}F / {nr}R)")
        self.target_combo.blockSignals(False)

        if not self._groups:
            self._placeholder("No primer hits to map." if results else
                              "Run in silico PCR, then choose a target to map its primers.")
            return
        note = (f"  ·  showing {MAP_SEQ_CAP:,} of the sequences with hits"
                if capped else "")
        self.target_combo.setToolTip(f"{len(self._groups):,} target sequence(s){note}")
        self.target_combo.setCurrentIndex(0)
        self._render_current()

    def _assign_colors(self) -> None:
        names = sorted({h.primer_name for g in self._groups
                        for h in (g.forward + g.reverse)})
        self._colors = {n: QColor(_PALETTE[i % len(_PALETTE)])
                        for i, n in enumerate(names)}

    # -- rendering ---------------------------------------------------------

    def _render_current(self) -> None:
        i = self.target_combo.currentIndex()
        if not (0 <= i < len(self._groups)):
            return
        g = self._groups[i]
        self.view.resetTransform()
        self.scene.clear()

        hits = g.forward + g.reverse
        if not hits:
            self._placeholder("This sequence has no primer hits.")
            return
        lo = min(h.start for h in hits)
        hi = max(h.end for h in hits)
        span = max(hi - lo, 1)
        vw = max(self.view.viewport().width() - 60, 500)
        sx = min(40.0, max(0.25, vw / span))   # px per bp before user zoom

        def x(pos: int) -> float:
            return (pos - lo) * sx

        self._draw_axis(lo, hi, x)
        self._draw_lanes(g.forward, x, up=True)
        self._draw_lanes(g.reverse, x, up=False)

        self.scene.setSceneRect(self.scene.itemsBoundingRect().adjusted(-30, -30, 30, 30))
        self._update_info(g)

    def _draw_axis(self, lo: int, hi: int, x) -> None:
        pen = QPen(QColor("#666"))
        self.scene.addLine(x(lo), 0, x(hi), 0, pen)
        step = _nice_step(hi - lo)
        tick = (lo // step) * step
        if tick < lo:
            tick += step
        tick_pen = QPen(QColor("#999"))
        while tick <= hi:
            self.scene.addLine(x(tick), -3, x(tick), 3, tick_pen)
            self._text(f"{tick:,}", x(tick), 6, "#666", 8, anchor_center=True)
            tick += step

    def _draw_lanes(self, hits: list[PrimerHit], x, up: bool) -> None:
        for h, lane in _lanes(hits):
            base = _AXIS_GAP + lane * _LANE_STEP
            y = -(base + _LANE_H) if up else base
            self._draw_primer(h, x(h.start), x(h.end + 1), y, up)

    def _draw_primer(self, h: PrimerHit, x0: float, x1: float, y: float,
                     forward: bool) -> None:
        color = self._colors.get(h.primer_name, QColor("#888"))
        head = min(_ARROW, max(1.0, (x1 - x0) * 0.5))
        yb = y + _LANE_H
        ymid = y + _LANE_H / 2
        if forward:   # arrow points right (5'->3')
            pts = [(x0, y), (x1 - head, y), (x1, ymid), (x1 - head, yb), (x0, yb)]
        else:         # reverse: arrow points left
            pts = [(x1, y), (x0 + head, y), (x0, ymid), (x0 + head, yb), (x1, yb)]
        poly = QPolygonF([QPointF(px, py) for px, py in pts])
        item = self.scene.addPolygon(
            poly, QPen(color.darker(140)), QBrush(color.lighter(115)))
        item.setToolTip(
            f"{h.primer_name}  ({'forward' if forward else 'reverse'})\n"
            f"pos {h.start:,}–{h.end:,}  ·  {h.identity:g}% id  ·  Tm {h.tm:g}°C\n"
            f"{h.mismatches} mismatch(es)")
        label = self._text(h.primer_name, x0 + 1, y - 11 if forward else yb + 1,
                           color.darker(150).name(), 8)
        label.setToolTip(item.toolTip())

    def _text(self, s: str, sx: float, sy: float, color: str, pt: int,
              anchor_center: bool = False) -> QGraphicsSimpleTextItem:
        t = self.scene.addSimpleText(s, QFont("Segoe UI", pt))
        t.setBrush(QBrush(QColor(color)))
        # Text keeps its pixel size whatever the horizontal zoom does.
        t.setFlag(QGraphicsSimpleTextItem.GraphicsItemFlag.ItemIgnoresTransformations)
        w = t.boundingRect().width()
        t.setPos(sx - (w / 2 if anchor_center else 0), sy)
        return t

    def _update_info(self, g: _Group) -> None:
        def phrase(hits: list[PrimerHit], role: str) -> str:
            if not hits:
                return f"<b>{role}:</b> none"
            n = len({h.primer_name for h in hits})
            sites = _sites(hits)
            overl = [s for s in sites if len(s.names) > 1]
            base = f"<b>{role}:</b> {n} primer(s) over {len(sites)} site(s)"
            if overl:
                where = "; ".join(f"{s.start:,}–{s.end:,}: " + ", ".join(sorted(s.names))
                                  for s in overl)
                return base + f" — <span style='color:#c0392b'>overlap</span> at {where}"
            return base + " — no overlap"

        self.info.setText(phrase(g.forward, "Forward") + "<br>"
                          + phrase(g.reverse, "Reverse"))

    def _placeholder(self, msg: str) -> None:
        self.scene.clear()
        self.info.setText("")
        t = self.scene.addSimpleText(msg, QFont("Segoe UI", 10))
        t.setBrush(QBrush(QColor("#888")))
        self.scene.setSceneRect(QRectF(-10, -10, 400, 40))
