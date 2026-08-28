"""Alignment tab: the primers drawn onto the amplified region.

virtualPCR already emits, per binding site, the primer strand, a match line and
the template context (ShowPrimerAlignment=true). The parser keeps all three on
every PrimerHit, and the extracted product on every Amplicon, so nothing has to
be recomputed here -- this is presentation over data the run already produced.

For a chosen (pair, target sequence) the tab shows, per predicted amplicon:
  * the product sequence with the forward and reverse primer footprints tinted;
  * the detailed primer/template alignment for each flanking primer, with
    mismatches coloured and 3'-terminal mismatches (the ones that actually cost
    amplification) emphasised.

Selection is explicit: pick a pair, then one of its amplified target sequences.
A run can hold millions of sequences, so the sequence list is capped.
"""
from __future__ import annotations

from html import escape
from pathlib import Path

from PyQt6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QSizePolicy,
                             QTextBrowser, QVBoxLayout, QWidget)

from vpcr.core.model import Amplicon, PairResult, PrimerHit, SequenceResult

# A single pair can amplify more sequences than a combo box should hold; past
# this the list is truncated with a note. The full detail is in the Detail tab
# and the export regardless.
SEQ_COMBO_CAP = 2000

# A mismatch this close to the 3' end (0 = terminal base) is flagged: the
# polymerase will not extend off an unpaired 3' base, so these are what break an
# assay, while a 5' mismatch is usually tolerated.
THREE_PRIME_WINDOW = 5

# The JAR prints the primer as '5-<seq>->' (a two-character lead) with the
# target context flush left, so the aligned region of the context begins two
# bases in. Used as the default when the offset cannot be pinned by content.
DEFAULT_TARGET_OFFSET = 2

# -- IUPAC pairing, for locating the context under the primer ----------------

_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _pairs_with(primer_base: str, target_base: str) -> bool:
    """Whether a primer base can Watson-Crick pair with a target base (IUPAC)."""
    ps = set(_IUPAC.get(primer_base, "ACGT"))
    tc = {_COMP[b] for b in _IUPAC.get(target_base, "ACGT") if b in _COMP}
    return bool(ps & tc)


def _best_target_offset(primer: str, match: str, target: str) -> int:
    """Index in the context where the primer's binding starts.

    The parser strips the report's leading whitespace, dropping the offset, so it
    is recovered by content: the offset that best pairs the primer against the
    context at the match line's '|' columns wins. Ties prefer the structural
    default (2), which is what a well-formed block always uses.
    """
    L = len(primer)
    if not target or L == 0:
        return 0
    strong = [j for j, c in enumerate(match[:L]) if c == "|"]
    best_off, best_score = DEFAULT_TARGET_OFFSET, -1
    for off in range(0, min(len(target), 7)):
        if off + L > len(target):
            break
        score = sum(1 for j in strong if _pairs_with(primer[j], target[off + j]))
        # +0.5 nudge toward the structural default on a tie.
        adj = score + (0.5 if off == DEFAULT_TARGET_OFFSET else 0.0)
        if adj > best_score:
            best_off, best_score = off, adj
    return best_off


# -- HTML colouring ----------------------------------------------------------

_STYLE = {
    "match": "",
    "degen": "color:#c77d00;font-weight:bold;",       # ':' partial/degenerate
    "mism": "color:#d11a1a;font-weight:bold;",         # blank: hard mismatch
    "mism3": "color:#d11a1a;font-weight:bold;background:#ffe08a;",  # 3'-terminal
    "flank": "color:#9a9a9a;",                         # context outside the primer
    "fwd": "background:#cde7d3;color:#10331b;",         # forward footprint
    "rev": "background:#cfe0f5;color:#0e2a45;",         # reverse footprint
    "plain": "",
}


def _spanify(cells: list[tuple[str, str]]) -> str:
    """Render (char, kind) cells to HTML, coalescing runs of one kind."""
    out: list[str] = []
    i = 0
    while i < len(cells):
        kind = cells[i][1]
        j = i
        while j < len(cells) and cells[j][1] == kind:
            j += 1
        text = escape("".join(c for c, _ in cells[i:j]))
        style = _STYLE.get(kind, "")
        out.append(f"<span style='{style}'>{text}</span>" if style else text)
        i = j
    return "".join(out)


def _symbol_kind(symbol: str, dist_from_3prime: int) -> str:
    if symbol == "|":
        return "match"
    if symbol in ":.":
        return "degen"
    return "mism3" if dist_from_3prime < THREE_PRIME_WINDOW else "mism"


# -- per-hit alignment block -------------------------------------------------

def _alignment_block(role: str, hit: PrimerHit) -> str:
    reverse = hit.orientation == "reverse"
    primer = (hit.primer_seq or "").upper()
    # virtualPCR prints the reverse primer 3'->5' (the reverse of the declared
    # 5'->3' sequence), and the match line follows that printed strand. Render
    # the strand as printed so primer, match and context line up column for
    # column; the 5'/3' labels then say which way it reads.
    strand = primer[::-1] if reverse else primer
    match = (hit.match_line or "").rstrip()
    target = (hit.target_context or "").upper()
    n = min(len(strand), len(match))
    off = _best_target_offset(strand[:n], match[:n], target)
    width = max(len(target), off + n)

    primer_row: list[tuple[str, str]] = [(" ", "plain")] * width
    match_row: list[tuple[str, str]] = [(" ", "plain")] * width
    target_row: list[tuple[str, str]] = [(" ", "plain")] * width

    for k, ch in enumerate(target):
        target_row[k] = (ch, "flank")

    offsets_3p: list[int] = []
    for j in range(n):
        col = off + j
        dist3 = j if reverse else n - 1 - j     # the 3' end is leftmost on reverse
        kind = _symbol_kind(match[j], dist3)
        if match[j] != "|":
            offsets_3p.append(dist3)
        primer_row[col] = (strand[j], kind)
        match_row[col] = (match[j], kind)
        if col < len(target):
            target_row[col] = (target[col], kind)

    # Column 0 of each row is the leftmost base printed; label its 5'/3' sense so
    # the reversed reverse-primer strand is not mistaken for the declared one.
    left_end, right_end = ("3'", "5'") if reverse else ("5'", "3'")

    def line(label: str, cells: list[tuple[str, str]]) -> str:
        return f"{label:<10}{_spanify(cells)}"

    body = "\n".join([
        line(left_end, primer_row) + f"  {right_end}   primer",
        line("", match_row),
        line("", target_row) + "   template",
    ])
    n_mm = len(offsets_3p)
    mm_note = "perfect match" if n_mm == 0 else (
        f"{n_mm} mismatch{'es' if n_mm != 1 else ''}"
        + (f", nearest {min(offsets_3p)} nt from 3'" if offsets_3p else ""))
    head = (f"<b>{escape(role)}</b> &nbsp; {escape(hit.primer_name or '?')} &nbsp; "
            f"pos {hit.start:,}&ndash;{hit.end:,} &nbsp; {hit.identity:g}% id &nbsp; "
            f"Tm {hit.tm:g}&deg;C &nbsp; <span style='color:#555'>{mm_note}</span>")
    return (f"<div style='margin:6px 0 2px'>{head}</div>"
            f"<pre style='margin:0 0 10px;font-family:Consolas,\"DejaVu Sans Mono\","
            f"monospace;font-size:12px;line-height:1.35'>{body}</pre>")


# -- product map -------------------------------------------------------------

def _wrap_sequence(cells: list[tuple[str, str]], start: int, width: int = 60) -> str:
    """A sequence of (char, kind) cells, wrapped with a 1-based position gutter."""
    lines = []
    for i in range(0, len(cells), width):
        chunk = cells[i:i + width]
        lines.append(f"<span style='color:#9a9a9a'>{start + i:>7} </span>"
                     f"{_spanify(chunk)}")
    return "\n".join(lines)


def _product_map(amp: Amplicon, fwd: PrimerHit | None, rev: PrimerHit | None) -> str:
    seq = (amp.sequence or "").upper()
    ta = f", Ta={amp.ta:g}&deg;C" if amp.ta is not None else ""
    head = (f"<div style='margin:4px 0'><b>Amplicon</b> &nbsp; {amp.size:,} bp"
            f"{ta} &nbsp; <span style='color:#555'>pos {amp.start:,}&ndash;"
            f"{amp.end:,}</span></div>")
    if not seq:
        return (head + "<p style='color:#777'><i>No product sequence available "
                "(Extract amplicon sequences was off).</i></p>")

    lf = len(fwd.primer_seq) if fwd else 0
    lr = len(rev.primer_seq) if rev else 0
    kinds = ["plain"] * len(seq)
    for k in range(min(lf, len(seq))):
        kinds[k] = "fwd"
    for k in range(max(0, len(seq) - lr), len(seq)):
        if kinds[k] == "plain":          # never let reverse overwrite forward
            kinds[k] = "rev"

    legend = ("<div style='margin:2px 0 4px;font-size:12px'>"
              "<span style='" + _STYLE["fwd"] + "'>&nbsp;forward&nbsp;</span> &nbsp; "
              "<span style='" + _STYLE["rev"] + "'>&nbsp;reverse&nbsp;</span></div>")

    # Elide a long interior: the primer footprints sit at the two ends, so a
    # multi-kb product need not be shown whole to make the point.
    keep = 90
    if len(seq) > 2 * keep + 20:
        head_cells = [(seq[k], kinds[k]) for k in range(keep)]
        tail_cells = [(seq[k], kinds[k]) for k in range(len(seq) - keep, len(seq))]
        gap = len(seq) - 2 * keep
        body = (_wrap_sequence(head_cells, 1)
                + f"\n<span style='color:#9a9a9a'>        … {gap:,} bp …</span>\n"
                + _wrap_sequence(tail_cells, len(seq) - keep + 1))
    else:
        body = _wrap_sequence([(seq[k], kinds[k]) for k in range(len(seq))], 1)

    return (head + legend
            + f"<pre style='margin:0 0 10px;font-family:Consolas,\"DejaVu Sans Mono\","
            f"monospace;font-size:12px;line-height:1.35'>{body}</pre>")


def _flanking(seq: SequenceResult, amp: Amplicon) -> tuple[PrimerHit | None, PrimerHit | None]:
    """The forward hit that opens the amplicon and the reverse hit that closes it."""
    fwd = next((h for h in seq.hits_for("forward") if h.start == amp.start), None)
    rev = next((h for h in seq.hits_for("reverse") if h.end == amp.end), None)
    return fwd, rev


def _render(pair_name: str, res: PairResult, seq: SequenceResult) -> str:
    parts = [f"<div style='font-size:13px;color:#333;margin-bottom:6px'>"
             f"<b>{escape(pair_name)}</b> on <code>{escape(res.pair.forward)}</code> / "
             f"<code>{escape(res.pair.reverse)}</code><br>"
             f"target <code>{escape(seq.header)}</code> "
             f"<span style='color:#777'>({escape(Path(res.target_file).name)})</span></div>"]

    if not seq.amplicons:
        parts.append("<p style='color:#777'>No predicted amplicon for this sequence.</p>")
    for n, amp in enumerate(seq.amplicons, 1):
        fwd, rev = _flanking(seq, amp)
        if len(seq.amplicons) > 1:
            parts.append(f"<h3 style='margin:14px 0 2px'>Product {n} of "
                         f"{len(seq.amplicons)}</h3>")
        parts.append("<hr style='border:none;border-top:1px solid #ddd'>")
        parts.append(_product_map(amp, fwd, rev))
        # When the sequences were not extracted the amplicon carries no
        # coordinates, so it cannot be tied to specific hits; show every hit.
        if fwd is None and rev is None and amp.start == 0:
            for h in seq.hits_for("forward"):
                parts.append(_alignment_block("Forward", h))
            for h in seq.hits_for("reverse"):
                parts.append(_alignment_block("Reverse", h))
        else:
            if fwd is not None:
                parts.append(_alignment_block("Forward", fwd))
            if rev is not None:
                parts.append(_alignment_block("Reverse", rev))
    return "<div style='padding:4px 8px'>" + "".join(parts) + "</div>"


class AlignmentTab(QWidget):
    """Pick a pair and one of its amplified sequences; see the primers on it."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._results: list[PairResult] = []
        # Parallel to seq_combo's items: the (result, sequence) each row means.
        self._seq_entries: list[tuple[PairResult, SequenceResult]] = []

        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Pair:"))
        self.pair_combo = QComboBox()
        self.pair_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.pair_combo.setMinimumContentsLength(16)
        self.pair_combo.currentIndexChanged.connect(self._on_pair_changed)
        row.addWidget(self.pair_combo)
        row.addWidget(QLabel("Target:"))
        self.seq_combo = QComboBox()
        self.seq_combo.setSizePolicy(QSizePolicy.Policy.Expanding,
                                     QSizePolicy.Policy.Fixed)
        self.seq_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.seq_combo.setMinimumContentsLength(24)
        self.seq_combo.currentIndexChanged.connect(self._on_seq_changed)
        row.addWidget(self.seq_combo, 1)
        lay.addLayout(row)

        self.count_label = QLabel("")
        self.count_label.setStyleSheet("color:#777")
        lay.addWidget(self.count_label)

        self.view = QTextBrowser()
        self.view.setOpenExternalLinks(False)
        lay.addWidget(self.view, 1)
        self._show_placeholder("Run in silico PCR, then choose a pair and target here.")

    # -- population --------------------------------------------------------

    def set_results(self, results: list[PairResult]) -> None:
        self._results = results or []
        self.pair_combo.blockSignals(True)
        self.pair_combo.clear()
        names: list[str] = []
        seen: set[str] = set()
        for r in self._results:
            if r.pair.name not in seen and any(s.amplified for s in r.sequences):
                seen.add(r.pair.name)
                names.append(r.pair.name)
        self.pair_combo.addItems(names)
        self.pair_combo.blockSignals(False)

        if not names:
            self.seq_combo.clear()
            self._seq_entries = []
            self.count_label.setText("")
            self._show_placeholder(
                "No amplified sequences to show." if self._results else
                "Run in silico PCR, then choose a pair and target here.")
            return
        self.pair_combo.setCurrentIndex(0)
        self._on_pair_changed()

    def _on_pair_changed(self) -> None:
        name = self.pair_combo.currentText()
        self._seq_entries = []
        for r in self._results:
            if r.pair.name != name:
                continue
            for s in r.sequences:
                if s.amplified:
                    self._seq_entries.append((r, s))

        total = len(self._seq_entries)
        shown = self._seq_entries[:SEQ_COMBO_CAP]
        self.seq_combo.blockSignals(True)
        self.seq_combo.clear()
        for r, s in shown:
            self.seq_combo.addItem(f"{Path(r.target_file).name} ▸ {s.header}")
        self.seq_combo.blockSignals(False)

        if total > SEQ_COMBO_CAP:
            self.count_label.setText(
                f"{total:,} amplified sequences; showing first {SEQ_COMBO_CAP:,}. "
                f"Use the Detail tab or the export for the full set.")
        else:
            self.count_label.setText(f"{total:,} amplified sequence(s)")

        if shown:
            self.seq_combo.setCurrentIndex(0)
            self._on_seq_changed()
        else:
            self._show_placeholder("This pair amplified no sequences.")

    def _on_seq_changed(self) -> None:
        i = self.seq_combo.currentIndex()
        if not (0 <= i < len(self._seq_entries)):
            return
        res, seq = self._seq_entries[i]
        self.view.setHtml(_render(self.pair_combo.currentText(), res, seq))

    # -- helpers -----------------------------------------------------------

    def _show_placeholder(self, msg: str) -> None:
        self.view.setHtml(f"<div style='padding:24px;color:#888'>{escape(msg)}</div>")
