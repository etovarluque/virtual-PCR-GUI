"""Parser for virtualPCR .out reports.

Written against the output shape produced with ShowOnlyAmplicons=true, which is
the only one that reports Ta and the amplicon counters. The parser never matches
on primer names or on the target file name, so it works for any assay.

Report grammar (one block per target sequence that had at least one hit):

    __________________________________________________
    In silico PCR, primers search for: <target path>
    <FASTA header>:

    <primer name>\t<primer seq>

    Position: 15374->15393\t84%\tTm=58.3
    5-gayaaartyccvttycaycc->
      |||||||||||:||||||||
    acctattttaaggcaaggtggggatg

    [more Position: groups for the same primer]
    [more primer sections]

    >15374-15549
    <amplicon sequence, wrapped>

    Total number of amplicons = 1
    Number of amplicons of unique size = 1
    176bp\tTa=50

The trailing per-primer statistics table and the timing line are appended after
the last block:

    PrimerID\tSequence\tHits
    L15411F\tgayaaartyccvttycaycc\t44

    Time taken: 5 seconds
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

from .model import Amplicon, PairResult, PrimerHit, PrimerPair, SequenceResult

BLOCK_SEP = "_" * 50
_SEARCH_LINE = "In silico PCR, primers search for:"
_STATS_HEADER = "PrimerID\tSequence\tHits"

_NUM = r"[\d,]+"


def _int(s: str) -> int:
    return int(s.replace(",", ""))


# Coordinates are plain digits in older builds and thousands-separated in newer
# ones, which also append '(176 bp)' to the amplicon header.
_POSITION = re.compile(
    rf"^Position:\s+({_NUM})\s*(->|<-)\s*({_NUM})\s+([\d.]+)%\s+Tm=(-?[\d.]+)"
)
_AMPLICON_HDR = re.compile(rf"^>({_NUM})-({_NUM})(?:\s*\(\s*\d+\s*bp\s*\))?\s*$")
_TOTAL = re.compile(r"^Total number of amplicons\s*=\s*(\d+)")
_UNIQUE = re.compile(r"^Number of amplicons of unique size\s*=\s*(\d+)")
# Newer builds name the primer pair after Ta: '176bp\tTa=50\tCYTB-F / CYTB-R'
_SIZE_TA = re.compile(rf"^({_NUM})\s*bp(?:\s+Ta=(-?[\d.]+))?")
_PRIMER_DECL = re.compile(r"^(\S.*?)\t([acgtunmrwsykvhdbiefjl]+)\s*$", re.IGNORECASE)
_TIME = re.compile(r"^Time taken:\s*(.+)$")
# Alignment line: only the primer/target pairing symbols, always indented.
_MATCH_LINE = re.compile(r"^\s{2}[|:. ]+$")


def _stream_report(fh: Iterator[str]) -> Iterator[tuple[str, str]]:
    """Yield ('block', text) per report block, then ('trailer', text) last.

    The file is walked line by line so a multi-GB .out never sits in memory
    whole (nor a second time as a split copy): only the current block's lines
    are held, then handed off and dropped. Blocks are delimited by BLOCK_SEP;
    everything from the statistics header on is the trailer, emitted once at the
    end (empty string when the report has none).
    """
    buf: list[str] = []
    trailer: list[str] | None = None
    for line in fh:
        if trailer is not None:
            trailer.append(line)
            continue
        if _STATS_HEADER in line:
            if buf:
                yield "block", "".join(buf)
                buf = []
            trailer = [line]
            continue
        if BLOCK_SEP in line:
            if buf:
                yield "block", "".join(buf)
                buf = []
            continue
        buf.append(line)
    if buf:
        yield "block", "".join(buf)
    yield "trailer", "".join(trailer) if trailer else ""


def _parse_stats(trailer: str) -> tuple[dict[str, int], str]:
    counts: dict[str, int] = {}
    elapsed = ""
    for line in trailer.splitlines():
        if (m := _TIME.match(line.strip())):
            elapsed = m.group(1).strip()
            continue
        parts = line.split("\t")
        if len(parts) == 3 and parts[2].strip().isdigit():
            if parts[0].strip() == "PrimerID":
                continue
            # A pair name may appear twice (once per primer); keep both by
            # keying on name+sequence when they collide.
            key = parts[0].strip()
            if key in counts:
                key = f"{key}\t{parts[1].strip()}"
            counts[key] = int(parts[2])
    return counts, elapsed


def _parse_block(block: str) -> SequenceResult | None:
    lines = block.splitlines()

    # Locate the header: first non-empty line after the "search for:" line.
    header = ""
    start_i = 0
    for i, line in enumerate(lines):
        if line.startswith(_SEARCH_LINE):
            for j in range(i + 1, len(lines)):
                if lines[j].strip():
                    header = lines[j].strip().rstrip(":")
                    start_i = j + 1
                    break
            break
    if not header:
        return None

    res = SequenceResult(header=header)
    cur_name = cur_seq = ""
    pending: PrimerHit | None = None
    amp_start = amp_end = 0
    amp_seq: list[str] = []
    sizes: list[tuple[int, float | None]] = []

    def flush_hit() -> None:
        nonlocal pending
        if pending is not None:
            res.hits.append(pending)
            pending = None

    def flush_amplicon() -> None:
        nonlocal amp_start, amp_end, amp_seq
        if amp_start:
            seq = "".join(amp_seq)
            res.amplicons.append(
                Amplicon(
                    start=amp_start,
                    end=amp_end,
                    size=amp_end - amp_start + 1,
                    sequence=seq,
                )
            )
        amp_start = amp_end = 0
        amp_seq = []

    i = start_i
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if (m := _POSITION.match(stripped)):
            flush_hit()
            a, arrow, b, pct, tm = m.groups()
            a, b = _int(a), _int(b)
            pending = PrimerHit(
                primer_name=cur_name,
                primer_seq=cur_seq,
                start=min(a, b),
                end=max(a, b),
                orientation="forward" if arrow == "->" else "reverse",
                identity=float(pct),
                tm=float(tm),
            )
            # Next three lines: primer strand, match line, target context.
            if i + 2 < len(lines) and _MATCH_LINE.match(lines[i + 2]):
                pending.match_line = lines[i + 2][2:]
                if i + 3 < len(lines):
                    pending.target_context = lines[i + 3].strip()
                i += 3
            i += 1
            continue

        if (m := _AMPLICON_HDR.match(stripped)):
            flush_hit()
            flush_amplicon()
            amp_start, amp_end = _int(m.group(1)), _int(m.group(2))
            i += 1
            continue

        if amp_start and re.fullmatch(r"[acgtnACGTN]+", stripped or "-"):
            amp_seq.append(stripped)
            i += 1
            continue

        if (m := _TOTAL.match(stripped)):
            flush_hit()
            flush_amplicon()
            res.total_amplicons = int(m.group(1))
            i += 1
            continue

        if (m := _UNIQUE.match(stripped)):
            res.unique_sizes = int(m.group(1))
            i += 1
            continue

        if (m := _SIZE_TA.match(stripped)):
            ta = float(m.group(2)) if m.group(2) else None
            sizes.append((_int(m.group(1)), ta))
            i += 1
            continue

        if (m := _PRIMER_DECL.match(line)):
            flush_hit()
            cur_name, cur_seq = m.group(1).strip(), m.group(2).strip()
            i += 1
            continue

        i += 1

    flush_hit()
    flush_amplicon()

    # Attach the reported Ta. The JAR lists one line per unique size, so keying
    # on the span is unambiguous even when two amplicons share a size.
    by_size = {s: ta for s, ta in sizes}
    for amp in res.amplicons:
        amp.ta = by_size.get(amp.size)
    if not res.amplicons and sizes:
        # SequenceExtract=false: only sizes were reported.
        res.amplicons = [Amplicon(start=0, end=0, size=s, ta=ta) for s, ta in sizes]
    if not res.total_amplicons:
        res.total_amplicons = len(res.amplicons)
    if not res.unique_sizes:
        res.unique_sizes = len({a.size for a in res.amplicons})
    return res


def parse_report(
    path: str | Path,
    pair: PrimerPair,
    target_file: str = "",
    all_headers: list[str] | None = None,
) -> PairResult:
    """Parse one virtualPCR report into a PairResult.

    The report is read block by block rather than loaded whole: a run over
    millions of sequences yields a multi-GB .out, and holding it (plus the split
    copy the old parser made) in memory was the parser's dominant cost.

    `all_headers` lets the caller reconcile against the source FASTA: sequences
    with no primer hit at all produce no block, and would otherwise vanish.
    """
    sequences: list[SequenceResult] = []
    trailer = ""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for kind, chunk in _stream_report(fh):
            if kind == "trailer":
                trailer = chunk
            elif _SEARCH_LINE in chunk and (r := _parse_block(chunk)):
                sequences.append(r)
    counts, elapsed = _parse_stats(trailer)

    silent: list[str] = []
    if all_headers:
        seen = {s.header for s in sequences}
        silent = [h for h in all_headers if h not in seen]

    return PairResult(
        pair=pair,
        target_file=target_file or str(path),
        sequences=sequences,
        primer_hit_counts=counts,
        elapsed=elapsed,
        silent_headers=silent,
    )


def fasta_headers(path: str | Path) -> list[str]:
    """Headers of a FASTA file, verbatim and without the leading '>'."""
    out: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].strip())
    return out


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs, streaming one record at a time."""
    header, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def repair_sequences(result: PairResult, fasta_path: str | Path) -> int:
    """Restore the final base the JAR drops from every extracted amplicon.

    virtualPCR reports the amplicon span as end - start + 1 but writes only
    end - start bases. Returns the number of amplicons repaired.
    """
    needed = {
        s.header: s
        for s in result.sequences
        if any(a.sequence and not a.sequence_complete for a in s.amplicons)
    }
    if not needed:
        return 0

    fixed = 0
    for header, seq in iter_fasta(fasta_path):
        target = needed.pop(header, None)
        if target is None:
            continue
        for amp in target.amplicons:
            if not amp.sequence or amp.sequence_complete:
                continue
            missing = amp.size - len(amp.sequence)
            if missing != 1 or amp.end > len(seq):
                continue
            amp.sequence += seq[amp.end - 1].lower()
            fixed += 1
        if not needed:
            break
    return fixed
