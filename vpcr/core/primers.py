"""Loading primer pairs from whatever the assay happens to use.

Primer files differ between assays, so nothing here assumes a fixed layout.
Three shapes are recognised automatically:

  * FASTA        -- records taken two at a time (forward, reverse)
  * table        -- TSV/CSV/XLSX with a name column and F/R sequence columns
  * pair table   -- name<TAB>forward<TAB>reverse  (the JAR's own FRpairs format)

For tables the columns are located by header heuristics, and the caller can
always override them explicitly.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable, NamedTuple, Sequence

from .model import PAIR_NAME_SEP, Primer, PrimerPair

IUPAC = set("ACGTUNMRWSYKVHDBIEFJL")
_UNAMBIGUOUS = set("ACGTU")

# The JAR itself discards anything shorter than 12 nt.
MIN_PRIMER_LEN = 12


def clean_seq(s: str) -> str:
    return re.sub(r"[^A-Za-z]", "", str(s or "")).upper()


def is_dna(s: str) -> bool:
    """Whether a cleaned string could be a primer.

    Length and base composition both matter: a pair name like '12Sa-12Sh'
    cleans to 'SASH', which is all-IUPAC but is obviously not a primer.
    """
    s = s.strip().upper()
    if len(s) < MIN_PRIMER_LEN or not all(c in IUPAC for c in s):
        return False
    plain = sum(1 for c in s if c in _UNAMBIGUOUS)
    return plain / len(s) >= 0.5


class PrimerFileError(ValueError):
    pass


# --------------------------------------------------------------------------
# column detection
# --------------------------------------------------------------------------

_F_HINT = re.compile(r"(forward|fwd|\bf\d*\b|primer\s*f)", re.I)
_R_HINT = re.compile(r"(reverse|\brev\b|\br\d*\b|primer\s*r)", re.I)


def _dna_columns(n_cols: int, rows: Sequence[Sequence[str]]) -> list[int]:
    """Columns whose non-empty cells are mostly primer-like sequences.

    Judged over non-empty cells only: an optional column such as a second
    reverse primer is populated in just one row out of ten.
    """
    out = []
    for i in range(n_cols):
        vals = [clean_seq(r[i]) for r in rows if i < len(r) and str(r[i]).strip()]
        if not vals:
            continue
        if sum(is_dna(v) for v in vals) / len(vals) >= 0.6:
            out.append(i)
    return out


# Columns that hold amplicon lengths look like names but are not.
_SIZE_HINT = re.compile(r"(size|tama|length|longitud|\bbp\b)", re.I)


def _looks_like_names(values: Sequence[str], header: str = "") -> bool:
    if _SIZE_HINT.search(header):
        return False
    vals = [str(v).strip() for v in values if str(v).strip()]
    if not vals:
        return False
    return sum(1 for v in vals if len(v) <= 40) / len(vals) >= 0.8


def detect_columns(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> dict[str, int]:
    """Locate the F/R sequence columns and the name column that labels each.

    Sequence columns are found by content. Their labels are found by position:
    across primer tables in the wild, a primer's name sits in the column
    immediately to its left ('PrimerF' | 'Secuencia Primer F').
    """
    hdr = [str(h or "") for h in headers]
    dna_cols = _dna_columns(len(hdr), rows)
    if len(dna_cols) < 2:
        raise PrimerFileError(
            "could not find two columns of primer sequence; map the columns manually"
        )

    def hint_score(i: int, hint: re.Pattern[str]) -> int:
        return 1 if hint.search(hdr[i]) else 0

    # Prefer header hints; fall back to left-to-right order.
    fwd = min(dna_cols, key=lambda i: (-hint_score(i, _F_HINT), i))
    rest = [i for i in dna_cols if i != fwd]
    rev = min(rest, key=lambda i: (-hint_score(i, _R_HINT), i))
    cols: dict[str, int] = {"forward": fwd, "reverse": rev}

    extra = [i for i in rest if i != rev]
    if extra:
        cols["reverse_extra"] = extra[0]

    def name_col_left_of(c: int) -> int | None:
        # Judged purely on content: a blank header ('Unnamed: 5' in pandas) is
        # exactly what the second-reverse name column looks like.
        j = c - 1
        if j < 0 or j in dna_cols:
            return None
        vals = [r[j] for r in rows if j < len(r)]
        return j if _looks_like_names(vals, hdr[j]) else None

    for key, col in (("forward", fwd), ("reverse", rev),
                     ("reverse_extra", extra[0] if extra else -1)):
        if col >= 0 and (nc := name_col_left_of(col)) is not None:
            cols[f"{key}_name"] = nc

    # A single label column (name | F | R) names the pair, not one primer.
    if "reverse_name" not in cols and "forward_name" in cols:
        cols["name"] = cols.pop("forward_name")
    if "name" not in cols and "forward_name" not in cols:
        for j in range(len(hdr)):
            if j in dna_cols:
                continue
            if _looks_like_names([r[j] for r in rows if j < len(r)], hdr[j]):
                cols["name"] = j
                break
    return cols


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------

def _from_fasta(path: Path) -> list[PrimerPair]:
    names, seqs = [], []
    cur: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
            names.append(line[1:].strip())
        elif line:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    if len(names) != len(seqs) or len(names) < 2:
        raise PrimerFileError("FASTA primer file must contain an even number of records")
    if len(names) % 2:
        raise PrimerFileError(
            f"FASTA has {len(names)} records; pairs are read two at a time"
        )

    pairs = []
    for i in range(0, len(names), 2):
        fn, rn = names[i], names[i + 1]
        pairs.append(PrimerPair(
            name=_common_name(fn, rn),
            forward=clean_seq(seqs[i]), reverse=clean_seq(seqs[i + 1]),
            forward_name=fn, reverse_name=rn,
        ))
    return pairs


def _common_name(fn: str, rn: str) -> str:
    """Join both primer names verbatim: 'MiFish-U-F|MiFish-U-R'.

    Names are never shortened to a shared stem. Collapsing 'MiFish-U-F' and
    'MiFish-U-R' into 'MiFish-U' loses which primers were actually used, and
    two different pairs can share a stem.
    """
    fn, rn = fn.strip().lstrip(">").strip(), rn.strip().lstrip(">").strip()
    if not fn or not rn:
        return fn or rn
    return f"{fn}{PAIR_NAME_SEP}{rn}"


def _rows_from_delimited(path: Path) -> tuple[list[str], list[list[str]]]:
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    sample = text[:4096]
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    rows = [r for r in csv.reader(text.splitlines(), delimiter=delim) if any(c.strip() for c in r)]
    if not rows:
        raise PrimerFileError("empty primer table")
    # A header row is one whose cells are not DNA.
    header = rows[0] if not all(is_dna(clean_seq(c)) for c in rows[0][1:] if c.strip()) else []
    body = rows[1:] if header else rows
    if not header:
        header = [""] * max(len(r) for r in body)
    return header, body


def _rows_from_excel(path: Path, sheet: str | int = 0) -> tuple[list[str], list[list[str]]]:
    try:
        import pandas as pd
    except ImportError as e:  # pragma: no cover
        raise PrimerFileError("reading .xlsx requires pandas + openpyxl") from e
    df = pd.read_excel(path, sheet_name=sheet, dtype=str).fillna("")
    return [str(c) for c in df.columns], df.astype(str).values.tolist()


def _pairs_from_rows(
    header: list[str], rows: list[list[str]], cols: dict[str, int] | None = None
) -> list[PrimerPair]:
    cols = cols or detect_columns(header, rows)
    fi, ri = cols["forward"], cols["reverse"]
    r2i = cols.get("reverse_extra")
    fni, rni, r2ni = (cols.get(k) for k in ("forward_name", "reverse_name",
                                            "reverse_extra_name"))
    ni = cols.get("name")

    pairs: list[PrimerPair] = []
    for n, row in enumerate(rows, 1):
        get = lambda i: (str(row[i]) if i is not None and 0 <= i < len(row) else "")  # noqa: E731
        f, r = clean_seq(get(fi)), clean_seq(get(ri))
        if not (is_dna(f) and is_dna(r)):
            continue

        fname = get(fni).strip().lstrip(">")
        rname = get(rni).strip().lstrip(">")
        if fname and rname:
            name = _common_name(fname, rname)
        elif ni is not None and get(ni).strip():
            name = get(ni).strip().lstrip(">")
        else:
            name = f"pair_{n}"
        pairs.append(PrimerPair(name=name, forward=f, reverse=r,
                                forward_name=fname, reverse_name=rname))

        r2 = clean_seq(get(r2i))
        if is_dna(r2):
            r2name = get(r2ni).strip().lstrip(">")
            name2 = _common_name(fname, r2name) if (fname and r2name) else f"{name}_R2"
            pairs.append(PrimerPair(name=name2, forward=f, reverse=r2,
                                    forward_name=fname, reverse_name=r2name))
    if not pairs:
        raise PrimerFileError("no valid primer pairs found")
    return _dedupe(pairs)


def _dedupe(pairs: list[PrimerPair]) -> list[PrimerPair]:
    """Make pair names unique; identical F/R sequences collapse to one pair."""
    seen_seq: dict[tuple[str, str], PrimerPair] = {}
    used: set[str] = set()
    out: list[PrimerPair] = []
    for p in pairs:
        key = (p.forward, p.reverse)
        if key in seen_seq:
            continue
        name, base, i = p.name, p.name, 2
        while name in used:
            name, i = f"{base}_{i}", i + 1
        used.add(name)
        p2 = PrimerPair(name, p.forward, p.reverse, p.forward_name, p.reverse_name,
                        p.marker, p.expected_size)
        seen_seq[key] = p2
        out.append(p2)
    return out


def is_fasta(path: str | Path) -> bool:
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls", ".xlsm"):
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    return any(line.lstrip().startswith(">") for line in text.splitlines())


def read_table(path: str | Path, sheet: str | int = 0) -> tuple[list[str], list[list[str]]]:
    """Header row and body rows of a primer table, for manual column mapping."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls", ".xlsm"):
        return _rows_from_excel(path, sheet)
    return _rows_from_delimited(path)


def pairs_from_rows(
    header: list[str], rows: list[list[str]], cols: dict[str, int] | None = None
) -> list[PrimerPair]:
    """Build pairs from a table, using an explicit column map when given."""
    return _pairs_from_rows(header, rows, cols)


def load_pairs(
    path: str | Path,
    *,
    sheet: str | int = 0,
    cols: dict[str, int] | None = None,
) -> list[PrimerPair]:
    """Load primer pairs, detecting the file shape."""
    path = Path(path)
    if not path.exists():
        raise PrimerFileError(f"primer file not found: {path}")

    if is_fasta(path):
        return _from_fasta(path)

    header, rows = read_table(path, sheet)
    return _pairs_from_rows(header, rows, cols)


def write_pair_fasta(pair: PrimerPair, path: str | Path) -> None:
    """Write one pair as the two-record FASTA the JAR consumes."""
    Path(path).write_text(
        f">{pair.fwd_name}\n{pair.forward}\n>{pair.rev_name}\n{pair.reverse}\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# individual primers and free recombination
# ---------------------------------------------------------------------------

class PrimerSet(NamedTuple):
    """What a primer file holds: the oligos, and the pairing it suggests."""
    forwards: list[Primer]
    reverses: list[Primer]
    pairs: list[PrimerPair]


def pair_of(f: Primer, r: Primer) -> PrimerPair:
    """Combine any forward with any reverse, whatever file they came from."""
    return PrimerPair(name=_common_name(f.name, r.name),
                      forward=f.sequence, reverse=r.sequence,
                      forward_name=f.name, reverse_name=r.name)


def combine(forwards: Sequence[Primer], reverses: Sequence[Primer]) -> list[PrimerPair]:
    """Every forward against every reverse, forward-major order."""
    return [pair_of(f, r) for f in forwards for r in reverses]


def _dedupe_primers(primers: Iterable[Primer]) -> list[Primer]:
    """Keep first occurrence; a primer reused across pairs is listed once."""
    seen: dict[tuple[str, str], Primer] = {}
    for p in primers:
        seen.setdefault(p.key, p)
    return list(seen.values())


def primer_set_from_pairs(pairs: Sequence[PrimerPair]) -> PrimerSet:
    """Split pairs into their individual primers, preserving file order."""
    fwd = _dedupe_primers(Primer(p.fwd_name, p.forward, "forward") for p in pairs)
    rev = _dedupe_primers(Primer(p.rev_name, p.reverse, "reverse") for p in pairs)
    return PrimerSet(fwd, rev, list(pairs))


class PrimerConflicts(NamedTuple):
    """Suspicious overlaps found among loaded primers, after exact dedup.

    A clean primer file names each oligo once and each name once. Two mistakes
    break that: the same sequence entered twice under different names, or one
    name reused for two different sequences. Neither is an error the loader can
    resolve, so they are reported to the user rather than silently merged.
    """
    same_seq_diff_name: list[tuple[str, list[str]]]   # sequence -> its names
    same_name_diff_seq: list[tuple[str, list[str]]]   # name -> its sequences

    @property
    def any(self) -> bool:
        return bool(self.same_seq_diff_name or self.same_name_diff_seq)


def find_primer_conflicts(primers: Iterable[Primer]) -> PrimerConflicts:
    """Group primers by sequence and by name to surface ambiguous overlaps.

    Exact (name, sequence) repeats are dropped first: a forward reused across
    many pairs is not a conflict. What remains is one sequence under several
    names, or one name over several sequences. Sequences are compared
    case-insensitively so 'acgt' and 'ACGT' count as the same oligo.
    """
    names_of: dict[str, list[str]] = {}   # sequence -> names seen for it
    seqs_of: dict[str, list[str]] = {}    # name -> sequences seen for it
    seen: set[tuple[str, str]] = set()
    for p in primers:
        seq = p.sequence.upper()
        if (p.name, seq) in seen:
            continue
        seen.add((p.name, seq))
        names_of.setdefault(seq, []).append(p.name)
        seqs_of.setdefault(p.name, []).append(seq)
    return PrimerConflicts(
        same_seq_diff_name=[(s, ns) for s, ns in names_of.items() if len(ns) > 1],
        same_name_diff_seq=[(n, ss) for n, ss in seqs_of.items() if len(ss) > 1],
    )


def load_primer_set(path: str | Path, *, sheet: str | int = 0,
                    cols: dict[str, int] | None = None) -> PrimerSet:
    """Load a primer file as individual primers plus its suggested pairing."""
    return primer_set_from_pairs(load_pairs(path, sheet=sheet, cols=cols))
