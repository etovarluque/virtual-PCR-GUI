"""Data model for in silico PCR runs.

Nothing here knows about a particular assay: primer names, FASTA header layout
and biological grouping are all treated as opaque strings.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

Orientation = Literal["forward", "reverse"]

# A pair is named by joining both primer names: 'L15411F|H15546R'.
PAIR_NAME_SEP = "|"


@dataclass(frozen=True)
class Primer:
    """A single primer, kept apart from whatever pairing it arrived with.

    Primer files ship pairs, but the pairing is the assay author's choice, not
    a property of the oligo: the same forward is routinely tried against
    several reverses. Loading primers individually lets the user recombine them.
    """
    name: str
    sequence: str
    role: Orientation

    @property
    def key(self) -> tuple[str, str]:
        return (self.name, self.sequence)


@dataclass(frozen=True)
class PrimerPair:
    """A forward/reverse primer pair. `name` identifies the pair, not the primer."""
    name: str
    forward: str
    reverse: str
    forward_name: str = ""
    reverse_name: str = ""
    marker: str = ""          # optional free-text label (12S, COI, ...)
    expected_size: str = ""   # optional, as reported by the assay author

    def __post_init__(self) -> None:
        if not self.forward or not self.reverse:
            raise ValueError(f"pair {self.name!r}: forward and reverse are required")

    def _split_name(self) -> tuple[str, str] | None:
        """The two primer names encoded in the pair name, when it holds both.

        Callers that build a pair from a name/forward/reverse table (the GUI's
        primer grid, a headerless TSV) have no separate name columns, so the
        individual names have to be recovered from the pair name itself.
        """
        parts = [p.strip() for p in self.name.split(PAIR_NAME_SEP)]
        return (parts[0], parts[1]) if len(parts) == 2 and all(parts) else None

    @property
    def fwd_name(self) -> str:
        if self.forward_name:
            return self.forward_name
        pair = self._split_name()
        return pair[0] if pair else f"{self.name}-F"

    @property
    def rev_name(self) -> str:
        if self.reverse_name:
            return self.reverse_name
        pair = self._split_name()
        return pair[1] if pair else f"{self.name}-R"


@dataclass
class PrimerHit:
    """One binding site of one primer on one target sequence."""
    primer_name: str
    primer_seq: str
    start: int                # 1-based, as reported by virtualPCR
    end: int
    orientation: Orientation
    identity: float           # percent, as reported
    tm: float                 # degrees C
    match_line: str = ""      # the '|||:|| ' line of the alignment
    target_context: str = ""  # target sequence around the binding site

    def _footprint(self) -> str:
        """Match symbols aligned one-per-base to the primer, in printed order.

        The parser stores the match line already left-aligned to the first
        printed primer base, so it is used without stripping: a blank in the
        first or last column is a real 5'/3' terminal mismatch, not padding, and
        stripping would hide exactly the terminal mismatches that matter most. It
        is clamped to the primer length to drop any trailing pad.
        """
        L = len(self.primer_seq)
        return self.match_line[:L] if L else self.match_line.strip()

    @property
    def mismatches(self) -> int:
        """Positions the aligner did not call a full match ('|')."""
        return sum(1 for c in self._footprint() if c != "|")

    @property
    def hard_mismatches(self) -> int:
        """Blanks only. ':' marks a degenerate/partial match, which is milder."""
        return self._footprint().count(" ")

    def mismatch_offsets_from_3prime(self) -> list[int]:
        """Distance of each mismatch from the primer 3' end (0 = terminal base).

        virtualPCR prints a forward primer 5'->3' (3' end rightmost) but a
        reverse primer 3'->5' (3' end leftmost), so the 3' column depends on
        orientation. The match line follows the printed strand either way.
        """
        m = self._footprint()
        if self.orientation == "reverse":
            return [i for i, c in enumerate(m) if c != "|"]
        return [len(m) - 1 - i for i, c in enumerate(m) if c != "|"]


@dataclass
class Amplicon:
    """A predicted PCR product.

    `size` is the span reported by virtualPCR (end - start + 1). The JAR's
    extracted `sequence` is one base short of that: it omits the last base of
    the reverse primer binding site. Use `sequence_complete` to detect it and
    `repair_sequence` (see parser) to restore the missing base from the target.
    """
    start: int
    end: int
    size: int
    ta: float | None = None   # annealing temperature, when reported
    sequence: str = ""

    @property
    def sequence_complete(self) -> bool:
        return not self.sequence or len(self.sequence) == self.size


@dataclass
class SequenceResult:
    """Outcome for one target sequence under one primer pair."""
    header: str                                   # full FASTA header, verbatim
    hits: list[PrimerHit] = field(default_factory=list)
    amplicons: list[Amplicon] = field(default_factory=list)
    total_amplicons: int = 0
    unique_sizes: int = 0

    @property
    def amplified(self) -> bool:
        return bool(self.amplicons)

    def hits_for(self, orientation: Orientation) -> list[PrimerHit]:
        return [h for h in self.hits if h.orientation == orientation]


@dataclass
class PairResult:
    """Outcome for one primer pair over one target file."""
    pair: PrimerPair
    target_file: str
    sequences: list[SequenceResult] = field(default_factory=list)
    primer_hit_counts: dict[str, int] = field(default_factory=dict)
    elapsed: str = ""
    # Headers present in the FASTA that produced no output block at all.
    silent_headers: list[str] = field(default_factory=list)

    @property
    def n_amplified(self) -> int:
        return sum(1 for s in self.sequences if s.amplified)


@dataclass
class RunConfig:
    """virtualPCR parameters. Names mirror the JAR's config keys."""
    min_len: int = 30
    max_len: int = 3000
    n_3prime_errors: int = 1
    probe: bool = False
    circular: bool = False
    ct_conversion: bool = False
    sequence_extract: bool = True
    show_primer_alignment: bool = True
    show_pcr_products: bool = True
    # Kept True: this is the output shape the parser expects (gives Ta + counts).
    show_only_amplicons: bool = True
    primer_statistic: bool = True

    def to_config_text(self, targets: str, primers: str, output: str) -> str:
        b = lambda v: "true" if v else "false"  # noqa: E731
        return "\n".join([
            f"targets_path={targets}",
            f"primers_path={primers}",
            f"output_path={output}",
            f"type={'probe' if self.probe else 'primer'}",
            f"molecular={'circle' if self.circular else 'linear'}",
            "linkedsearch=false",
            "FRpairs=false",
            f"CTconversion={b(self.ct_conversion)}",
            f"number3errors={self.n_3prime_errors}",
            f"minlen={self.min_len}",
            f"maxlen={self.max_len}",
            f"primerstatistic={b(self.primer_statistic)}",
            f"SequenceExtract={b(self.sequence_extract)}",
            f"ShowPrimerAlignment={b(self.show_primer_alignment)}",
            f"ShowPCRProducts={b(self.show_pcr_products)}",
            f"ShowOnlyAmplicons={b(self.show_only_amplicons)}",
            "ShowPrimerAlignmentPCRproduct=true",
            "",
        ])
