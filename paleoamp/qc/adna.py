"""
Ancient DNA damage quality assessment.

Computes C→T frequencies at the 5' end and G→A frequencies at the 3' end
of sequencing reads — the hallmark of post-mortem cytosine deamination.
Also collects standard QC metrics (read lengths, mean Phred scores).

Works directly on FASTQ/FASTQ.gz files without needing alignment,
using read-terminal nucleotide composition as a proxy for damage.
"""

from __future__ import annotations

import gzip
import math
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator


@dataclass
class DamageProfile:
    """Per-file ancient DNA damage and basic QC metrics."""

    path: str
    total_reads: int = 0
    passed_length_filter: int = 0

    # Positional base counts: [pos] → {"A": n, "C": n, "G": n, "T": n, "N": n}
    # 5′ end (left side of read)
    five_prime_counts: list[dict[str, int]] = field(default_factory=list)
    # 3′ end (right side of read, reversed so index 0 = last base)
    three_prime_counts: list[dict[str, int]] = field(default_factory=list)

    read_lengths: list[int] = field(default_factory=list)
    mean_phred_scores: list[float] = field(default_factory=list)

    # Populated by finalise()
    ct_rates_5prime: list[float] = field(default_factory=list)
    ga_rates_3prime: list[float] = field(default_factory=list)
    mean_read_length: float = 0.0
    median_read_length: float = 0.0
    mean_phred: float = 0.0
    pass_verdict: bool = False
    fail_reasons: list[str] = field(default_factory=list)

    def finalise(self, thresholds: "QCThresholds") -> None:
        """Compute derived metrics and apply quality thresholds."""
        if not self.read_lengths:
            self.fail_reasons.append("no reads")
            return

        self.mean_read_length = statistics.mean(self.read_lengths)
        self.median_read_length = statistics.median(self.read_lengths)

        if self.mean_phred_scores:
            self.mean_phred = statistics.mean(self.mean_phred_scores)

        # C→T rates at each 5' position
        self.ct_rates_5prime = _compute_substitution_rates(
            self.five_prime_counts, from_base="C", to_base="T"
        )
        # G→A rates at each 3' position
        self.ga_rates_3prime = _compute_substitution_rates(
            self.three_prime_counts, from_base="G", to_base="A"
        )

        self._apply_thresholds(thresholds)

    def _apply_thresholds(self, t: "QCThresholds") -> None:
        reasons: list[str] = []

        if self.mean_phred < t.min_mean_phred:
            reasons.append(
                f"mean Phred {self.mean_phred:.1f} < {t.min_mean_phred}"
            )

        if self.mean_read_length < t.min_read_length:
            reasons.append(
                f"mean read length {self.mean_read_length:.0f} < {t.min_read_length}"
            )

        if self.mean_read_length > t.max_read_length:
            reasons.append(
                f"mean read length {self.mean_read_length:.0f} > {t.max_read_length} bp "
                "(possible modern contamination — aDNA is typically short)"
            )

        frac_passing = (
            self.passed_length_filter / self.total_reads
            if self.total_reads > 0 else 0.0
        )
        if frac_passing < t.min_passing_fraction:
            reasons.append(
                f"only {frac_passing:.1%} of reads pass length filter "
                f"(min {t.min_passing_fraction:.1%})"
            )

        ct_pos1 = self.ct_rates_5prime[0] if self.ct_rates_5prime else 0.0
        ga_pos1 = self.ga_rates_3prime[0] if self.ga_rates_3prime else 0.0

        ct_ok = ct_pos1 >= t.min_ct_rate_5prime
        ga_ok = ga_pos1 >= t.min_ga_rate_3prime

        if t.require_both_termini:
            if not ct_ok:
                reasons.append(
                    f"5′ C→T rate {ct_pos1:.3f} < {t.min_ct_rate_5prime} (no aDNA damage signal)"
                )
            if not ga_ok:
                reasons.append(
                    f"3′ G→A rate {ga_pos1:.3f} < {t.min_ga_rate_3prime} (no aDNA damage signal)"
                )
        else:
            if not (ct_ok or ga_ok):
                reasons.append(
                    f"neither 5′ C→T ({ct_pos1:.3f}) nor 3′ G→A ({ga_pos1:.3f}) "
                    "meet aDNA damage thresholds"
                )

        self.fail_reasons = reasons
        self.pass_verdict = len(reasons) == 0


@dataclass
class QCThresholds:
    min_mean_phred: float = 20.0
    min_read_length: int = 30
    max_read_length: int = 300
    min_passing_fraction: float = 0.5
    min_ct_rate_5prime: float = 0.10
    min_ga_rate_3prime: float = 0.10
    terminal_positions: int = 10
    require_both_termini: bool = True

    @classmethod
    def from_config(cls, cfg: dict) -> "QCThresholds":
        q = cfg.get("quality", {})
        d = cfg.get("adna_damage", {})
        return cls(
            min_mean_phred=q.get("min_mean_phred", 20.0),
            min_read_length=q.get("min_read_length", 30),
            max_read_length=q.get("max_read_length", 300),
            min_passing_fraction=q.get("min_passing_fraction", 0.5),
            min_ct_rate_5prime=d.get("min_ct_rate_5prime", 0.05),
            min_ga_rate_3prime=d.get("min_ga_rate_3prime", 0.05),
            terminal_positions=d.get("terminal_positions", 10),
            require_both_termini=d.get("require_both_termini", True),
        )


def assess_fastq(
    fastq_path: Path,
    thresholds: QCThresholds | None = None,
    max_reads: int | None = None,
) -> DamageProfile:
    """
    Parse *fastq_path* and return a DamageProfile with aDNA quality metrics.

    Parameters
    ----------
    fastq_path:
        Path to a FASTQ or FASTQ.gz file.
    thresholds:
        QCThresholds instance; defaults to QCThresholds() if not provided.
    max_reads:
        Cap the number of reads processed (useful for large files).
    """
    thresholds = thresholds or QCThresholds()
    n = thresholds.terminal_positions

    profile = DamageProfile(path=str(fastq_path))
    profile.five_prime_counts = [_empty_counts() for _ in range(n)]
    profile.three_prime_counts = [_empty_counts() for _ in range(n)]

    for i, (seq, quals) in enumerate(_iter_fastq(fastq_path)):
        if max_reads and i >= max_reads:
            break

        profile.total_reads += 1
        read_len = len(seq)
        profile.read_lengths.append(read_len)

        if quals:
            profile.mean_phred_scores.append(
                sum(quals) / len(quals)
            )

        if read_len < thresholds.min_read_length:
            continue
        profile.passed_length_filter += 1

        # 5′ terminal positions
        for pos in range(min(n, read_len)):
            base = seq[pos].upper()
            if base in profile.five_prime_counts[pos]:
                profile.five_prime_counts[pos][base] += 1

        # 3′ terminal positions (reversed: index 0 = last base of read)
        for pos in range(min(n, read_len)):
            base = seq[-(pos + 1)].upper()
            if base in profile.three_prime_counts[pos]:
                profile.three_prime_counts[pos][base] += 1

    profile.finalise(thresholds)
    return profile


def assess_fastq_directory(
    directory: Path,
    thresholds: QCThresholds | None = None,
    max_reads: int | None = None,
    recursive: bool = False,
) -> list[DamageProfile]:
    """Assess all FASTQ / FASTQ.gz files in *directory*.

    Set *recursive=True* to search subdirectories (e.g. data/reads/ containing
    per-accession subdirs). Files are sorted so output order is deterministic.
    """
    directory = Path(directory)
    patterns = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    if recursive:
        files = sorted(f for pat in patterns for f in directory.rglob(pat))
    else:
        files = sorted(f for pat in patterns for f in directory.glob(pat))
    if not files:
        loc = "in or under" if recursive else "in"
        raise FileNotFoundError(f"No FASTQ files found {loc} {directory}")

    return [
        assess_fastq(f, thresholds=thresholds, max_reads=max_reads)
        for f in files
    ]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _empty_counts() -> dict[str, int]:
    return {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}


def _compute_substitution_rates(
    positional_counts: list[dict[str, int]],
    from_base: str,
    to_base: str,
) -> list[float]:
    """
    For each position, compute to_base / (from_base + to_base).

    This is the misincorporation rate assuming all from_base → to_base changes
    are damage events. Returns 0.0 when coverage is zero.
    """
    rates = []
    for counts in positional_counts:
        f = counts.get(from_base, 0)
        t = counts.get(to_base, 0)
        total = f + t
        rates.append(t / total if total > 0 else 0.0)
    return rates


def _iter_fastq(path: Path) -> Iterator[tuple[str, list[int]]]:
    """Yield (sequence, phred_scores) tuples from a FASTQ or FASTQ.gz file."""
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open

    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip()
            fh.readline()  # + line
            qual_str = fh.readline().rstrip()

            if not seq:
                continue

            phred = [ord(c) - 33 for c in qual_str] if qual_str else []
            yield seq, phred
