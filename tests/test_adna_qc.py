"""Tests for ancient DNA damage quality assessment."""

import gzip
import tempfile
from pathlib import Path

import pytest

from paleoamp.qc.adna import (
    QCThresholds,
    DamageProfile,
    assess_fastq,
    _compute_substitution_rates,
    _empty_counts,
)


def _make_fastq(records: list[tuple[str, str]], gz: bool = False) -> Path:
    """Write a temporary FASTQ file and return its path."""
    tmp = tempfile.NamedTemporaryFile(
        suffix=".fastq.gz" if gz else ".fastq",
        delete=False,
        mode="wb" if gz else "w",
    )
    content = ""
    for i, (seq, qual) in enumerate(records):
        content += f"@read_{i}\n{seq}\n+\n{qual}\n"

    if gz:
        tmp.write(gzip.compress(content.encode()))
    else:
        tmp.write(content)
    tmp.flush()
    return Path(tmp.name)


class TestComputeSubstitutionRates:
    def test_basic_ct_rate(self):
        counts = [{"A": 0, "C": 80, "G": 0, "T": 20, "N": 0}]
        rates = _compute_substitution_rates(counts, "C", "T")
        assert abs(rates[0] - 0.20) < 1e-9

    def test_zero_coverage(self):
        counts = [{"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}]
        rates = _compute_substitution_rates(counts, "C", "T")
        assert rates[0] == 0.0

    def test_multiple_positions(self):
        counts = [
            {"A": 0, "C": 90, "G": 0, "T": 10, "N": 0},
            {"A": 0, "C": 95, "G": 0, "T": 5, "N": 0},
        ]
        rates = _compute_substitution_rates(counts, "C", "T")
        assert abs(rates[0] - 0.10) < 1e-9
        assert abs(rates[1] - 0.05) < 1e-9


class TestAssessFastq:
    def _make_damaged_reads(self, n: int = 200) -> list[tuple[str, str]]:
        """
        Generate synthetic reads with aDNA-like damage:
        first base is biased toward T (from C→T), last base toward A (from G→A).
        Mean length ~60 bp (ancient-like).
        """
        records = []
        for i in range(n):
            # Alternate a couple patterns to produce measurable signal
            if i % 5 == 0:
                seq = "T" + "ACGT" * 14 + "A"  # 5' T, 3' A → mimics C→T and G→A
            else:
                seq = "C" + "ACGT" * 14 + "G"  # normal
            qual = "I" * len(seq)  # Phred 40
            records.append((seq, qual))
        return records

    def test_basic_pass(self):
        # High damage signal: all reads start T (C→T) and end A (G→A)
        records = [("T" + "ACGT" * 10 + "A", "I" * 42)] * 300
        path = _make_fastq(records)
        t = QCThresholds(
            min_ct_rate_5prime=0.01,
            min_ga_rate_3prime=0.01,
            min_mean_phred=5,
            min_read_length=20,
            min_passing_fraction=0.1,
        )
        profile = assess_fastq(path, thresholds=t)
        assert profile.total_reads == 300
        assert profile.pass_verdict, profile.fail_reasons

    def test_short_reads_fail_length_filter(self):
        records = [("ACGT", "IIII")] * 100
        path = _make_fastq(records)
        t = QCThresholds(min_read_length=30)
        profile = assess_fastq(path, thresholds=t)
        assert profile.passed_length_filter == 0
        assert not profile.pass_verdict

    def test_low_damage_fails(self):
        # Reads with no C→T or G→A at terminals
        records = [("A" + "ACGT" * 10 + "T", "I" * 42)] * 300
        path = _make_fastq(records)
        t = QCThresholds(
            min_ct_rate_5prime=0.05,
            min_ga_rate_3prime=0.05,
            require_both_termini=True,
        )
        profile = assess_fastq(path, thresholds=t)
        assert not profile.pass_verdict
        assert any("C→T" in r or "G→A" in r for r in profile.fail_reasons)

    def test_gzipped_fastq(self):
        records = [("T" + "ACGT" * 10 + "A", "I" * 42)] * 100
        path = _make_fastq(records, gz=True)
        t = QCThresholds(
            min_ct_rate_5prime=0.01,
            min_ga_rate_3prime=0.01,
            min_mean_phred=5,
            min_read_length=20,
            min_passing_fraction=0.1,
        )
        profile = assess_fastq(path, thresholds=t)
        assert profile.total_reads == 100

    def test_max_reads_cap(self):
        records = [("TACGTA" * 8, "I" * 48)] * 1000
        path = _make_fastq(records)
        profile = assess_fastq(path, max_reads=50)
        assert profile.total_reads == 50

    def test_ct_and_ga_rates_are_computed(self):
        # 50 % of reads start with T at pos 0 (half of C→T + C pool)
        records = (
            [("T" + "ACGT" * 10, "I" * 41)] * 50
            + [("C" + "ACGT" * 10, "I" * 41)] * 50
        )
        path = _make_fastq(records)
        t = QCThresholds(min_read_length=10, min_passing_fraction=0.0,
                         min_mean_phred=0)
        profile = assess_fastq(path, thresholds=t)
        ct = profile.ct_rates_5prime[0]
        assert abs(ct - 0.5) < 0.05


class TestQCThresholds:
    def test_from_config_defaults(self):
        t = QCThresholds.from_config({})
        assert t.min_mean_phred == 20.0
        assert t.min_ct_rate_5prime == 0.05

    def test_from_config_override(self):
        cfg = {
            "quality": {"min_mean_phred": 25},
            "adna_damage": {"min_ct_rate_5prime": 0.10},
        }
        t = QCThresholds.from_config(cfg)
        assert t.min_mean_phred == 25
        assert t.min_ct_rate_5prime == 0.10
