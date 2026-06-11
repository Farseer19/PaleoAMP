"""Hydrophobicity-based pre-filter for candidate AMP ORFs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Kyte-Doolittle hydrophobicity scale
_KD: dict[str, float] = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Residues considered hydrophobic for fraction calculation
_HYDROPHOBIC = frozenset("ACFILMVWY")


@dataclass
class HydrophobicityResult:
    record: SeqRecord
    mean_kd: float
    hydrophobic_fraction: float
    passed: bool


def score_record(record: SeqRecord) -> HydrophobicityResult:
    seq = str(record.seq).upper().rstrip("*")
    if not seq:
        return HydrophobicityResult(record, 0.0, 0.0, False)

    kd_values = [_KD[aa] for aa in seq if aa in _KD]
    mean_kd = sum(kd_values) / len(kd_values) if kd_values else 0.0
    hydrophobic_fraction = sum(1 for aa in seq if aa in _HYDROPHOBIC) / len(seq)

    return HydrophobicityResult(
        record=record,
        mean_kd=mean_kd,
        hydrophobic_fraction=hydrophobic_fraction,
        passed=False,
    )


def filter_by_hydrophobicity(
    input_fasta: Path,
    output_fasta: Path,
    min_kd: float = 0.0,
    min_hydrophobic_fraction: float = 0.30,
) -> tuple[int, int]:
    """
    Filter ORFs that are unlikely to be membrane-active AMPs based on
    Kyte-Doolittle mean hydrophobicity and hydrophobic residue fraction.

    Both thresholds must be satisfied. Returns (n_input, n_passed).

    Defaults (min_kd=0.0, min_hydrophobic_fraction=0.30) are deliberately
    permissive — AMPs span a wide hydrophobicity range and false negatives
    here are hard to recover downstream.
    """
    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    passed: list[SeqRecord] = []

    for rec in records:
        r = score_record(rec)
        if r.mean_kd >= min_kd and r.hydrophobic_fraction >= min_hydrophobic_fraction:
            passed.append(rec)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(passed, str(output_fasta), "fasta")
    return len(records), len(passed)


def score_fasta(input_fasta: Path) -> list[HydrophobicityResult]:
    """Return hydrophobicity scores for all records without filtering."""
    return [score_record(r) for r in SeqIO.parse(str(input_fasta), "fasta")]
