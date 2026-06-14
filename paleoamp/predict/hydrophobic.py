"""Hydrophobicity-based pre-filter for candidate AMP ORFs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Eisenberg consensus hydrophobicity scale — used for hydrophobic moment (µH).
# This is the correct scale for amphipathicity; Kyte-Doolittle mean is not
# because it averages both faces of an amphipathic helix to near zero.
_EISENBERG: dict[str, float] = {
    "A":  0.25, "R": -1.80, "N": -0.64, "D": -0.72, "C":  0.04,
    "Q": -0.69, "E": -0.62, "G":  0.16, "H": -0.40, "I":  0.73,
    "L":  0.53, "K": -1.10, "M":  0.26, "F":  0.61, "P": -0.07,
    "S": -0.26, "T": -0.18, "W":  0.37, "Y":  0.02, "V":  0.54,
}

# Residues considered hydrophobic for fraction calculation
_HYDROPHOBIC = frozenset("ACFILMVWY")


@dataclass
class HydrophobicityResult:
    record: SeqRecord
    h_moment: float
    hydrophobic_fraction: float
    passed: bool


def _hydrophobic_moment(seq: str, angle: float = 100.0) -> float:
    """Eisenberg hydrophobic moment for an α-helix (100° per residue)."""
    rad = math.radians(angle)
    sin_sum = sum(_EISENBERG.get(aa, 0.0) * math.sin(i * rad) for i, aa in enumerate(seq))
    cos_sum = sum(_EISENBERG.get(aa, 0.0) * math.cos(i * rad) for i, aa in enumerate(seq))
    return math.sqrt(sin_sum ** 2 + cos_sum ** 2) / len(seq)


def score_record(record: SeqRecord) -> HydrophobicityResult:
    seq = str(record.seq).upper().rstrip("*")
    if not seq:
        return HydrophobicityResult(record, 0.0, 0.0, False)

    h_moment = _hydrophobic_moment(seq)
    hydrophobic_fraction = sum(1 for aa in seq if aa in _HYDROPHOBIC) / len(seq)

    return HydrophobicityResult(
        record=record,
        h_moment=h_moment,
        hydrophobic_fraction=hydrophobic_fraction,
        passed=False,
    )


def filter_by_hydrophobicity(
    input_fasta: Path,
    output_fasta: Path,
    min_h_moment: float = 0.10,
    min_hydrophobic_fraction: float = 0.30,
) -> tuple[int, int]:
    """
    Filter ORFs that are unlikely to be membrane-active AMPs based on
    Eisenberg hydrophobic moment (µH) and hydrophobic residue fraction.

    Both thresholds must be satisfied. Returns (n_input, n_passed).

    µH captures amphipathicity — the spatial segregation of hydrophobic and
    hydrophilic faces — rather than net hydrophobicity.  A cationic AMP with
    many K/R residues can have a negative mean KD but a high µH.

    Defaults (min_h_moment=0.10, min_hydrophobic_fraction=0.30) are
    deliberately permissive — false negatives here are hard to recover.
    """
    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    passed: list[SeqRecord] = []

    for rec in records:
        r = score_record(rec)
        if r.h_moment >= min_h_moment and r.hydrophobic_fraction >= min_hydrophobic_fraction:
            passed.append(rec)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(passed, str(output_fasta), "fasta")
    return len(records), len(passed)


def score_fasta(input_fasta: Path) -> list[HydrophobicityResult]:
    """Return hydrophobicity scores for all records without filtering."""
    return [score_record(r) for r in SeqIO.parse(str(input_fasta), "fasta")]
