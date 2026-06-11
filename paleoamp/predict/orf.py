"""Wrap Prodigal for ORF prediction on metagenomic contigs."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class ORFResult:
    sample_id: str
    protein_fasta: Path
    coords_gff: Path
    n_orfs: int
    skipped: bool = False
    skip_reason: str = ""


def _check_prodigal() -> str:
    exe = shutil.which("prodigal")
    if exe is None:
        raise RuntimeError(
            "prodigal not found in PATH. Install with: conda install -c bioconda prodigal"
        )
    return exe


def predict_orfs(
    contigs: Path,
    output_dir: Path,
    sample_id: str,
) -> ORFResult:
    """
    Run Prodigal in metagenome mode on *contigs* and write protein FASTA + GFF.

    Uses -p meta so no per-sample training is needed — appropriate for
    short, fragmented ancient metagenome assemblies.
    """
    exe = _check_prodigal()
    output_dir.mkdir(parents=True, exist_ok=True)

    protein_fa = output_dir / f"{sample_id}.faa"
    coords_gff = output_dir / f"{sample_id}.gff"

    cmd = [
        exe,
        "-i", str(contigs),
        "-a", str(protein_fa),
        "-f", "gff",
        "-o", str(coords_gff),
        "-p", "meta",
        "-q",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Prodigal failed for {sample_id}:\n{result.stderr.strip()}"
        )

    n_orfs = sum(1 for _ in SeqIO.parse(str(protein_fa), "fasta"))

    return ORFResult(
        sample_id=sample_id,
        protein_fasta=protein_fa,
        coords_gff=coords_gff,
        n_orfs=n_orfs,
    )


def predict_all(
    assembly_dir: Path,
    output_dir: Path,
) -> list[ORFResult]:
    """
    Run ORF prediction on every sample assembly found under *assembly_dir*.

    Expects the layout produced by assemble_all(): one subdirectory per
    sample containing final.contigs.fa.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(d for d in assembly_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        raise FileNotFoundError(f"No sample directories found in {assembly_dir}")

    results: list[ORFResult] = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        contigs = sample_dir / "final.contigs.fa"

        if not contigs.exists():
            results.append(ORFResult(
                sample_id=sample_id,
                protein_fasta=Path(),
                coords_gff=Path(),
                n_orfs=0,
                skipped=True,
                skip_reason="final.contigs.fa not found",
            ))
            continue

        result = predict_orfs(
            contigs=contigs,
            output_dir=output_dir / sample_id,
            sample_id=sample_id,
        )
        results.append(result)

    return results


def load_orfs(protein_fasta: Path) -> list[SeqRecord]:
    return list(SeqIO.parse(str(protein_fasta), "fasta"))
