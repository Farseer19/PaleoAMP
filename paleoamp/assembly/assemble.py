"""Wrap MEGAHIT for metagenomic assembly with aDNA QC gating."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass
class AssemblyResult:
    sample_id: str
    contigs_path: Path
    n_contigs: int
    total_length_bp: int
    n50: int
    skipped: bool = False
    skip_reason: str = ""


def load_qc_verdicts(qc_report: Path) -> dict[str, bool]:
    """Return {sample_stem: pass_verdict} from a QC JSON report."""
    with open(qc_report) as fh:
        records = json.load(fh)
    return {
        _strip_fastq_extensions(Path(r["path"]).name): r["pass_verdict"]
        for r in records
    }


def _strip_fastq_extensions(filename: str) -> str:
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def _check_megahit() -> str:
    exe = shutil.which("megahit")
    if exe is None:
        raise RuntimeError(
            "megahit not found in PATH. Install with: conda install -c bioconda megahit"
        )
    return exe


def _find_reads(sample_dir: Path) -> tuple[list[Path], bool]:
    """
    Locate FASTQ files in a sample directory.

    Returns (read_paths, is_paired). Paired detection looks for _1/_2 or
    _R1/_R2 suffixes; everything else is treated as single-end.
    """
    fastqs = sorted(
        list(sample_dir.glob("*.fastq.gz"))
        + list(sample_dir.glob("*.fastq"))
        + list(sample_dir.glob("*.fq.gz"))
        + list(sample_dir.glob("*.fq"))
    )
    if not fastqs:
        return [], False

    r1 = [f for f in fastqs if any(_strip_fastq_extensions(f.name).endswith(s) for s in ("_1", "_R1"))]
    r2 = [f for f in fastqs if any(_strip_fastq_extensions(f.name).endswith(s) for s in ("_2", "_R2"))]
    if r1 and r2 and len(r1) == len(r2):
        return sorted(r1) + sorted(r2), True

    return fastqs, False


def _contig_stats(fasta: Path) -> dict:
    lengths = sorted(
        (len(r.seq) for r in SeqIO.parse(str(fasta), "fasta")),
        reverse=True,
    )
    if not lengths:
        return {"n": 0, "total": 0, "n50": 0}
    total = sum(lengths)
    cumulative = 0
    n50 = 0
    for ln in lengths:
        cumulative += ln
        if cumulative >= total / 2:
            n50 = ln
            break
    return {"n": len(lengths), "total": total, "n50": n50}


def assemble_sample(
    sample_dir: Path,
    output_dir: Path,
    sample_id: str,
    threads: int = 4,
    min_contig_len: int = 200,
) -> AssemblyResult:
    """Run MEGAHIT on a single sample directory and return assembly stats."""
    exe = _check_megahit()

    reads, paired = _find_reads(sample_dir)
    if not reads:
        return AssemblyResult(
            sample_id=sample_id,
            contigs_path=Path(),
            n_contigs=0,
            total_length_bp=0,
            n50=0,
            skipped=True,
            skip_reason="no FASTQ files found",
        )

    out_path = output_dir / sample_id
    if out_path.exists():
        shutil.rmtree(out_path)

    cmd = [exe]
    if paired:
        mid = len(reads) // 2
        cmd += ["-1", ",".join(str(r) for r in reads[:mid])]
        cmd += ["-2", ",".join(str(r) for r in reads[mid:])]
    else:
        cmd += ["-r", ",".join(str(r) for r in reads)]

    cmd += [
        "-o", str(out_path),
        "--min-contig-len", str(min_contig_len),
        "-t", str(threads),
        "--no-mercy",
    ]

    subprocess.run(cmd, check=True, capture_output=True, text=True)

    contigs = out_path / "final.contigs.fa"
    if not contigs.exists():
        raise RuntimeError(f"MEGAHIT completed but {contigs} not found")

    stats = _contig_stats(contigs)
    return AssemblyResult(
        sample_id=sample_id,
        contigs_path=contigs,
        n_contigs=stats["n"],
        total_length_bp=stats["total"],
        n50=stats["n50"],
    )


def assemble_all(
    reads_dir: Path,
    output_dir: Path,
    qc_report: Path,
    threads: int = 4,
    min_contig_len: int = 200,
) -> list[AssemblyResult]:
    """
    Assemble all samples under reads_dir that passed aDNA QC.

    Expects each sample in its own subdirectory named by accession (the layout
    produced by `paleoamp sra download`). Falls back to treating reads_dir
    itself as a single sample if no subdirectories are found.
    """
    verdicts = load_qc_verdicts(qc_report)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(d for d in reads_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        sample_dirs = [reads_dir]

    results: list[AssemblyResult] = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name

        fastqs = (
            list(sample_dir.glob("*.fastq.gz"))
            + list(sample_dir.glob("*.fastq"))
            + list(sample_dir.glob("*.fq.gz"))
            + list(sample_dir.glob("*.fq"))
        )
        stems = [_strip_fastq_extensions(f.name) for f in fastqs]
        passed = any(verdicts.get(s, False) for s in stems) or verdicts.get(sample_id, False)

        if not passed:
            results.append(AssemblyResult(
                sample_id=sample_id,
                contigs_path=Path(),
                n_contigs=0,
                total_length_bp=0,
                n50=0,
                skipped=True,
                skip_reason="failed aDNA QC",
            ))
            continue

        result = assemble_sample(
            sample_dir=sample_dir,
            output_dir=output_dir,
            sample_id=sample_id,
            threads=threads,
            min_contig_len=min_contig_len,
        )
        results.append(result)

    return results
