"""BLASTP-based screening of candidate ORFs against the merged AMP database."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

_BLAST_FMT6_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


@dataclass
class BlastResult:
    sample_id: str
    hits: pd.DataFrame      # all hits above threshold
    best_hits: pd.DataFrame # top hit per query
    n_queries: int
    n_hits: int
    tsv_path: Path


def _check_blast(tool: str) -> str:
    exe = shutil.which(tool)
    if exe is None:
        raise RuntimeError(
            f"{tool} not found in PATH. Install with: conda install -c bioconda blast"
        )
    return exe


def build_blast_db(fasta: Path, db_dir: Path | None = None) -> Path:
    """
    Run makeblastdb on *fasta* (protein) if the index files don't exist.

    Returns the database path prefix (passed to -db in blastp).
    """
    _check_blast("makeblastdb")
    db_path = (db_dir or fasta.parent) / fasta.stem
    index_file = db_path.with_suffix(".phr")  # created by makeblastdb

    if not index_file.exists():
        cmd = [
            "makeblastdb",
            "-in", str(fasta),
            "-dbtype", "prot",
            "-out", str(db_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"makeblastdb failed:\n{result.stderr.strip()}")

    return db_path


def run_blastp(
    query: Path,
    db: Path | str,
    output_dir: Path,
    sample_id: str,
    evalue: float = 1e-3,
    pident_min: float = 30.0,
    qcov_min: float = 50.0,
    threads: int = 4,
    remote: bool = False,
) -> BlastResult:
    """
    Run blastp of *query* against *db* and return filtered hits.

    remote=True  — query NCBI's remote servers; db should be an NCBI database
                   name such as "swissprot" or "nr". -num_threads is omitted
                   (remote jobs are server-side) and timeouts are longer.
    remote=False — local blastp against a makeblastdb-indexed file.
    """
    _check_blast("blastp")
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_tsv = output_dir / f"{sample_id}_blast_raw.tsv"
    filtered_tsv = output_dir / f"{sample_id}_blast_hits.tsv"

    cmd = [
        "blastp",
        "-query", str(query),
        "-db", str(db),
        "-out", str(raw_tsv),
        "-outfmt", "6",
        "-evalue", str(evalue),
        "-max_target_seqs", "5",
    ]
    if remote:
        cmd.append("-remote")
    else:
        cmd += ["-num_threads", str(threads)]

    timeout = 600 if remote else 120
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(f"blastp failed for {sample_id}:\n{result.stderr.strip()}")

    if not raw_tsv.exists() or raw_tsv.stat().st_size == 0:
        empty = pd.DataFrame(columns=_BLAST_FMT6_COLS + ["qcov"])
        return BlastResult(
            sample_id=sample_id,
            hits=empty,
            best_hits=empty,
            n_queries=0,
            n_hits=0,
            tsv_path=filtered_tsv,
        )

    df = pd.read_csv(raw_tsv, sep="\t", names=_BLAST_FMT6_COLS)
    df["qcov"] = (df["length"] / df["qend"].max()).clip(upper=1.0) * 100

    # Apply post-search filters
    df = df[(df["pident"] >= pident_min) & (df["qcov"] >= qcov_min)]

    # Best hit per query: lowest evalue, then highest bitscore
    best = (
        df.sort_values(["evalue", "bitscore"], ascending=[True, False])
        .drop_duplicates(subset=["qseqid"], keep="first")
        .reset_index(drop=True)
    )

    df.to_csv(filtered_tsv, sep="\t", index=False)

    # Count queries that had at least one hit
    n_queries = df["qseqid"].nunique()

    return BlastResult(
        sample_id=sample_id,
        hits=df,
        best_hits=best,
        n_queries=n_queries,
        n_hits=len(df),
        tsv_path=filtered_tsv,
    )
