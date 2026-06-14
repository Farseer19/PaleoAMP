"""
Redundancy reduction and cluster-based train/test split.

Clusters the combined dataset at 40% sequence identity using MMseqs2
(preferred) or CD-HIT, then assigns whole clusters to train/test so
no sequence in test is >40% identical to any sequence in train.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _check_tool(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        raise RuntimeError(
            f"{name} not found. Install with: conda install -c bioconda {name}"
        )
    return exe


def _df_to_fasta(df: pd.DataFrame, fasta_path: Path) -> None:
    records = [
        SeqRecord(
            seq=__import__("Bio.Seq", fromlist=["Seq"]).Seq(row.sequence),
            id=str(i),
            description="",
        )
        for i, row in df.iterrows()
    ]
    SeqIO.write(records, str(fasta_path), "fasta")


def cluster_mmseqs2(
    df: pd.DataFrame,
    work_dir: Path,
    min_seq_id: float = 0.40,
    threads: int = 4,
) -> pd.Series:
    """
    Cluster sequences using MMseqs2 easy-cluster at *min_seq_id* identity.

    Returns a pd.Series indexed by df.index with integer cluster IDs.
    """
    _check_tool("mmseqs")
    work_dir.mkdir(parents=True, exist_ok=True)

    fasta = work_dir / "input.fasta"
    _df_to_fasta(df, fasta)

    prefix = work_dir / "clusters"
    tmp = work_dir / "tmp"
    tmp.mkdir(exist_ok=True)

    cmd = [
        "mmseqs", "easy-cluster",
        str(fasta), str(prefix), str(tmp),
        "--min-seq-id", str(min_seq_id),
        "--cov-mode", "0",
        "-c", "0.8",
        "--threads", str(threads),
        "-v", "0",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MMseqs2 failed:\n{result.stderr.strip()}")

    # Parse cluster TSV: rep_id \t member_id
    cluster_tsv = work_dir / "clusters_cluster.tsv"
    cluster_map: dict[str, int] = {}
    rep_to_id: dict[str, int] = {}

    with open(cluster_tsv) as fh:
        for line in fh:
            rep, member = line.strip().split("\t")
            if rep not in rep_to_id:
                rep_to_id[rep] = len(rep_to_id)
            cluster_map[member] = rep_to_id[rep]

    return pd.Series([cluster_map[str(i)] for i in df.index], index=df.index, name="cluster_id")


def cluster_cdhit(
    df: pd.DataFrame,
    work_dir: Path,
    min_seq_id: float = 0.40,
    threads: int = 4,
) -> pd.Series:
    """
    Cluster sequences using CD-HIT at *min_seq_id* identity.

    Returns a pd.Series indexed by df.index with integer cluster IDs.
    """
    _check_tool("cd-hit")
    work_dir.mkdir(parents=True, exist_ok=True)

    fasta = work_dir / "input.fasta"
    _df_to_fasta(df, fasta)
    out = work_dir / "cdhit_out"

    # CD-HIT word size: 2 for <0.5, 3 for 0.5-0.6, 4 for 0.6-0.7, 5 for ≥0.7
    word_size = 2 if min_seq_id < 0.5 else (3 if min_seq_id < 0.6 else (4 if min_seq_id < 0.7 else 5))

    cmd = [
        "cd-hit",
        "-i", str(fasta),
        "-o", str(out),
        "-c", str(min_seq_id),
        "-n", str(word_size),
        "-T", str(threads),
        "-M", "0",
        "-d", "0",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"CD-HIT failed:\n{result.stderr.strip()}")

    # Parse .clstr file
    clstr_path = Path(str(out) + ".clstr")
    cluster_map: dict[str, int] = {}
    cluster_id = -1

    with open(clstr_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id += 1
            else:
                m = __import__("re").search(r">(\S+?)\.\.\..*", line)
                if m:
                    cluster_map[m.group(1)] = cluster_id

    # Sequences absent from the clstr file (filtered by CD-HIT) get singleton clusters.
    next_id = cluster_id + 1
    result = []
    for i in df.index:
        key = str(i)
        if key not in cluster_map:
            cluster_map[key] = next_id
            next_id += 1
        result.append(cluster_map[key])

    return pd.Series(result, index=df.index, name="cluster_id")


def cluster_sequences(
    df: pd.DataFrame,
    work_dir: Path,
    min_seq_id: float = 0.40,
    threads: int = 4,
) -> pd.Series:
    """Try MMseqs2 first, fall back to CD-HIT."""
    if shutil.which("mmseqs"):
        print(f"[cluster] Using MMseqs2 at {min_seq_id:.0%} identity")
        return cluster_mmseqs2(df, work_dir, min_seq_id, threads)
    elif shutil.which("cd-hit"):
        print(f"[cluster] Using CD-HIT at {min_seq_id:.0%} identity")
        return cluster_cdhit(df, work_dir, min_seq_id, threads)
    else:
        raise RuntimeError(
            "Neither mmseqs nor cd-hit found. Install one with:\n"
            "  conda install -c bioconda mmseqs2\n"
            "  conda install -c bioconda cd-hit"
        )


def cluster_split(
    df: pd.DataFrame,
    cluster_ids: pd.Series,
    test_fraction: float = 0.20,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Assign whole clusters to train (80%) or test (20%) splits.

    This guarantees no sequence in the test set shares >40% identity
    with any training sequence — the step most published tools skip.

    Adds a 'split' column ('train' or 'test') to a copy of *df*.
    """
    rng = np.random.default_rng(random_seed)

    df = df.copy()
    df["cluster_id"] = cluster_ids

    unique_clusters = df["cluster_id"].unique()
    rng.shuffle(unique_clusters)

    n_test = max(1, int(len(unique_clusters) * test_fraction))
    test_clusters = set(unique_clusters[:n_test])

    df["split"] = df["cluster_id"].apply(lambda c: "test" if c in test_clusters else "train")

    n_train = (df["split"] == "train").sum()
    n_test_seqs = (df["split"] == "test").sum()
    n_clusters = len(unique_clusters)
    print(
        f"[split] {n_clusters:,} clusters → "
        f"train: {n_train:,} seqs | test: {n_test_seqs:,} seqs"
    )
    pos_train = ((df["split"] == "train") & (df["label"] == 1)).sum()
    pos_test = ((df["split"] == "test") & (df["label"] == 1)).sum()
    print(f"        positives — train: {pos_train:,} | test: {pos_test:,}")

    return df
