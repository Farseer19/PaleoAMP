"""
Novelty screening of high-scoring AMP candidates via MMseqs2.

Sequences with no hit in the modern AMP database at ≥ 40% identity
and ≥ 80% query coverage are flagged as candidate novel ancient AMPs.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_MMSEQS_RESULT_COLS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
    "qlen", "tlen",
]


@dataclass
class NoveltyResult:
    sample_id: str
    n_candidates: int     # sequences passed to MMseqs2
    n_hits: int           # sequences with a match above threshold
    n_novel: int          # sequences with no match (novel)
    hits: pd.DataFrame    # all hits above threshold
    novel: pd.DataFrame   # novel candidates (no hit)
    tsv_path: Path


def _check_mmseqs() -> str:
    exe = shutil.which("mmseqs")
    if exe is None:
        raise RuntimeError(
            "mmseqs not found. Install with: conda install -c bioconda mmseqs2"
        )
    return exe


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    seqrecs = [SeqRecord(Seq(seq), id=sid, description="") for sid, seq in records]
    SeqIO.write(seqrecs, str(path), "fasta")


def screen_novelty(
    scores_tsv: Path,
    amp_db_fasta: Path,
    output_dir: Path,
    sample_id: str,
    score_threshold: float = 0.5,
    pident_threshold: float = 40.0,
    qcov_threshold: float = 80.0,
    threads: int = 4,
) -> NoveltyResult:
    """
    Screen high-scoring AMP candidates against *amp_db_fasta* with MMseqs2.

    Candidates are taken from *scores_tsv* (output of 'paleoamp ml score')
    where amp_score >= *score_threshold*.

    A candidate is "novel" if it has no hit with pident >= *pident_threshold*
    AND qcov >= *qcov_threshold*.  These thresholds are intentionally strict
    (40% identity, 80% coverage) to ensure genuine novelty relative to all
    known modern AMP sequences.
    """
    _check_mmseqs()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load candidates above score threshold
    df_scores = pd.read_csv(scores_tsv, sep="\t")
    candidates = df_scores[df_scores["amp_score"] >= score_threshold][["seq_id", "sequence", "amp_score"]]

    if candidates.empty:
        empty = pd.DataFrame(columns=["seq_id", "sequence", "amp_score"])
        return NoveltyResult(
            sample_id=sample_id,
            n_candidates=0, n_hits=0, n_novel=0,
            hits=empty, novel=empty,
            tsv_path=output_dir / f"{sample_id}_novel.tsv",
        )

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        query_fa  = tmp / "query.fasta"
        db_path   = tmp / "ampdb"
        query_db  = tmp / "querydb"
        result_db = tmp / "resultdb"
        result_m8 = output_dir / f"{sample_id}_mmseqs_hits.tsv"

        _write_fasta(list(candidates[["seq_id", "sequence"]].itertuples(index=False, name=None)), query_fa)

        # Build MMseqs2 target database from AMP FASTA
        subprocess.run(
            ["mmseqs", "createdb", str(amp_db_fasta), str(db_path)],
            capture_output=True, check=True,
        )
        subprocess.run(
            ["mmseqs", "createdb", str(query_fa), str(query_db)],
            capture_output=True, check=True,
        )

        # Search
        subprocess.run(
            [
                "mmseqs", "search",
                str(query_db), str(db_path), str(result_db), str(tmp / "searchtmp"),
                "--min-seq-id", str(pident_threshold / 100),
                "-c", str(qcov_threshold / 100),
                "--cov-mode", "0",
                "--threads", str(threads),
                "-v", "0",
            ],
            capture_output=True, check=True,
        )

        # Convert to tabular m8 format with qlen/tlen
        subprocess.run(
            [
                "mmseqs", "convertalis",
                str(query_db), str(db_path), str(result_db), str(result_m8),
                "--format-output", "query,target,pident,alnlen,mismatch,gapopen,"
                                   "qstart,qend,tstart,tend,evalue,bits,qlen,tlen",
            ],
            capture_output=True, check=True,
        )

    # Parse hits
    if result_m8.exists() and result_m8.stat().st_size > 0:
        hits = pd.read_csv(result_m8, sep="\t", names=_MMSEQS_RESULT_COLS)
        hits["qcov"] = (hits["qend"] - hits["qstart"] + 1) / hits["qlen"] * 100
        hits = hits[(hits["pident"] * 100 >= pident_threshold) & (hits["qcov"] >= qcov_threshold)]
    else:
        hits = pd.DataFrame(columns=_MMSEQS_RESULT_COLS + ["qcov"])

    hit_ids = set(hits["query"].unique()) if not hits.empty else set()

    novel = candidates[~candidates["seq_id"].isin(hit_ids)].copy()
    novel = novel.sort_values("amp_score", ascending=False).reset_index(drop=True)

    novel_path = output_dir / f"{sample_id}_novel.tsv"
    novel.to_csv(novel_path, sep="\t", index=False)

    hits_with_scores = hits.merge(
        candidates[["seq_id", "amp_score"]],
        left_on="query", right_on="seq_id", how="left",
    )
    hits_path = output_dir / f"{sample_id}_known_hits.tsv"
    hits_with_scores.to_csv(hits_path, sep="\t", index=False)

    return NoveltyResult(
        sample_id=sample_id,
        n_candidates=len(candidates),
        n_hits=len(hit_ids),
        n_novel=len(novel),
        hits=hits_with_scores,
        novel=novel,
        tsv_path=novel_path,
    )
