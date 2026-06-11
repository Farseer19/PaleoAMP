"""Download and manage AMP reference databases (DRAMP, APD3, dbAMP)."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Iterator

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


_DATABASES: dict[str, dict] = {
    "uniprot_amp": {
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=keyword:KW-0929+reviewed:true&format=fasta&size=500"
        ),
        "filename": "uniprot_amp_reviewed.fasta",
        "description": "UniProt reviewed antimicrobial peptides (KW-0929)",
        "paginated": True,
    },
    "dramp": {
        # DRAMP 3.0 bulk download — verify URL at https://dramp.cpu-bioinfor.org/downloads/
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=DRAMP3.0_general.fasta",
        "filename": "dramp_general.fasta",
        "description": "DRAMP 3.0 — general AMP sequences",
    },
    "apd3": {
        # APD3 moved to https://aps.unmc.edu — check for current release filename
        "url": "https://aps.unmc.edu/AP/lib/APD_sequence_release_09142020.fasta",
        "filename": "apd3.fasta",
        "description": "Antimicrobial Peptide Database (APD3)",
    },
    "dbamp": {
        "url": "https://awi.cuhk.edu.cn/dbAMP/download/dbAMP_v2.0.fasta",
        "filename": "dbamp_v2.fasta",
        "description": "dbAMP v2.0",
    },
}


def download_amp_databases(
    output_dir: Path,
    databases: list[str] | None = None,
    force: bool = False,
    timeout: int = 60,
) -> dict[str, Path]:
    """
    Download AMP reference databases to *output_dir*.

    Parameters
    ----------
    databases:
        Names from {"dramp", "apd3", "dbamp"}. Downloads all if None.
    force:
        Re-download even if the file already exists.

    Returns a mapping of database name → local FASTA path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    targets = databases or list(_DATABASES.keys())
    unknown = set(targets) - set(_DATABASES.keys())
    if unknown:
        raise ValueError(f"Unknown database(s): {unknown}. Choose from {list(_DATABASES)}")

    paths: dict[str, Path] = {}
    for name in targets:
        info = _DATABASES[name]
        dest = output_dir / info["filename"]

        if dest.exists() and not force:
            print(f"[skip] {name}: {dest} already exists (use --force to re-download)")
            paths[name] = dest
            continue

        print(f"[download] {name}: {info['description']}")
        try:
            _fetch_file(info["url"], dest, timeout=timeout)
            paths[name] = dest
            count = _count_sequences(dest)
            print(f"[ok] {name}: {count:,} sequences → {dest}")
        except Exception as exc:
            print(f"[error] {name}: {exc}", file=sys.stderr)

    return paths


def _fetch_file(url: str, dest: Path, timeout: int = 60) -> None:
    with requests.get(url, stream=True, timeout=timeout) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=65536):
                fh.write(chunk)


def _count_sequences(fasta_path: Path) -> int:
    return sum(1 for _ in SeqIO.parse(str(fasta_path), "fasta"))


def load_amp_sequences(
    db_paths: dict[str, Path],
) -> dict[str, list[SeqRecord]]:
    """Load all sequences from downloaded database files."""
    loaded: dict[str, list[SeqRecord]] = {}
    for name, path in db_paths.items():
        if not path.exists():
            print(f"[warn] {name}: file not found at {path}", file=sys.stderr)
            continue
        records = list(SeqIO.parse(str(path), "fasta"))
        loaded[name] = records
    return loaded


def merge_amp_databases(
    db_paths: dict[str, Path],
    output_path: Path,
    deduplicate: bool = True,
) -> Path:
    """
    Merge multiple AMP databases into a single non-redundant FASTA.

    Deduplication is by MD5 of the uppercase sequence string.
    """
    output_path = Path(output_path)
    seen: set[str] = set()
    written = 0

    with open(output_path, "w") as out:
        for db_name, path in db_paths.items():
            if not path.exists():
                continue
            for rec in SeqIO.parse(str(path), "fasta"):
                seq = str(rec.seq).upper()
                key = hashlib.md5(seq.encode()).hexdigest()
                if deduplicate and key in seen:
                    continue
                seen.add(key)
                rec.id = f"{db_name}|{rec.id}"
                rec.description = ""
                SeqIO.write(rec, out, "fasta")
                written += 1

    print(f"[merge] {written:,} unique AMP sequences → {output_path}")
    return output_path


def iter_local_fasta(path: Path) -> Iterator[SeqRecord]:
    """Yield SeqRecords from a local FASTA or gzipped FASTA file."""
    import gzip
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            yield from SeqIO.parse(fh, "fasta")
    else:
        yield from SeqIO.parse(str(path), "fasta")
