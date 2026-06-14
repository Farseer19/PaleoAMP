"""Download and manage AMP reference databases for the novelty screen."""

from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path
from typing import Iterator

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


_STANDARD_AAS = frozenset("ACDEFGHIKLMNPQRSTVWY")

_DATABASES: dict[str, dict] = {
    "apd": {
        "url": "https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta",
        "filename": "apd_natural.fasta",
        "description": "APD2024a — natural antimicrobial peptides",
    },
    "dramp_general": {
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta",
        "filename": "dramp_general.fasta",
        "description": "DRAMP 3.0 — general AMP sequences",
    },
    "dramp_specific": {
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/specific_amps.fasta",
        "filename": "dramp_specific.fasta",
        "description": "DRAMP 3.0 — specific AMPs (target-activity annotations)",
    },
    "dramp_predicted": {
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/predicted_amps.fasta",
        "filename": "dramp_predicted.fasta",
        "description": "DRAMP 3.0 — computationally predicted AMPs (~60k sequences, broader novelty coverage)",
    },
    "uniprot_amp": {
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=reviewed:true+keyword:KW-0929+length:[5+TO+200]&format=fasta&size=500"
        ),
        "filename": "uniprot_amp_reviewed.fasta",
        "description": "UniProt reviewed antimicrobial peptides (KW-0929)",
        "paginated": True,
        "paginate": 10000,
    },
    "dbamp": {
        "url": "https://awi.cuhk.edu.cn/dbAMP/download/dbAMP_v2.0.fasta",
        "filename": "dbamp_v2.fasta",
        "description": "dbAMP v2.0",
        "ssl_verify": False,  # server SSL certificate issue — skip verification
    },
    "dbaasp": {
        "url": "https://dbaasp.org/api/peptides/?format=fasta",
        "filename": "dbaasp.fasta",
        "description": "DBAASP v3 — database of antimicrobial activity and structure",
    },
    "bactibase": {
        "url": "http://bactibase.hammamilab.org/peptides.fasta",
        "filename": "bactibase.fasta",
        "description": "BACTIBASE — experimentally validated bacteriocins",
    },
    "uniprot_amp_bacteria": {
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=reviewed:true+keyword:KW-0929+taxonomy_id:2+length:[5+TO+200]"
            "&format=fasta&size=500"
        ),
        "filename": "uniprot_amp_bacteria.fasta",
        "description": "UniProt reviewed bacterial AMP peptides (KW-0929 + taxonomy:Bacteria)",
        "paginated": True,
        "paginate": 10000,
    },
    "uniprot_amp_archaea": {
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=reviewed:true+keyword:KW-0929+taxonomy_id:2157+length:[5+TO+200]"
            "&format=fasta&size=500"
        ),
        "filename": "uniprot_amp_archaea.fasta",
        "description": "UniProt reviewed archaeal AMP peptides (KW-0929 + taxonomy:Archaea)",
        "paginated": True,
        "paginate": 5000,
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
            if info.get("paginated"):
                _fetch_paginated(info["url"], dest,
                                 total=info.get("paginate", 10000),
                                 timeout=timeout)
            else:
                _fetch_file(info["url"], dest, timeout=timeout,
                            ssl_verify=info.get("ssl_verify", True))
            paths[name] = dest
            count = _count_sequences(dest)
            print(f"[ok] {name}: {count:,} sequences → {dest}")
        except Exception as exc:
            print(f"[error] {name}: {exc}", file=sys.stderr)

    return paths


def _fetch_file(url: str, dest: Path, timeout: int = 60, ssl_verify: bool = True) -> None:
    with requests.get(url, stream=True, timeout=timeout, verify=ssl_verify) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=65536):
                fh.write(chunk)


def _fetch_paginated(url: str, dest: Path, total: int, timeout: int = 60,
                     max_retries: int = 5) -> None:
    """Fetch a UniProt FASTA result set across multiple pages until *total* seqs collected."""
    import re as _re
    collected = 0
    next_url: str | None = url
    with open(dest, "w") as fh:
        while next_url and collected < total:
            for attempt in range(max_retries):
                try:
                    r = requests.get(next_url, timeout=timeout)
                    r.raise_for_status()
                    break
                except requests.RequestException as exc:
                    if attempt == max_retries - 1:
                        raise
                    wait = 2 ** attempt
                    print(f"\n  [retry {attempt + 1}/{max_retries} after {wait}s: {exc}]",
                          flush=True)
                    time.sleep(wait)
            text = r.text
            fh.write(text)
            collected += text.count("\n>") + (1 if text.startswith(">") else 0)
            link = r.headers.get("Link", "")
            m = _re.search(r'<([^>]+)>;\s*rel="next"', link)
            next_url = m.group(1) if m else None
            print(f"  collected ~{collected:,} sequences …", end="\r", flush=True)
    print()


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
    Sequences containing non-standard amino acids (B, O, U, Z, hyphens, etc.)
    are silently skipped — they cause BLASTP/MMseqs2 errors.
    """
    output_path = Path(output_path)
    seen: set[str] = set()
    written = 0
    skipped = 0

    with open(output_path, "w") as out:
        for db_name, path in db_paths.items():
            if not path.exists():
                continue
            with open(path, encoding="utf-8", errors="replace") as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    try:
                        seq = str(rec.seq).upper()
                    except UnicodeDecodeError:
                        seq = rec.seq._data.decode("ascii", errors="ignore").upper()
                    seq = seq.rstrip("*")
                    if not seq or not all(aa in _STANDARD_AAS for aa in seq):
                        skipped += 1
                        continue
                    rec.seq = Seq(seq)
                    key = hashlib.md5(seq.encode()).hexdigest()
                    if deduplicate and key in seen:
                        continue
                    seen.add(key)
                    rec.id = f"{db_name}|{rec.id}"
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    written += 1

    if skipped:
        print(f"[merge] {skipped:,} sequences skipped (non-standard amino acids)")
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
