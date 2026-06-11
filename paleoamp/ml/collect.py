"""
Collect, filter, and deduplicate AMP training data.

Positives: APD6, DRAMP 3.0, DBAASP
Negatives: AMPlify negative set + UniProt non-AMP peptides
"""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

_STANDARD_AAS = frozenset("ACDEFGHIKLMNPQRSTVWY")

_SOURCES: dict[str, dict] = {
    # Positives
    "apd": {
        # APD2024a — all natural AMPs; verify latest filename at aps.unmc.edu/downloads
        "url": "https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta",
        "filename": "apd_natural.fasta",
        "label": 1,
        "description": "APD2024a — natural antimicrobial peptides",
    },
    "dramp_general": {
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta",
        "filename": "dramp_general.fasta",
        "label": 1,
        "description": "DRAMP 3.0 — general AMPs (experimentally validated)",
    },
    "dramp_specific": {
        "url": "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/specific_amps.fasta",
        "filename": "dramp_specific.fasta",
        "label": 1,
        "description": "DRAMP 3.0 — specific AMPs (target-activity annotations)",
    },
    "uniprot_amp": {
        # Reviewed SwissProt proteins with AMP keyword (KW-0929), length 5–200 aa.
        # Replaces DBAASP which no longer offers a bulk download endpoint.
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=reviewed:true+keyword:KW-0929+length:[5+TO+200]"
            "&format=fasta&size=500"
        ),
        "filename": "uniprot_amp.fasta",
        "label": 1,
        "description": "UniProt reviewed AMP peptides (KW-0929)",
        "paginate": 10000,
    },
    # Negatives
    "amplify_neg": {
        # AMPlify negative set — check https://github.com/bcgsc/AMPlify for current path
        "url": (
            "https://raw.githubusercontent.com/bcgsc/AMPlify/main/"
            "src/AMPlify/models/train_val_test_split/AMPlify_negative_dataset.csv"
        ),
        "filename": "amplify_negatives.csv",
        "label": 0,
        "description": "AMPlify negative training set",
        "format": "csv",
        "seq_col": "sequence",
    },
    "uniprot_non_amp": {
        # Reviewed SwissProt proteins without AMP keyword (KW-0929), length 5–200 aa.
        # Fetches up to 10,000 to give a ~3:1 neg:pos ratio vs APD.
        "url": (
            "https://rest.uniprot.org/uniprotkb/search"
            "?query=reviewed:true+NOT+keyword:KW-0929+length:[5+TO+200]"
            "&format=fasta&size=500"
        ),
        "filename": "uniprot_non_amp.fasta",
        "label": 0,
        "description": "UniProt reviewed non-AMP peptides",
        "paginate": 50000,  # total sequences to collect across pages
    },
}


def _download_paginated(url: str, dest: Path, total: int) -> None:
    """Fetch a UniProt FASTA result set across multiple pages until *total* seqs collected."""
    collected = 0
    next_url: str | None = url
    with open(dest, "w") as fh:
        while next_url and collected < total:
            r = requests.get(next_url, timeout=60)
            r.raise_for_status()
            text = r.text
            fh.write(text)
            collected += text.count("\n>") + (1 if text.startswith(">") else 0)
            # Parse Link header for next page
            link = r.headers.get("Link", "")
            import re as _re
            m = _re.search(r'<([^>]+)>;\s*rel="next"', link)
            next_url = m.group(1) if m else None
            print(f"  collected ~{collected:,} sequences …", end="\r", flush=True)
    print()


def download_sources(
    output_dir: Path,
    sources: list[str] | None = None,
    force: bool = False,
) -> dict[str, Path]:
    """Download raw source files. Returns {name: local_path} for successful downloads."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    targets = sources or list(_SOURCES.keys())
    paths: dict[str, Path] = {}

    for name in targets:
        info = _SOURCES[name]
        dest = output_dir / info["filename"]

        if dest.exists() and not force:
            print(f"[skip] {name}: already downloaded")
            paths[name] = dest
            continue

        print(f"[download] {name}: {info['description']}")
        try:
            if info.get("paginate"):
                _download_paginated(info["url"], dest, total=info["paginate"])
            else:
                r = requests.get(info["url"], stream=True, timeout=60)
                r.raise_for_status()
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(65536):
                        fh.write(chunk)
            print(f"[ok] {name} → {dest}")
            paths[name] = dest
        except Exception as exc:
            print(f"[error] {name}: {exc}", file=sys.stderr)

    return paths


def _is_standard(seq: str) -> bool:
    return bool(seq) and all(aa in _STANDARD_AAS for aa in seq.upper())


def _load_fasta(path: Path, label: int) -> list[tuple[str, str, int]]:
    """Return [(id, sequence, label)] from a FASTA file."""
    records = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            try:
                seq = str(rec.seq)
            except UnicodeDecodeError:
                # Some DRAMP sequences contain non-ASCII bytes (modified residues);
                # strip them — _is_standard will exclude the sequence downstream.
                seq = rec.seq._data.decode("ascii", errors="ignore")
            seq = seq.upper().rstrip("*")
            records.append((rec.id, seq, label))
    return records


def _load_csv(path: Path, seq_col: str, label: int) -> list[tuple[str, str, int]]:
    """Return [(id, sequence, label)] from a CSV with a sequence column."""
    df = pd.read_csv(path)
    if seq_col not in df.columns:
        raise ValueError(f"Column '{seq_col}' not found in {path}. Available: {list(df.columns)}")
    records = []
    for i, row in df.iterrows():
        seq = str(row[seq_col]).upper().strip().rstrip("*")
        records.append((f"neg_{i}", seq, label))
    return records


def build_dataset(
    raw_dir: Path,
    min_len: int = 5,
    max_len: int = 200,
) -> pd.DataFrame:
    """
    Load all downloaded source files, apply length + standard-AA filters,
    deduplicate by exact sequence, and return a DataFrame with columns:
      sequence, label, source, seq_id
    """
    all_records: list[tuple[str, str, int, str]] = []  # (id, seq, label, source)

    for name, info in _SOURCES.items():
        dest = raw_dir / info["filename"]
        if not dest.exists():
            print(f"[warn] {name}: file not found, skipping", file=sys.stderr)
            continue

        fmt = info.get("format", "fasta")
        label = info["label"]
        try:
            if fmt == "csv":
                rows = _load_csv(dest, info["seq_col"], label)
            else:
                rows = _load_fasta(dest, label)
            all_records.extend((sid, seq, lbl, name) for sid, seq, lbl in rows)
            print(f"[load] {name}: {len(rows):,} sequences")
        except Exception as exc:
            print(f"[error] loading {name}: {exc}", file=sys.stderr)

    df = pd.DataFrame(all_records, columns=["seq_id", "sequence", "label", "source"])

    # Filter
    before = len(df)
    df = df[df["sequence"].apply(_is_standard)]
    df = df[df["sequence"].str.len().between(min_len, max_len)]
    print(f"[filter] {len(df):,} / {before:,} sequences passed length+AA filter")

    # Deduplicate by exact sequence, keeping first occurrence
    before = len(df)
    df["_hash"] = df["sequence"].apply(lambda s: hashlib.md5(s.encode()).hexdigest())
    df = df.drop_duplicates(subset=["_hash"]).drop(columns=["_hash"])
    print(f"[dedup] {len(df):,} unique sequences ({before - len(df):,} duplicates removed)")

    df = df.reset_index(drop=True)
    print(f"\nPositives: {(df['label'] == 1).sum():,}   Negatives: {(df['label'] == 0).sum():,}")
    return df


def save_dataset(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"[save] {len(df):,} sequences → {path}")


def load_dataset(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)
