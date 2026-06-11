"""
Multi-criteria validation of novel AMP candidates to reduce false positives.

Three independent validation layers:
  1. Physicochemical filters — net charge, hydrophobic moment, length,
                               low-complexity, duplicate removal
  2. Prodigal GFF partial-ORF check — flags incomplete gene calls
  3. AAC logistic regression — second classifier trained on amino acid
                               composition, fully independent of ESM-2
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ---------------------------------------------------------------------------
# Scales
# ---------------------------------------------------------------------------

# Eisenberg consensus hydrophobicity (for hydrophobic moment)
_EISENBERG: dict[str, float] = {
    "A":  0.25, "R": -1.80, "N": -0.64, "D": -0.72, "C":  0.04,
    "Q": -0.69, "E": -0.62, "G":  0.16, "H": -0.40, "I":  0.73,
    "L":  0.53, "K": -1.10, "M":  0.26, "F":  0.61, "P": -0.07,
    "S": -0.26, "T": -0.18, "W":  0.37, "Y":  0.02, "V":  0.54,
}

# Boman index (protein-binding potential)
_BOMAN: dict[str, float] = {
    "L": -1.70, "A": -0.61, "V": -1.33, "I": -1.56, "P":  0.00,
    "W": -0.64, "F": -1.16, "M": -0.64, "G":  0.48, "T":  0.45,
    "S":  0.45, "C":  0.51, "Y":  0.12, "H":  0.81, "Q":  0.95,
    "N":  0.95, "K":  1.54, "R":  1.54, "D":  0.62, "E":  0.62,
}

# Net charge pKa contributions (simplified Henderson-Hasselbalch at pH 7)
_CHARGE_PH7: dict[str, float] = {
    "K":  1.0, "R":  1.0, "H":  0.1,
    "D": -1.0, "E": -1.0,
}
_NTERM_CHARGE =  1.0
_CTERM_CHARGE = -1.0

_STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# ---------------------------------------------------------------------------
# Physicochemical calculations
# ---------------------------------------------------------------------------

def net_charge(seq: str) -> float:
    seq = seq.upper()
    charge = _NTERM_CHARGE + _CTERM_CHARGE
    for aa in seq:
        charge += _CHARGE_PH7.get(aa, 0.0)
    return charge


def hydrophobic_moment(seq: str, angle: float = 100.0) -> float:
    """
    Eisenberg hydrophobic moment for a helix with *angle* degrees per residue.
    Higher values indicate greater amphipathicity.
    """
    seq = seq.upper()
    rad = math.radians(angle)
    sin_sum = sum(_EISENBERG.get(aa, 0.0) * math.sin(i * rad) for i, aa in enumerate(seq))
    cos_sum = sum(_EISENBERG.get(aa, 0.0) * math.cos(i * rad) for i, aa in enumerate(seq))
    return math.sqrt(sin_sum ** 2 + cos_sum ** 2) / len(seq)


def boman_index(seq: str) -> float:
    """Mean Boman index — values > 2.48 suggest protein-binding potential."""
    seq = seq.upper()
    vals = [_BOMAN.get(aa, 0.0) for aa in seq]
    return sum(vals) / len(vals) if vals else 0.0


def shannon_entropy(seq: str) -> float:
    """Per-residue Shannon entropy in bits. Low values indicate low complexity."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    counts = {aa: seq.count(aa) for aa in set(seq)}
    return -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)


def instability_index(seq: str) -> float:
    """Biopython ProteinAnalysis instability index. < 40 = stable."""
    seq = seq.upper().rstrip("*")
    clean = "".join(aa for aa in seq if aa in _STANDARD_AAS)
    if len(clean) < 2:
        return 999.0
    try:
        return ProteinAnalysis(clean).instability_index()
    except Exception:
        return 999.0


def compute_properties(seq: str) -> dict:
    seq = seq.upper().rstrip("*")
    return {
        "length":            len(seq),
        "net_charge":        net_charge(seq),
        "hydrophobic_moment": hydrophobic_moment(seq),
        "boman_index":       boman_index(seq),
        "instability_index": instability_index(seq),
        "shannon_entropy":   shannon_entropy(seq),
        "is_low_complexity": shannon_entropy(seq) < 1.5,
        "has_repeat":        _has_repeat(seq),
        "starts_with_met":   seq.startswith("M"),
        "frac_standard_aa":  sum(1 for aa in seq if aa in _STANDARD_AAS) / max(len(seq), 1),
    }


def _has_repeat(seq: str, min_unit: int = 3, min_copies: int = 3) -> bool:
    """Return True if any short motif repeats ≥ min_copies times consecutively."""
    for unit_len in range(min_unit, min(8, len(seq) // min_copies + 1)):
        for start in range(len(seq) - unit_len * min_copies + 1):
            unit = seq[start: start + unit_len]
            pattern = f"({re.escape(unit)}){{{min_copies},}}"
            if re.search(pattern, seq[start:]):
                return True
    return False


# ---------------------------------------------------------------------------
# Prodigal GFF partial flag
# ---------------------------------------------------------------------------

def parse_partial_flags(gff_path: Path) -> dict[str, str]:
    """
    Parse Prodigal GFF and return {orf_id: partial_flag} where partial_flag
    is one of '00' (complete), '10' (no start), '01' (no stop), '11' (both).
    """
    flags: dict[str, str] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            m_id = re.search(r"ID=([^;]+)", line)
            m_partial = re.search(r"partial=(\d{2})", line)
            m_contig = re.match(r"(\S+)\t", line)
            if m_id and m_partial and m_contig:
                contig = m_contig.group(1)
                orf_num = m_id.group(1).split("_")[-1]
                orf_id = f"{contig}_{orf_num}"
                flags[orf_id] = m_partial.group(1)
    return flags


# ---------------------------------------------------------------------------
# AAC second classifier (independent of ESM-2)
# ---------------------------------------------------------------------------

def _aac_features(seq: str) -> np.ndarray:
    """20-dim amino acid composition vector (fractions, sum to 1)."""
    seq = seq.upper()
    n = len(seq)
    return np.array([seq.count(aa) / n for aa in _AA_ORDER], dtype=np.float32)


def train_aac_classifier(dataset_csv: Path):
    """
    Train a logistic regression on amino acid composition using the same
    dataset used for ESM-2 training.  Returns a fitted sklearn estimator.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline

    df = pd.read_csv(dataset_csv)
    df = df[df["sequence"].apply(lambda s: all(aa in _STANDARD_AAS for aa in s.upper()))]

    X = np.vstack([_aac_features(s) for s in df["sequence"]])
    y = df["label"].values

    clf = Pipeline([
        ("scaler", StandardScaler()),
        ("lr", LogisticRegression(
            C=1.0, max_iter=1000, class_weight="balanced", random_state=42
        )),
    ])
    clf.fit(X, y)
    return clf


def aac_predict_proba(clf, sequences: list[str]) -> np.ndarray:
    X = np.vstack([_aac_features(s) for s in sequences])
    return clf.predict_proba(X)[:, 1]


# ---------------------------------------------------------------------------
# Main validation entry point
# ---------------------------------------------------------------------------

_PHYSCOCHEM_CRITERIA = {
    "length_ok":          lambda p: 10 <= p["length"] <= 100,
    "charge_ok":          lambda p: p["net_charge"] >= 0.0,
    "amphipathic_ok":     lambda p: p["hydrophobic_moment"] >= 0.10,
    "not_low_complexity": lambda p: not p["is_low_complexity"],
    "no_repeat":          lambda p: not p["has_repeat"],
    "standard_aa_ok":     lambda p: p["frac_standard_aa"] == 1.0,
}


def validate_candidates(
    novel_tsv: Path,
    gff_path: Path,
    dataset_csv: Path,
    output_dir: Path,
    aac_threshold: float = 0.4,
) -> pd.DataFrame:
    """
    Run all three validation layers and return a ranked DataFrame with
    per-criterion flags and a final PASS/WARN/FAIL verdict.

    aac_threshold — minimum AAC classifier probability to count as a vote
                    (set lower than 0.5 since AAC is weaker than ESM-2)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(novel_tsv, sep="\t")

    # --- Layer 1: duplicates by exact sequence ---
    df["is_duplicate"] = df["sequence"].duplicated(keep="first")

    # --- Layer 1: physicochemical ---
    props = [compute_properties(s) for s in df["sequence"]]
    for key in _PHYSCOCHEM_CRITERIA:
        df[key] = [_PHYSCOCHEM_CRITERIA[key](p) for p in props]
    for prop in ("length", "net_charge", "hydrophobic_moment",
                 "boman_index", "instability_index", "shannon_entropy",
                 "starts_with_met"):
        df[prop] = [p[prop] for p in props]
    df["physcochem_pass"] = df[[k for k in _PHYSCOCHEM_CRITERIA]].all(axis=1)

    # --- Layer 2: Prodigal partial flag ---
    partial_flags = parse_partial_flags(gff_path)
    df["partial_flag"] = df["seq_id"].map(partial_flags).fillna("??")
    df["is_complete_orf"] = df["partial_flag"] == "00"

    # --- Layer 3: AAC second classifier ---
    print("[validate] Training AAC classifier …")
    clf = train_aac_classifier(dataset_csv)
    df["aac_score"] = aac_predict_proba(clf, list(df["sequence"]))
    df["aac_pass"] = df["aac_score"] >= aac_threshold

    # --- Consensus verdict ---
    n_pass_criteria = (
        df["physcochem_pass"].astype(int)
        + df["aac_pass"].astype(int)
        + df["is_complete_orf"].astype(int)
        + (~df["is_duplicate"]).astype(int)
    )

    def _verdict(is_dup, n):
        if is_dup:
            return "DUPLICATE"
        if n == 4:
            return "PASS"
        if n >= 2:
            return "WARN"
        return "FAIL"

    df["verdict"] = [_verdict(dup, n) for dup, n in zip(df["is_duplicate"], n_pass_criteria)]
    df["n_criteria_met"] = n_pass_criteria

    df = df.sort_values(["n_criteria_met", "amp_score"], ascending=[False, False])

    out_path = output_dir / f"{Path(novel_tsv).stem}_validated.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    return df
