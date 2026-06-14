"""
Generate ESM-2 per-sequence embeddings via mean pooling over residues.

Model: facebook/esm2_t12_35M_UR50D  (480-dim, frozen during inference)
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmModel

_MODEL_ID = "facebook/esm2_t12_35M_UR50D"
_EMBED_DIM = 480


def _sequence_hash(sequences: list[str]) -> str:
    """MD5 of the newline-joined sequence list — used to detect stale caches."""
    return hashlib.md5("\n".join(sequences).encode()).hexdigest()


def _resolve_device(device: str) -> str:
    if device == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return device


def _load_model(device: str) -> tuple:
    device = _resolve_device(device)
    print(f"[embed] Loading {_MODEL_ID} ...")
    tokenizer = AutoTokenizer.from_pretrained(_MODEL_ID)
    model = EsmModel.from_pretrained(_MODEL_ID)
    model.eval()
    model.to(device)
    return tokenizer, model


@torch.no_grad()
def embed_sequences(
    sequences: list[str],
    batch_size: int = 64,
    device: str | None = None,
    cache_path: Path | None = None,
) -> np.ndarray:
    """
    Embed *sequences* with frozen ESM-2 and return an (N, 480) float32 array.

    Mean-pools over the residue dimension (excluding special tokens).
    If *cache_path* is provided and the file exists, loads from cache instead
    of re-running the model.
    """
    if not sequences:
        return np.empty((0, _EMBED_DIM), dtype=np.float32)

    if cache_path and Path(cache_path).exists():
        hash_path = Path(str(cache_path) + ".hash")
        current_hash = _sequence_hash(sequences)
        if hash_path.exists() and hash_path.read_text().strip() == current_hash:
            print(f"[embed] Loading cached embeddings from {cache_path}")
            return np.load(cache_path)
        elif hash_path.exists():
            print(f"[embed] Cache hash mismatch — dataset changed, re-embedding …")
        else:
            print(f"[embed] No hash sidecar found — re-embedding to verify freshness …")

    device = _resolve_device(device or "auto")
    print(f"[embed] Device: {device}  |  Sequences: {len(sequences):,}")

    tokenizer, model = _load_model(device)
    all_embeddings: list[np.ndarray] = []

    for start in range(0, len(sequences), batch_size):
        batch = sequences[start : start + batch_size]
        inputs = tokenizer(
            batch,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=202,  # 200aa + 2 special tokens
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}

        outputs = model(**inputs)
        hidden = outputs.last_hidden_state  # (B, L, 480)

        # Mean-pool over residue positions, excluding padding and special tokens
        attention_mask = inputs["attention_mask"]  # (B, L)
        # Zero out special tokens (first and last non-pad position)
        mask = attention_mask.clone().float()
        mask[:, 0] = 0  # [CLS]
        for i, seq in enumerate(batch):
            # find last non-pad token = len(seq) + 1 (for [CLS])
            mask[i, len(seq) + 1 :] = 0  # [EOS] and padding

        mask_expanded = mask.unsqueeze(-1)  # (B, L, 1)
        summed = (hidden * mask_expanded).sum(dim=1)  # (B, 480)
        counts = mask.sum(dim=1, keepdim=True).clamp(min=1e-9)  # (B, 1)
        pooled = (summed / counts).cpu().numpy()  # (B, 480)
        all_embeddings.append(pooled)

        if (start // batch_size) % 10 == 0:
            pct = min(start + batch_size, len(sequences)) / len(sequences) * 100
            print(f"  {pct:.0f}%", end="\r", flush=True)

    print()
    embeddings = np.vstack(all_embeddings).astype(np.float32)

    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, embeddings)
        Path(str(cache_path) + ".hash").write_text(_sequence_hash(sequences))
        print(f"[embed] Saved {embeddings.shape} → {cache_path}")

    return embeddings


def embed_dataframe(
    df: pd.DataFrame,
    output_dir: Path,
    batch_size: int = 64,
    device: str | None = None,
    split: str | None = None,
) -> np.ndarray:
    """
    Embed all sequences in *df* and save to *output_dir*.

    If *split* is given (e.g. 'train', 'test'), only sequences with that
    split value are embedded and the cache file is named accordingly.
    """
    if split:
        subset = df[df["split"] == split].reset_index(drop=True)
        tag = split
    else:
        subset = df
        tag = "all"

    cache = Path(output_dir) / f"embeddings_{tag}.npy"
    labels_path = Path(output_dir) / f"labels_{tag}.npy"

    embeddings = embed_sequences(
        list(subset["sequence"]),
        batch_size=batch_size,
        device=device,
        cache_path=cache,
    )

    np.save(labels_path, subset["label"].values.astype(np.int8))
    print(f"[embed] Labels → {labels_path}")

    return embeddings
