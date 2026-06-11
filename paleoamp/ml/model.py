"""
MLP classifier trained on ESM-2 embeddings to predict AMP probability.

Architecture: 480 → 256 → 64 → 1
Loss:         BCEWithLogitsLoss(pos_weight=3.0)
Optimizer:    AdamW, lr=1e-3, weight_decay=1e-4
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


class AMPClassifier(nn.Module):
    def __init__(self, input_dim: int = 480, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


@dataclass
class TrainConfig:
    input_dim: int = 480
    dropout: float = 0.3
    lr: float = 1e-3
    weight_decay: float = 1e-4
    pos_weight: float = 3.0       # 1:3 pos:neg ratio per AmPEP recommendation
    batch_size: int = 256
    max_epochs: int = 50
    patience: int = 7             # early stopping patience (epochs without AUPR gain)
    device: str = "auto"


@dataclass
class TrainResult:
    best_epoch: int
    best_val_aupr: float
    train_loss_history: list[float] = field(default_factory=list)
    val_aupr_history: list[float] = field(default_factory=list)


def _aupr(labels: np.ndarray, scores: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(labels, scores))


def train(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
    config: TrainConfig,
    checkpoint_dir: Path,
) -> tuple[AMPClassifier, TrainResult]:
    """
    Train the MLP classifier and save the best checkpoint by val AUPR.

    Returns the best model and training statistics.
    """
    device = (
        torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if config.device == "auto"
        else torch.device(config.device)
    )
    print(f"[train] Device: {device}")

    model = AMPClassifier(config.input_dim, config.dropout).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.lr, weight_decay=config.weight_decay
    )
    loss_fn = nn.BCEWithLogitsLoss(
        pos_weight=torch.tensor([config.pos_weight], device=device)
    )

    train_ds = TensorDataset(
        torch.from_numpy(X_train).float(),
        torch.from_numpy(y_train).float(),
    )
    train_loader = DataLoader(train_ds, batch_size=config.batch_size, shuffle=True)

    val_X = torch.from_numpy(X_val).float().to(device)
    val_y = y_val

    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    best_ckpt = checkpoint_dir / "best_model.pt"

    result = TrainResult(best_epoch=0, best_val_aupr=0.0)
    no_improve = 0

    for epoch in range(1, config.max_epochs + 1):
        model.train()
        epoch_loss = 0.0
        for X_b, y_b in train_loader:
            X_b, y_b = X_b.to(device), y_b.to(device)
            optimizer.zero_grad()
            logits = model(X_b)
            loss = loss_fn(logits, y_b)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(X_b)

        epoch_loss /= len(train_ds)
        result.train_loss_history.append(epoch_loss)

        # Validation
        model.eval()
        with torch.no_grad():
            val_logits = model(val_X).cpu().numpy()
        val_scores = 1 / (1 + np.exp(-val_logits))  # sigmoid
        val_aupr = _aupr(val_y, val_scores)
        result.val_aupr_history.append(val_aupr)

        print(
            f"  Epoch {epoch:3d}/{config.max_epochs}  "
            f"loss={epoch_loss:.4f}  val_AUPR={val_aupr:.4f}",
            end="",
        )

        if val_aupr > result.best_val_aupr:
            result.best_val_aupr = val_aupr
            result.best_epoch = epoch
            torch.save(model.state_dict(), best_ckpt)
            no_improve = 0
            print("  ✓ best")
        else:
            no_improve += 1
            print()
            if no_improve >= config.patience:
                print(f"[train] Early stopping at epoch {epoch}")
                break

    print(
        f"\n[train] Best epoch: {result.best_epoch}  "
        f"val_AUPR: {result.best_val_aupr:.4f}  "
        f"checkpoint: {best_ckpt}"
    )
    model.load_state_dict(torch.load(best_ckpt, map_location=device))
    return model, result


def load_model(checkpoint: Path, input_dim: int = 480, device: str = "cpu") -> AMPClassifier:
    model = AMPClassifier(input_dim=input_dim)
    model.load_state_dict(torch.load(checkpoint, map_location=device))
    model.eval()
    return model


def predict_proba(
    model: AMPClassifier,
    embeddings: np.ndarray,
    batch_size: int = 256,
    device: str = "cpu",
) -> np.ndarray:
    """Return AMP probability scores (0–1) for each embedding."""
    model = model.to(device)
    model.eval()
    all_scores: list[np.ndarray] = []

    with torch.no_grad():
        for start in range(0, len(embeddings), batch_size):
            batch = torch.from_numpy(embeddings[start : start + batch_size]).float().to(device)
            logits = model(batch).cpu().numpy()
            scores = 1 / (1 + np.exp(-logits))
            all_scores.append(scores)

    return np.concatenate(all_scores)
