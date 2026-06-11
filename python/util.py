"""Shared training/evaluation helpers."""

import random

import numpy as np
import torch

import config


def set_seed(seed=config.SEED):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def get_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def split_indices(n, val_fraction, rng):
    """Shuffled (train_idx, val_idx) split of range(n)."""
    idx = rng.permutation(n)
    n_val = int(round(n * val_fraction))
    return idx[n_val:], idx[:n_val]


_LOG_EPS = 1e-12


def cycle_batches(loader):
    """Endless batch iterator (the negatives pool is smaller/larger than
    the beam pool, so the two loaders never align per epoch)."""
    while True:
        for batch in loader:
            yield batch


def per_trace_mse(model, xb):
    return torch.mean((model(xb) - xb)**2, dim=1)


def oe_hinge(err_neg, beam_level, margin_ratio):
    """Mean outlier-exposure hinge over a batch of negative per-trace MSEs.

    max(0, log(R * beam_level) - log(err_neg)) per trace: zero once the
    negative reconstructs at least margin_ratio times worse than the
    (detached) beam reconstruction level, log-space so the pull is
    scale-free and bounded.
    """
    floor = torch.log(margin_ratio * beam_level + _LOG_EPS)
    return torch.clamp(floor - torch.log(err_neg + _LOG_EPS), min=0.0).mean()


def roc_points(pos_scores, neg_scores):
    """(fpr, tpr) arrays, descending-threshold sweep (higher score = signal)."""
    scores = np.concatenate([pos_scores, neg_scores])
    labels = np.concatenate(
        [np.ones(pos_scores.size),
         np.zeros(neg_scores.size)])
    order = np.argsort(-scores, kind="stable")
    tp = np.cumsum(labels[order])
    fp = np.cumsum(1.0 - labels[order])
    tpr = tp / max(pos_scores.size, 1)
    fpr = fp / max(neg_scores.size, 1)
    return fpr, tpr


def auc(pos_scores, neg_scores):
    """ROC AUC via the Mann-Whitney rank statistic (ties get half credit)."""
    scores = np.concatenate([pos_scores, neg_scores])
    order = scores.argsort(kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, scores.size + 1)
    # average ranks over ties
    sorted_scores = scores[order]
    i = 0
    while i < sorted_scores.size:
        j = i
        while (j + 1 < sorted_scores.size
               and sorted_scores[j + 1] == sorted_scores[i]):
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = 0.5 * (i + 1 + j + 1)
        i = j + 1
    n_pos = pos_scores.size
    n_neg = neg_scores.size
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    rank_sum = ranks[:n_pos].sum()
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
