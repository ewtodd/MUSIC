"""Experimental event loading and preprocessing for the blind pipeline.

Per subfile the events tree's raw ADC arrays are calibrated per channel
(GainLeft/GainRight from the same file's one-row calibration tree) and
assembled into the per-strip totals exactly like EnergyView::Decode:
total[s] = left[s] + right[s] (both ends summed), guard strips 0/17 kept
unless INCLUDE_GUARD_STRIPS is off, split strips kept long-side-only when
IGNORE_SHORT_STRIPS. The calibration already puts the beam at ~1 per strip,
so downstream code uses UNIT gains. NO event selection is applied here --
the blind clustering pipeline's step 1 is the first filter.
"""

import contextlib
import os

import numpy as np
from scipy.signal import savgol_filter

import config


def _assemble_totals(left, right):
    """Per-strip totals view from (n, 18) left/right arrays, mirroring
    EnergyView::Decode exactly.

    total[s] = left[s] + right[s] -- both ends summed, the C++ default
    (Constants::IGNORE_SHORT_STRIPS = false). Guard strips 0/17 are
    single-ended (one of left/right is zero). With config.IGNORE_SHORT_STRIPS
    the split strips 1-16 instead keep only their long end (L_odd / R_even).
    Guard columns are dropped when INCLUDE_GUARD_STRIPS is off.
    """
    total = left + right
    if config.IGNORE_SHORT_STRIPS:
        for s in range(1, 17):
            total[:, s] = left[:, s] if s % 2 == 1 else right[:, s]
    if not config.INCLUDE_GUARD_STRIPS:
        total = total[:, 1:17]
    return total


def list_event_files():
    """Sorted experimental Events_Run*.root paths under root_files."""
    base = config.ROOT_FILES_DIR
    if not base.is_dir():
        raise FileNotFoundError(f"no root_files dir: {base}")
    files = sorted(
        str(base / entry) for entry in os.listdir(base)
        if entry.startswith("Events_Run") and entry.endswith(".root"))
    if not files:
        raise FileNotFoundError(f"no Events_Run*.root files in {base}")
    return files


def _load_calibrated_lr(path, max_events_per_file=None):
    """Calibrated (left, right) arrays for one events file, per channel
    like EnergyView::Decode (the per-file one-row calibration tree)."""
    from analysis_utilities.io import load_leaf_array_data
    cache = str(config.CACHE_DIR)
    # Silence load_leaf_array_data's per-file "Loading cached leaf arrays"
    # line (no quiet flag upstream): redirect its stdout to devnull.
    with open(os.devnull, "w") as _devnull, \
            contextlib.redirect_stdout(_devnull):
        ev = load_leaf_array_data(path,
                                  "events", ["Left_0_17_dE", "RightdE"],
                                  max_events=max_events_per_file,
                                  cache_dir=cache)
        cal = load_leaf_array_data(path,
                                   "calibration", ["GainLeft", "GainRight"],
                                   cache_dir=cache)
    left = ev["Left_0_17_dE"].astype(np.float32) * \
        cal["GainLeft"][0].astype(np.float32)
    right = ev["RightdE"].astype(np.float32) * \
        cal["GainRight"][0].astype(np.float32)
    # Both-channel firing from the RAW ADC, not the calibrated ends: the
    # short-end gains are 0 (uncalibrated -- no sim anchor), so calibrated
    # short ends are always zero. We only care whether a channel fired, and
    # the raw ADC carries that regardless of gain.
    both = _both_fired(ev["Left_0_17_dE"], ev["RightdE"],
                       config.BLIND_MULT_THRESH)
    return left, right, both


def _both_fired(left, right, thresh):
    """Per-strip both-channel fired mask (n, 18): True where BOTH raw ADC ends
    of a split strip (1-16) exceed thresh -- firing, not energy. Guard strips
    0/17 are single-ended, so always False. The per-event multiplicity is the
    row sum over strips 1-16; the full per-strip mask is kept so a firing can
    be localized to the reaction strip."""
    fired = (left.astype(np.float64) > thresh) & \
        (right.astype(np.float64) > thresh)
    fired[:, 0] = False
    fired[:, 17] = False
    return fired


def load_experimental_totals(max_files=None, max_events_per_file=None):
    """Calibrated per-strip totals from the experimental events files, with
    ZERO selection -- every event is returned (the blind pipeline's own
    first stage is the filter). Returns (totals, both): totals is the summed
    per-strip view (already in a.u., beam ~ 1, so downstream uses UNIT
    gains); both is the per-strip both-channel fired mask (n, 18)."""
    files = list_event_files()
    if max_files is not None:
        files = files[:max_files]
    out = []
    both_out = []
    for path in files:
        left, right, both = _load_calibrated_lr(path, max_events_per_file)
        print(f"  {os.path.basename(path)}: {left.shape[0]} events "
              "(zero cuts)")
        out.append(_assemble_totals(left, right))
        both_out.append(both)
    if not out:
        raise RuntimeError("no experimental events loaded")
    return np.concatenate(out), np.concatenate(both_out)


def max_normalize(totals, gain):
    """Flatten per strip by gain, then drop rows whose flattened max is <= 0.

    Returns (X, kept): X is float32 (n_kept, n_cols); kept is the boolean
    row mask. Flattened a.u. (beam == 1), absolute amplitude preserved (no
    per-trace peak division)."""
    flat = totals.astype(np.float64) * gain[np.newaxis, :]
    peak = flat.max(axis=1)
    kept = peak > 0.0
    return flat[kept].astype(np.float32), kept


def derivative(X):
    """Adjacent-strip differences of the trace, same scale as the trace."""
    return np.diff(X, axis=1).astype(np.float32)


def savgol_filter_trace(X):
    """Savitzky-Golay smooth each row of X (n, 18), matching
    StripSumScatter::SavitzkyGolay in C++: 5-point window, cubic polynomial,
    coefficients [-3, 12, 17, 12, -3] / 35. At edges (strips 0, 1, 16, 17)
    the window is clipped and coefficients are renormalized by the sum of
    included coefficients — identical to the C++ val / wsum path.

    Uses scipy.signal.savgol_filter to derive the coefficients (rather than
    hard-coding them); the convolution itself is manual so that the edge
    handling exactly matches the C++ clip-and-renormalize rather than scipy's
    default interpolation mode.

    Returns float64 array of the same shape as X.
    """
    Xd = np.asarray(X, dtype=np.float64)
    n, m = Xd.shape
    out = np.empty((n, m), dtype=np.float64)
    # Derive SG coefficients from scipy (5-point, cubic).
    c = savgol_filter(np.array([0, 0, 1, 0, 0], dtype=np.float64),
                      window_length=5,
                      polyorder=3)
    K = 2  # half-window
    for s in range(m):
        lo = max(0, s - K)
        hi = min(m - 1, s + K)
        offset = lo - s + K  # start index in c
        seg = c[offset:offset + (hi - lo + 1)]
        inv_w = 1.0 / seg.sum()
        out[:, s] = (Xd[:, lo:hi + 1] * seg).sum(axis=1) * inv_w
    return out
