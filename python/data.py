"""Experimental event loading and preprocessing for the blind pipeline.

Per subfile the events tree's raw ADC is calibrated per channel (OffsetLeft /
OffsetRight subtracted from an end that fired, then GainLeft / GainRight /
GainStrip0 / GainStrip17, from the same file's one-row calibration tree) and
assembled into the per-strip totals exactly like EnergyView::Decode:
total[s] = left[s] + right[s] (both ends summed), the unsegmented strips 0/17
kept unless INCLUDE_UNSEGMENTED_STRIPS is off, split strips kept long-side-only
when IGNORE_SHORT_STRIPS. After gain, the per-strip two-point alignment
(StripFactor and StripOffset from the calibration tree) is applied:
total[s] = total[s] * strip_factor[s] + strip_offset[s], which puts the beam
peak at 1 and the two-particle pile-up peak at 2 on every strip, so downstream
code uses UNIT gains. NO event selection is applied here
-- the blind clustering pipeline's step 1 is the first filter.
"""

import contextlib
import os

import numpy as np
from scipy.signal import savgol_filter

import config


def _assemble_totals(left, right, strip_factor=None, strip_offset=None):
    """Per-strip totals view from (n, 18) left/right arrays, mirroring
    EnergyView::Decode exactly.

    total[s] = left[s] + right[s] -- both ends summed, the C++ default
    (Constants::IGNORE_SHORT_STRIPS = false). The unsegmented strips 0/17
    sit wholly in the left column (right is zero there). With
    config.IGNORE_SHORT_STRIPS the split strips 1-16 instead keep only their
    long end (L_odd / R_even). Columns 0 and 17 are dropped when
    INCLUDE_UNSEGMENTED_STRIPS is off.

    When strip_factor / strip_offset are provided (18-element arrays from
    the calibration tree's StripFactor / StripOffset branches), they are
    applied after gain and IGNORE_SHORT_STRIPS:
    total[s] = total[s] * strip_factor[s] + strip_offset[s], matching
    EnergyView::Decode (which carries the offset on the long end, so the
    strip total is the same).
    """
    total = left + right
    if config.IGNORE_SHORT_STRIPS:
        for s in range(1, 17):
            total[:, s] = left[:, s] if s % 2 == 1 else right[:, s]
    # Per-strip two-point alignment (matches EnergyView::Decode: factor on
    # every end, offset on the long end, after IGNORE_SHORT_STRIPS).
    if strip_factor is not None:
        total *= strip_factor[np.newaxis, :]
    if strip_offset is not None:
        total += strip_offset[np.newaxis, :]
    if not config.INCLUDE_UNSEGMENTED_STRIPS:
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


def _with_unsegmented(seg, strip0, strip17):
    """Widen a 16-column (strips 1-16, index strip-1) block to the 18-strip
    view: column 0 is strip 0, column 17 is strip 17. seg is (n, 16) or (16,);
    strip0/strip17 are the matching (n,) columns or scalars."""
    seg = np.asarray(seg, dtype=np.float32)
    strip0 = np.asarray(strip0, dtype=np.float32)
    strip17 = np.asarray(strip17, dtype=np.float32)
    if seg.ndim == 1:
        return np.concatenate([strip0.reshape(1), seg, strip17.reshape(1)])
    n = seg.shape[0]
    return np.concatenate([
        np.broadcast_to(strip0, (n, )).reshape(n, 1), seg,
        np.broadcast_to(strip17, (n, )).reshape(n, 1)
    ],
                          axis=1)


def _load_calibrated_lr(path, max_events_per_file=None):
    """Calibrated (left, right) arrays (n, 18) + StripFactor / StripOffset
    for one events file, per channel like EnergyView::Decode; strip_factor
    and strip_offset are (18,) float32 (identity 1.0 / 0.0 if a branch is
    absent).

    The events tree stores strips 1-16 as LeftdE[16]/RightdE[16] (index
    strip-1) and the unsegmented strips 0 and 17 as the scalars Strip0dE /
    Strip17dE; the calibration tree matches with GainLeft[16]/GainRight[16]
    and GainStrip0/GainStrip17. Here they are widened to the 18-strip view
    the rest of the pipeline works in, with strips 0/17 in the left column
    and zero on the right."""
    from analysis_utilities.io import load_leaf_array_data, load_tree_data
    cache = str(config.CACHE_DIR)
    # Silence the loaders' per-file "Loading cached ..." lines (no quiet flag
    # upstream): redirect their stdout to devnull.
    with open(os.devnull, "w") as _devnull, \
            contextlib.redirect_stdout(_devnull):
        ev = load_leaf_array_data(path,
                                  "events", ["LeftdE", "RightdE"],
                                  max_events=max_events_per_file,
                                  cache_dir=cache)
        # Scalar branches come through the scalar loader; both loaders walk
        # the entries the same way, so their rows line up.
        ev_scalar = load_tree_data(path,
                                   "events",
                                   max_events=max_events_per_file,
                                   cache_dir=cache)
        cal = load_leaf_array_data(path,
                                   "calibration",
                                   ["GainLeft", "GainRight", "OffsetLeft",
                                    "OffsetRight", "StripFactor",
                                    "StripOffset"],
                                   cache_dir=cache)
        cal_scalar = load_tree_data(path, "calibration", cache_dir=cache)
    n = ev["LeftdE"].shape[0]
    if len(ev_scalar) != n:
        raise RuntimeError(
            f"{path}: {n} array rows but {len(ev_scalar)} scalar rows")
    raw_left = _with_unsegmented(ev["LeftdE"],
                                 ev_scalar["Strip0dE"].to_numpy(),
                                 ev_scalar["Strip17dE"].to_numpy())
    raw_right = _with_unsegmented(ev["RightdE"], 0.0, 0.0)
    gain_left = _with_unsegmented(cal["GainLeft"][0],
                                  cal_scalar["GainStrip0"].iloc[0],
                                  cal_scalar["GainStrip17"].iloc[0])
    gain_right = _with_unsegmented(cal["GainRight"][0], 0.0, 0.0)
    # Per-end offsets (short ends: the fixed amount a fired end reads over
    # its charge), subtracted only where the end fired and clamped at 0,
    # like EnergyView::Decode; old files lack them -> 0.
    if "OffsetLeft" in cal:
        off_left = _with_unsegmented(cal["OffsetLeft"][0], 0.0, 0.0)
        off_right = _with_unsegmented(cal["OffsetRight"][0], 0.0, 0.0)
    else:
        off_left = np.zeros(18, dtype=np.float32)
        off_right = np.zeros(18, dtype=np.float32)
    fired_l = raw_left > 0
    fired_r = raw_right > 0
    left = np.where(fired_l,
                    np.maximum(0.0, raw_left - off_left[np.newaxis, :]),
                    0.0) * gain_left[np.newaxis, :]
    right = np.where(fired_r,
                     np.maximum(0.0, raw_right - off_right[np.newaxis, :]),
                     0.0) * gain_right[np.newaxis, :]
    # StripFactor / StripOffset: per-strip two-point alignment (beam peak to
    # 1, pile-up peak to 2); old files lack them -> identity.
    if "StripFactor" in cal:
        strip_factor = cal["StripFactor"][0].astype(np.float32)
    else:
        strip_factor = np.ones(18, dtype=np.float32)
    if "StripOffset" in cal:
        strip_offset = cal["StripOffset"][0].astype(np.float32)
    else:
        strip_offset = np.zeros(18, dtype=np.float32)
    # Both-channel firing from the RAW ADC, not the calibrated ends: short-end
    # gains are 0 (no sim anchor), so calibrated short ends are always zero.
    both = _both_fired(raw_left, raw_right, config.BLIND_MULT_THRESH)
    return left, right, both, strip_factor, strip_offset


def _both_fired(left, right, thresh):
    """Per-strip both-channel fired mask (n, 18): True where BOTH raw ADC ends
    of a split strip (1-16) exceed thresh -- firing, not energy. The
    unsegmented strips 0/17 have one end, so always False. The per-event
    multiplicity is the row sum over strips 1-16; the full per-strip mask is
    kept so a firing can be localized to the reaction strip."""
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
        left, right, both, strip_factor, strip_offset = _load_calibrated_lr(
            path, max_events_per_file)
        print(f"  {os.path.basename(path)}: {left.shape[0]} events "
              "(zero cuts)")
        out.append(_assemble_totals(left, right, strip_factor, strip_offset))
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
