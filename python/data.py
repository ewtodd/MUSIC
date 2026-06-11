"""Sim trace loading and preprocessing.

Sim populations live in <dataset>/sim_root_files as traces_<ds>_<tag>.root
(tags like beam_eres, aa_eres_s4, an_eres_s7, mirroring RemixSim); the
events_MeV tree holds per-strip Left_0_17_dE[18] / RightdE[18] deposits in
MeV. Training uses the _eres (resolution-smeared) variants only.

Preprocessing chain (inference on data must apply the same normalization to
its already per-channel-calibrated totals):

  1. long-side trace (L_odd / R_even; guard-strip columns 0/17 dropped
     when INCLUDE_GUARD_STRIPS is off), plus, when INCLUDE_SHORT_STRIPS is
     on, the 16 short-side ends appended as a second column block (never
     summed -- see config.block_widths())
  2. per-strip flatten: total[s] *= gain[s], with gain[s] =
     FLAT_TARGET / <mean sim-beam total[s]>, mirroring
     StripSumScatter::SimBeamGains and the data's per-channel calibration
  4. with INCLUDE_DERIVATIVE, adjacent-strip differences of the trace
     (block-local), appended flat
"""

import array as _array
import os
import re
from dataclasses import dataclass

import numpy as np
import ROOT

import config

_STRIP_SUFFIX_RE = re.compile(r"^(.*)_s(\d+)$")


@dataclass
class SimSpec:
    tag: str  # full tag, e.g. "an_eres_s3"
    base: str  # population: "beam", "aa", "an", ...
    strip: int  # reaction strip; -1 for unreacted beam
    eres: bool  # resolution-smeared variant
    path: str


def list_sim_specs(eres_only=True):
    """Enumerate sim files by name, sorted by tag (RemixSim convention)."""
    prefix = f"traces_{config.DATASET}_"
    specs = []
    if not config.SIM_ROOT_DIR.is_dir():
        raise FileNotFoundError(
            f"no sim_root_files dir: {config.SIM_ROOT_DIR}")
    for entry in sorted(os.listdir(config.SIM_ROOT_DIR)):
        if not (entry.startswith(prefix) and entry.endswith(".root")):
            continue
        tag = entry[len(prefix):-len(".root")]
        m = _STRIP_SUFFIX_RE.match(tag)
        if m:
            base, strip = m.group(1), int(m.group(2))
        else:
            base, strip = tag, -1
        eres = base.endswith("_eres")
        if eres:
            base = base[:-len("_eres")]
        if eres_only and not eres:
            continue
        specs.append(
            SimSpec(tag=tag,
                    base=base,
                    strip=strip,
                    eres=eres,
                    path=str(config.SIM_ROOT_DIR / entry)))
    return specs


def load_left_right(path, tree_name="events_MeV"):
    """Per-strip left/right deposits (n_events, 18) each, npz-cached.

    Cached unsummed so the short-strip / guard-strip knobs never invalidate
    the cache; load_totals assembles the configured view on top.
    """
    config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    stem = os.path.splitext(os.path.basename(path))[0]
    cache = config.CACHE_DIR / f"{stem}_leftright.npz"
    if cache.exists() and cache.stat().st_mtime >= os.path.getmtime(path):
        loaded = np.load(cache)
        return loaded["left"], loaded["right"]

    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise IOError(f"cannot open {path}")
    tree = f.Get(tree_name)
    if not tree:
        f.Close()
        raise IOError(f"no {tree_name} tree in {path}")

    left_buf = _array.array("f", [0.0] * config.N_STRIPS)
    right_buf = _array.array("f", [0.0] * config.N_STRIPS)
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("Left_0_17_dE", 1)
    tree.SetBranchStatus("RightdE", 1)
    tree.SetBranchAddress("Left_0_17_dE", left_buf)
    tree.SetBranchAddress("RightdE", right_buf)

    n = tree.GetEntries()
    left = np.empty((n, config.N_STRIPS), dtype=np.float32)
    right = np.empty((n, config.N_STRIPS), dtype=np.float32)
    for i in range(n):
        tree.GetEntry(i)
        left[i] = np.frombuffer(left_buf, dtype=np.float32)
        right[i] = np.frombuffer(right_buf, dtype=np.float32)
    f.Close()

    np.savez(cache, left=left, right=right)
    return left, right


def _assemble_totals(left, right):
    """Configured totals view from (n, 18) left/right arrays.

    The first block is always the long-side trace: guards 0/17 are
    single-ended (left + right, one of which is zero) and split strips
    1-16 take their long end (L_odd / R_even). Guard columns are dropped
    when INCLUDE_GUARD_STRIPS is off. With INCLUDE_SHORT_STRIPS a second
    16-column block of the short ends (R_odd / L_even, strips 1-16)
    follows -- the two ends stay separate columns rather than being
    summed, so the per-strip L/R sharing the sum would erase remains
    visible to the model (config.block_widths() describes the layout).
    """
    long_side = left + right
    short_side = np.empty((left.shape[0], 16), dtype=left.dtype)
    for s in range(1, 17):
        if s % 2 == 1:
            long_side[:, s] = left[:, s]
            short_side[:, s - 1] = right[:, s]
        else:
            long_side[:, s] = right[:, s]
            short_side[:, s - 1] = left[:, s]
    if not config.INCLUDE_GUARD_STRIPS:
        long_side = long_side[:, 1:17]
    if config.INCLUDE_SHORT_STRIPS:
        return np.concatenate([long_side, short_side], axis=1)
    return long_side


def load_totals(path, tree_name="events_MeV"):
    """Configured per-strip sim totals (n_events, n_strip_inputs())."""
    left, right = load_left_right(path, tree_name)
    return _assemble_totals(left, right)


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


def load_experimental_totals(max_files=None, max_events_per_file=None):
    """Calibrated per-strip totals from the experimental events files.

    Per subfile: the events tree's raw ADC arrays plus the same file's
    one-row calibration tree (GainLeft/GainRight), applied per channel
    exactly like EnergyView::Decode, then the same configured totals view
    as the sim (_assemble_totals). The result is already in normalized
    a.u. (the calibration puts the beam at 1 per strip), so downstream
    code uses UNIT gains -- the sim-beam flatten does not apply here.
    """
    from analysis_utilities.io import load_leaf_array_data
    files = list_event_files()
    if max_files is not None:
        files = files[:max_files]
    cache = str(config.CACHE_DIR)
    out = []
    for path in files:
        ev = load_leaf_array_data(path,
                                  "events", ["Left_0_17_dE", "RightdE"],
                                  max_events=max_events_per_file,
                                  cache_dir=cache)
        cal = load_leaf_array_data(path,
                                   "calibration", ["GainLeft", "GainRight"],
                                   cache_dir=cache)
        gain_left = cal["GainLeft"][0].astype(np.float32)
        gain_right = cal["GainRight"][0].astype(np.float32)
        left = ev["Left_0_17_dE"].astype(np.float32) * gain_left
        right = ev["RightdE"].astype(np.float32) * gain_right
        out.append(_assemble_totals(left, right))
    return np.concatenate(out)


def find_beam_spec(specs):
    for spec in specs:
        if spec.base == "beam" and spec.strip == -1:
            return spec
    raise FileNotFoundError(
        "no beam sim file among specs (expected tag 'beam_eres')")


def beam_gains(specs=None):
    """Per-strip flatten gains from the sim beam file.

    gain[s] = FLAT_TARGET / <mean beam total[s]> over events with
    total[s] > 0; strips with no positive beam signal keep gain 0 and drop
    out, exactly like StripSumScatter::SimBeamGains.
    """
    if specs is None:
        specs = list_sim_specs(eres_only=True)
    totals = load_totals(find_beam_spec(specs).path)
    gain = np.zeros(totals.shape[1], dtype=np.float64)
    for s in range(totals.shape[1]):
        col = totals[:, s]
        pos = col[col > 0.0]
        if pos.size > 0 and pos.sum() > 0.0:
            gain[s] = 1 / pos.mean()
    return gain


def max_normalize(totals, gain):
    """Flatten per strip, then optionally max-deviation-normalize per trace.

    Returns (X, kept): X is float32 (n_kept, n_strip_inputs()); kept is the
    boolean row mask (rows whose flattened max is <= 0 are dropped either
    way). Flattened a.u. (beam == 1), preserving absolute amplitude.
    """
    flat = totals.astype(np.float64) * gain[np.newaxis, :]
    peak = flat.max(axis=1)
    kept = peak > 0.0
    X = flat[kept]
    return X.astype(np.float32), kept


def derivative(X):
    """Adjacent-strip differences of the trace, same scale as the trace.

    d[i] = X[i+1] - X[i] with no further normalization, so a flat beam keeps
    its derivative small (noise-sized) while a reaction jump stays large.
    Taken within each block of config.block_widths() (long-side trace,
    then the short-side block when present), so no difference spans the
    long/short boundary.
    """
    widths = config.block_widths()
    if len(widths) == 1:
        return np.diff(X, axis=1).astype(np.float32)
    parts = []
    lo = 0
    for w in widths:
        parts.append(np.diff(X[:, lo:lo + w], axis=1))
        lo += w
    return np.concatenate(parts, axis=1).astype(np.float32)


def ae_features(X):
    """Autoencoder input: the trace, plus the derivative appended flat.

    With SUBTRACT_BEAM_LEVEL the trace block is centered at 0 (beam noise
    around 0); the derivative is shift-invariant either way.
    """
    trace = X - np.float32(1.0) if config.SUBTRACT_BEAM_LEVEL else X
    if not config.INCLUDE_DERIVATIVE:
        return trace
    return np.concatenate([trace, derivative(X)], axis=1)


def load_population(spec, gain):
    """Max-normalized traces for one sim population."""
    X, _ = max_normalize(load_totals(spec.path), gain)
    return X
