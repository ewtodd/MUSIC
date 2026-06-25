"""Blind, transferable (a,n) clustering pipeline.

A fully data-driven look at MUSIC #DeltaE-vs-strip traces: no sim, no (a,n)
template, no hand-drawn cuts. Both stages share the continuous features
(energy above beam #Sigma(#DeltaE-1), trigger 3-sum peak3 -- 0 for no-trigger
beam, post-trigger plateau #Sigma_{trig+1..+POST}(#DeltaE-1), raw end strip #DeltaE(s17),
both-channel multiplicity):

  Step 0  every calibrated event (zero cuts, optional beam gate / pileup
          rejection in config), max-normalized (beam ~ 1/strip).
  Step *  the reaction strip is the NO-FALLBACK leading-edge onset (first strip
          crossing N-sigma of the beam noise RMS); flat beam gets no trigger
          (reac=-1, peak3=0) and drops out of the per-strip slices.
  Step 1  ONE global k=2 split -- keep the higher-plateau-excess cluster
          (real deposits), drop the flat background / beam / low-energy noise.
  Step 2  partition the kept half by leading-edge reaction strip, then
          auto-k cluster WITHIN each strip (the plateau-up / end-collapse
          (a,n) signature, plus the near-reaction both-channel multiplicity,
          resolve (a,n) vs (a,a') there).

cluster_per_reaction_strip is shared with blind_combined.py (the single-step
variant that skips step 1). Each step emits cluster scatters + one per-cluster
mean-trace overlay (StripSumScatter DrawRegionTraces style, beam reference
overlaid).

Run inside the dataset dev shell, from python/:

    python -m music_ml.blind_an
"""

import math
import random

import numpy as np
from sklearn.mixture import GaussianMixture

import config, data, mlplots

_EPS = 1.0e-12


def _fit_subsample(n):
    """Row indices to FIT on: a seeded random subsample capped at
    config.BLIND_FIT_CAP (assignment is then a single predict pass over all
    rows, so fit cost is decoupled from the reservoir size). cap None = fit on
    EVERY row (no subsample)."""
    cap = config.BLIND_FIT_CAP
    if cap is not None and n > cap:
        return np.random.default_rng(config.SEED).choice(n, cap, replace=False)
    return np.arange(n)


def _noise_fit_subsample(n):
    """Row indices to FIT on for noise-capable clustering. Uses
    config.BLIND_NOISE_FIT_CAP (separate from BLIND_FIT_CAP used by GMM).
    cap None = fit on EVERY row."""
    cap = config.BLIND_NOISE_FIT_CAP
    if cap is not None and n > cap:
        return np.random.default_rng(config.SEED).choice(n, cap, replace=False)
    return np.arange(n)


def _torch_device():
    """The torch device for the GPU DBSCAN backend, per config.BLIND_TORCH_DEVICE
    ("auto" -> cuda when available, else cpu)."""
    import torch
    dev = config.BLIND_TORCH_DEVICE
    if dev == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return dev


def _knee_index(k_sorted, min_samples, n):
    """Index of the k-distance elbow by the parameter-free max-distance-to-chord
    (Kneedle/triangle) method, over the WHOLE descending curve: draw the chord
    from the first point (rank 0, max k-distance) to the last (rank n-1, min),
    and take the rank whose perpendicular distance from that chord is largest --
    the point where the steep outlier drop meets the dense bulk plateau. No
    window and no floor (the old two-segment OLS searched only the top ~1000
    ranks, cutting off the real elbow and settling near min_samples+1). To keep a
    handful of extreme-isolated points from stretching the chord, the top end is
    anchored at the 99.5th percentile rather than the literal max. `k_sorted` is
    the per-point distance to the min_samples-th nearest neighbour, descending."""
    y = k_sorted.astype(np.float64)
    m = y.shape[0]
    if m < 3:
        return m - 1
    # Anchor the high end below the most extreme outliers so they don't dominate
    # the chord; clamp the curve to that ceiling.
    hi = float(np.percentile(y, 99.5))
    y = np.minimum(y, hi)
    x = np.arange(m, dtype=np.float64)
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    denom = math.hypot(dx, dy)
    if denom < 1.0e-15:  # flat curve -> no elbow
        return m - 1
    # Perpendicular distance of every point to the chord; the elbow is the max.
    dist = np.abs(dy * (x - x[0]) - dx * (y - y[0])) / denom
    return int(np.argmax(dist))


def _plot_kdist(k_sorted, knee_idx, eps, min_samples, tag, subdir):
    """Plot the descending k-distance graph with the chosen knee (eps) marked,
    mirroring the old auto-eps diagnostic."""
    n = k_sorted.shape[0]
    R = mlplots._root()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    g = R.TGraph(n)
    for i in range(n):
        g.SetPoint(i, float(i), float(k_sorted[i]))
    g.SetTitle(f"k-distance graph (torch_dbscan); k={min_samples}; "
               f"auto-eps={eps:.4f} (knee idx {knee_idx});"
               f" min={k_sorted[-1]:.4f} max={k_sorted[0]:.4f}")
    g.GetXaxis().SetTitle("rank (descending)")
    g.GetYaxis().SetTitle(f"k-distance (k={min_samples}) [z-scored]")
    g.SetLineColor(R.kBlue)
    g.SetLineWidth(2)
    g.Draw("ALP")
    knee = R.TMarker(float(knee_idx), float(k_sorted[knee_idx]), 4)
    knee.SetMarkerColor(R.kRed)
    knee.SetMarkerSize(2)
    knee.Draw()
    R.PlottingUtils.SaveFigure(c, f"blind_autoeps_{tag}", subdir,
                               R.PlotSaveOptions.kLINEAR)


def _torch_auto_eps(Z, min_samples, tag, subdir, device):
    """Auto-detect eps for torch_dbscan via the k-distance elbow, the whole
    k-distance computed on the GPU in tiles. Z is the z-scored fit subsample.
    Returns the chosen eps (falls back to 1.0 when too few points)."""
    import torch
    n = Z.shape[0]
    if n < 2 * min_samples:
        print(f"  auto-eps {tag}: n={n} < 2*min_samples={2*min_samples}, "
              f"using default eps=1.0")
        return 1.0
    t = torch.as_tensor(Z, device=device)
    tile = config.BLIND_TORCH_TILE
    kd = torch.empty(n, device=device)
    for i in range(0, n, tile):
        d = torch.cdist(t[i:i + tile], t)
        # k-th smallest distance INCLUDING self (the self 0 is the smallest), so
        # k=min_samples gives the (min_samples-1)-th true neighbour -- matching
        # sklearn's kneighbors(n_neighbors=min_samples)[:, -1].
        kd[i:i + tile] = torch.kthvalue(d, k=min_samples, dim=1).values
    k_sorted = torch.sort(kd, descending=True).values.cpu().numpy()
    knee_idx = _knee_index(k_sorted, min_samples, n)
    eps = float(k_sorted[knee_idx])
    print(f"  auto-eps {tag}: min_samples={min_samples}, eps={eps:.4f} "
          f"(knee at index {knee_idx}/{n}, "
          f"k-dist={k_sorted[knee_idx]:.4f})")
    _plot_kdist(k_sorted, knee_idx, eps, min_samples, tag, subdir)
    return eps


def _torch_dbscan_cores(Zfit, eps, min_samples, device):
    """Fit DBSCAN cores on the z-scored subsample Zfit: a point is a core if it
    has >= min_samples points (including itself) within eps. Cores within eps of
    one another are linked, and connected components over that core graph are the
    clusters. Returns (core_pts (nc, d) float32, core_labels (nc,) in 0..k-1,
    k). All pairwise work is tiled on the GPU; the components run on the CPU
    (nc <= the fit cap, so the core graph is small)."""
    import torch
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    t = torch.as_tensor(Zfit, device=device)
    n = t.shape[0]
    tile = config.BLIND_TORCH_TILE
    counts = torch.zeros(n, dtype=torch.int64, device=device)
    for i in range(0, n, tile):
        d = torch.cdist(t[i:i + tile], t)
        counts[i:i + tile] = (d <= eps).sum(dim=1)
    core = counts >= min_samples
    core_idx = torch.nonzero(core, as_tuple=False).squeeze(1)
    nc = int(core_idx.shape[0])
    if nc == 0:
        return (np.zeros((0, Zfit.shape[1]),
                         dtype=np.float32), np.zeros(0, dtype=np.int64), 0)
    tc = t[core_idx]
    rows, cols = [], []
    for i in range(0, nc, tile):
        d = torch.cdist(tc[i:i + tile], tc)
        r, c = torch.nonzero(d <= eps, as_tuple=True)
        rows.append((r + i).cpu().numpy())
        cols.append(c.cpu().numpy())
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    graph = coo_matrix((np.ones(rows.size, dtype=np.int8), (rows, cols)),
                       shape=(nc, nc))
    k, comp = connected_components(graph, directed=False)
    return tc.cpu().numpy(), comp.astype(np.int64), int(k)


def _torch_assign(Zall, core_pts, core_labels, eps, device):
    """Assign every z-scored row in Zall to its nearest core sample on the GPU:
    label = that core's cluster when the distance is <= eps, else -1 (noise).
    Tiled to cap peak memory; this single pass replaces any out-of-sample
    predict (DBSCAN has none)."""
    import torch
    tcore = torch.as_tensor(core_pts, device=device)
    clab = torch.as_tensor(core_labels, device=device)
    n = Zall.shape[0]
    tile = config.BLIND_TORCH_TILE
    out = np.empty(n, dtype=np.int64)
    for i in range(0, n, tile):
        q = torch.as_tensor(Zall[i:i + tile], device=device)
        d = torch.cdist(q, tcore)
        md, mi = d.min(dim=1)
        lab = torch.where(md <= eps, clab[mi], torch.full_like(clab[mi], -1))
        out[i:i + tile] = lab.cpu().numpy()
    return out


_GPU_FREED = False  # llama-swap eject runs at most once per process


def _free_gpu_via_llama_swap():
    """Best-effort: if a local llama-swap instance has a model resident on the
    GPU, ask it to unload so the VRAM is free for the torch_dbscan assign.
    llama-swap respawns the model on its next request, so FIM coding keeps
    working -- it just reloads. Silent no-op when llama-swap is unreachable or
    nothing is loaded. Uses only the stdlib (no requests dependency)."""
    import json
    import urllib.request
    url = config.BLIND_LLAMA_SWAP_URL
    if not url:
        return
    try:
        with urllib.request.urlopen(f"{url}/running", timeout=3) as resp:
            running = json.loads(resp.read()).get("running", [])
    except Exception:
        return  # not running / unreachable -> nothing to free
    active = [m for m in running if m.get("state") != "stopped"]
    if not active:
        return
    names = ", ".join(m.get("model", "?") for m in active)
    try:
        urllib.request.urlopen(f"{url}/unload", timeout=15).read()
    except Exception as exc:
        print(f"  WARNING: llama-swap unload failed ({exc}); the GPU may be "
              f"short on memory")
        return
    # /unload only INITIATES the stop; poll /running until the model is gone so
    # the VRAM is actually released before torch allocates (up to ~10 s).
    import time
    for _ in range(20):
        try:
            with urllib.request.urlopen(f"{url}/running", timeout=3) as resp:
                still = [
                    m for m in json.loads(resp.read()).get("running", [])
                    if m.get("state") != "stopped"
                ]
        except Exception:
            still = []
        if not still:
            break
        time.sleep(0.5)
    print(f"  free GPU: unloaded llama-swap model(s) [{names}] via "
          f"{url}/unload (reloads on its next request)")


def _ensure_gpu_free():
    """Run the llama-swap eject once per process, before the first GPU work,
    when config.BLIND_FREE_GPU is on."""
    global _GPU_FREED
    if _GPU_FREED:
        return
    _GPU_FREED = True
    if config.BLIND_FREE_GPU:
        _free_gpu_via_llama_swap()


def cluster_torch_dbscan(feats, tag="", subdir=config.PLOT_SUBDIR):
    """GPU DBSCAN backend for cluster_auto: z-score the features, fit cores on a
    capped subsample (config.BLIND_NOISE_FIT_CAP), then assign EVERY row to its
    nearest core within eps on the GPU (label -1 = noise). Returns (labels (n,),
    means in ORIGINAL feature units (k, d), k, []) -- the same shape as
    cluster_auto's GMM path; the BIC list is empty (not applicable)."""
    _ensure_gpu_free()  # free any resident llama-swap model before GPU work
    F = feats.reshape(feats.shape[0], -1).astype(np.float64)
    fit_idx = _noise_fit_subsample(F.shape[0])
    mu = F[fit_idx].mean(axis=0)
    sd = F[fit_idx].std(axis=0) + 1.0e-9
    Z = ((F - mu) / sd).astype(np.float32)
    Zfit = Z[fit_idx]
    device = _torch_device()
    min_samples = config.BLIND_DBSCAN_MIN_SAMPLES
    eps = config.BLIND_DBSCAN_EPS
    if eps is None:
        eps = _torch_auto_eps(Zfit, min_samples, tag, subdir, device)
    core_pts, core_labels, k = _torch_dbscan_cores(Zfit, eps, min_samples,
                                                   device)
    print(f"  torch_dbscan {tag}: device={device}, eps={eps:.4f}, "
          f"min_samples={min_samples}, {core_pts.shape[0]} core pts, "
          f"fit on {Zfit.shape[0]} of {F.shape[0]} rows -> k={k}")
    if k == 0:
        print("  torch_dbscan: no cores found; all rows -> noise (-1)")
        return (np.full(F.shape[0], -1, dtype=np.int64),
                np.zeros((0, F.shape[1]), dtype=np.float64), 0, [])
    labels = _torch_assign(Z, core_pts, core_labels, eps, device)
    # Cluster means in ORIGINAL feature units, from the fit subsample's own
    # assignment (cheap; avoids a second pass over the full reservoir).
    fit_lab = labels[fit_idx]
    means = np.zeros((k, F.shape[1]), dtype=np.float64)
    for j in range(k):
        m = fit_lab == j
        if m.any():
            means[j] = F[fit_idx][m].mean(axis=0)
    print(f"  torch_dbscan {tag}: assigned {int((labels >= 0).sum())} of "
          f"{labels.size} rows, {int((labels == -1).sum())} noise")
    return labels, means, k, []


def _trace(view):
    """Max-normalized trace view (unit gains: the data calibration already
    puts the beam at 1 per strip), plus the kept-row mask."""
    unit = np.ones(view.shape[1], dtype=np.float64)
    return data.max_normalize(view, unit)


def cluster_auto(feats,
                 k=None,
                 tag="",
                 subdir=config.PLOT_SUBDIR,
                 backend=None,
                 noise_pctl="config"):
    """Cluster a (n, d) numpy feature array. With k=None the cluster count is
    GUESSED (the k in 1..config.BLIND_MAX_K minimizing GMM BIC); pass an
    integer k to force that count.

    `backend` overrides config.BLIND_NOISE_CLUSTERING for this call (None =
    follow config) so a caller can pin its method -- e.g. the (a,n) sub-split
    pins 'gmm' regardless of the global switch. When the resolved backend is
    'torch_dbscan', delegates to cluster_torch_dbscan() (GPU DBSCAN that tags
    noise as label -1). Otherwise ('gmm'/'none') fits a full-covariance GMM on a
    capped subsample via sklearn, then assigns EVERY row in one predict pass --
    optionally tagging the lowest-density tail as noise. `noise_pctl` controls
    that: "config" = use config.BLIND_GMM_NOISE_PCTL, None = no tagging, a float
    = that bottom percentile. Returns (labels (n,), means in ORIGINAL feature
    units (k,d), k_chosen, [bic per k])."""
    method = backend if backend is not None else config.BLIND_NOISE_CLUSTERING
    if method == "torch_dbscan":
        return cluster_torch_dbscan(feats, tag=tag, subdir=subdir)
    pctl = config.BLIND_GMM_NOISE_PCTL if noise_pctl == "config" else noise_pctl
    F = feats.reshape(feats.shape[0], -1).astype(np.float64)
    fit_idx = _fit_subsample(F.shape[0])
    mu = F[fit_idx].mean(axis=0)
    sd = F[fit_idx].std(axis=0) + 1.0e-9
    Z = (F - mu) / sd
    Zfit = Z[fit_idx]
    Zall = Z
    n_fit, d = Zfit.shape
    if k is not None:
        candidates = [k]
    elif config.BLIND_MAX_K < 0:
        candidates = range(1, 20)  # unbounded cap — BIC penalises overfit
    else:
        candidates = range(1, config.BLIND_MAX_K + 1)
    bics = []
    best = None  # (bic, k, gmm)
    for k in candidates:
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                max_iter=200,
                n_init=1,  # deterministic init; no extra starts needed
                reg_covar=1.0e-5,
                random_state=config.SEED,
            )
            gmm.fit(Zfit)
            loglik = float(gmm.score(Zfit) * n_fit)  # total log-likelihood
        except Exception:
            bics.append(float("nan"))
            continue
        # BIC = p * log(n) - 2 * loglik
        # p = k*d (means) + k*d*(d+1)/2 (covariances) + (k-1) (weights)
        p = k * d + k * d * (d + 1) // 2 + (k - 1)
        bic = p * math.log(max(n_fit, 1)) - 2.0 * loglik
        bics.append(bic)
        if best is None or bic < best[0]:
            best = (bic, k, gmm)
    if best is None:  # every fit failed -> one cluster
        return np.zeros(F.shape[0], dtype=np.int64), mu[None, :], 1, bics
    _, k, best_gmm = best
    best_gmm = best_gmm.fit(Zfit)
    labels = best_gmm.predict(Zall).astype(np.int64)
    if pctl is not None:
        # Tag the lowest-density tail as noise: per-point log-likelihood below
        # the chosen percentile of the fit subsample's own logL distribution.
        thr = float(np.percentile(best_gmm.score_samples(Zfit), pctl))
        noise = best_gmm.score_samples(Zall) < thr
        labels[noise] = -1
        print(f"  GMM noise: tagged {int(noise.sum())} of {labels.size} rows "
              f"(logL < {pctl:g}th pctl = {thr:.2f}) as -1")
    return labels, best_gmm.means_ * sd + mu, k, bics


# Trace features.
def reaction_strip(X, baseline, margin):
    """Reaction strip per event by a cheap CONSTANT-FRACTION onset: the first
    strip (from strip 1 on) whose beam-baseline-subtracted excess reaches
    config.BLIND_REAC_ONSET_FRAC of the event's OWN peak excess. Unlike a
    leading-edge crossing of a fixed level this does NOT walk with amplitude (a
    bigger deposit no longer triggers earlier). The absolute `margin` (N-sigma
    beam noise) is only the trigger-EXISTS gate: an event whose peak excess
    does not clear it has NO trigger and gets -1 (flat beam), so it drops out
    of the per-strip slices. baseline is the per-strip beam level (mean
    pure-beam trace, _beam_reference) -- subtracting it removes the L/R
    sawtooth. Triggered strips are absolute detector indices via the guard
    offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    ex = X.astype(np.float64) - np.asarray(baseline, dtype=np.float64)
    peak = ex[:, 1:].max(axis=1)  # peak excess over strips 1+ (skip guard 0)
    has = peak > margin  # trigger exists only if the peak clears N-sigma
    thr = config.BLIND_REAC_ONSET_FRAC * np.maximum(peak, 0.0)
    above = ex >= thr[:, np.newaxis]
    above[:, 0] = False  # scan starts at strip 1, like the C++
    onset_col = np.argmax(above, axis=1)  # first strip reaching frac * peak
    return np.where(has, onset_col + first, -1)


def _reaction_strip_all(X, beam_ref, beam_sigma):
    """The constant-fraction reaction strip for every event: trigger-exists
    gate at config.BLIND_REAC_ONSET_NSIGMA * beam noise RMS, position at
    config.BLIND_REAC_ONSET_FRAC of each event's peak excess. Reported once;
    -1 marks no trigger (flat beam) -- those events drop out of the slices."""
    margin = config.BLIND_REAC_ONSET_NSIGMA * beam_sigma
    reac = reaction_strip(X, beam_ref, margin)
    print(f"  reaction onset: gate {config.BLIND_REAC_ONSET_NSIGMA:g}-sigma = "
          f"{margin:.4f}, CF fraction {config.BLIND_REAC_ONSET_FRAC:g}; "
          f"triggered {int((reac >= 0).sum())} of {X.shape[0]} "
          f"(no trigger: {int((reac < 0).sum())})")
    return reac


def _tail_value(X):
    """RAW #DeltaE at the last strip (strip 17), beam ~ 1. (a,n) collapses
    below beam at end of range so it drops clearly under 1. Strip maps to a
    column via the guard offset (falls back to the last strip when 17 is
    outside the view)."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    col = 17 - first
    if not 0 <= col < X.shape[1]:
        col = X.shape[1] - 1
    return X.astype(np.float64)[:, col]


def _plateau_excess(X, reac):
    """Beam-subtracted plateau over a SLIDING window that tracks the trigger:
    #Sigma(#DeltaE - 1) over the config.BLIND_PLATEAU_POST strips immediately
    after the reaction strip (reac+1 .. reac+POST). Beam ~ 0; an (a,n) rides
    above beam right after the reaction so it is positive there regardless of
    WHERE the reaction sat -- a fixed 2-10 band only caught early triggers. 0
    where there is no trigger (reac < 0). Window strips past the end of the view
    are dropped from the sum (slices where the window never fully fits are
    skipped upstream in cluster_per_reaction_strip). Strips map to columns via
    the guard offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    Xd = X.astype(np.float64)
    reac = np.asarray(reac)
    has = reac >= 0
    rows = np.arange(Xd.shape[0])
    out = np.zeros(Xd.shape[0], dtype=np.float64)
    for d in range(1, config.BLIND_PLATEAU_POST + 1):
        col = reac + d - first
        in_range = has & (col >= 0) & (col < Xd.shape[1])
        out += np.where(in_range,
                        Xd[rows, np.clip(col, 0, Xd.shape[1] - 1)] - 1.0, 0.0)
    return out


# Pipeline steps.
# Beam-gate fit (StripSumScatter FindBeamGate): histogram a strip pair,
# locate the beam peak, take threshold-weighted moments over the bins around
# it, accept within N-sigma. Beam is fit, not assumed at 1 (so the L_odd /
# R_even per-strip offset is handled).
_BEAM_BINS = 240
_BEAM_LO, _BEAM_HI = 0.0, 3.0
_BEAM_SEED_HALF = 40  # half-window (bins) for the moment sum around the peak
_BEAM_SEED_FRAC = 0.30  # ignore bins below this fraction of the peak height


def _fit_beam_ellipse(a, b):
    """(mu_a, mu_b, sigma_a, sigma_b) of the beam blob in the (a, b) plane
    via threshold-weighted moments around the 2-D histogram peak, or None if
    too sparse. Sigma is floored at two bin widths (like ComputeMoments)."""
    ok = (a > _BEAM_LO) & (b > _BEAM_LO)
    h, _, _ = np.histogram2d(a[ok],
                             b[ok],
                             bins=_BEAM_BINS,
                             range=((_BEAM_LO, _BEAM_HI), (_BEAM_LO,
                                                           _BEAM_HI)))
    if h.sum() < 100:
        return None
    bw = (_BEAM_HI - _BEAM_LO) / _BEAM_BINS
    centers = _BEAM_LO + bw * (np.arange(_BEAM_BINS) + 0.5)
    ba, bb = np.unravel_index(np.argmax(h), h.shape)
    lo_a, hi_a = max(0, ba - _BEAM_SEED_HALF), min(_BEAM_BINS - 1,
                                                   ba + _BEAM_SEED_HALF)
    lo_b, hi_b = max(0, bb - _BEAM_SEED_HALF), min(_BEAM_BINS - 1,
                                                   bb + _BEAM_SEED_HALF)
    w = h[lo_a:hi_a + 1, lo_b:hi_b + 1].copy()
    w[w < _BEAM_SEED_FRAC * h[ba, bb]] = 0.0
    weight = w.sum()
    if weight <= 0.0:
        return None
    ca = centers[lo_a:hi_a + 1][:, np.newaxis]
    cb = centers[lo_b:hi_b + 1][np.newaxis, :]
    mu_a = (w * ca).sum() / weight
    mu_b = (w * cb).sum() / weight
    sig_a = np.sqrt((w * (ca - mu_a)**2).sum() / weight)
    sig_b = np.sqrt((w * (cb - mu_b)**2).sum() / weight)
    return (float(mu_a), float(mu_b), max(float(sig_a), 2.0 * bw),
            max(float(sig_b), 2.0 * bw))


def _beam_gate(X, pairs):
    """Keep events inside the fitted beam ellipse of EVERY strip pair in
    `pairs` (e.g. ((0, 1),) or ((0, 1), (16, 17))). Strips map to columns via
    the guard offset; a pair with a missing column or failed fit is skipped."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    Xd = X.astype(np.float64)
    keep = np.ones(X.shape[0], dtype=bool)
    ns = config.BLIND_BEAM_NSIGMA
    for sa, sb in pairs:
        ca, cb = sa - first, sb - first
        if not (0 <= ca < X.shape[1] and 0 <= cb < X.shape[1]):
            continue
        fit = _fit_beam_ellipse(Xd[:, ca], Xd[:, cb])
        if fit is None:
            continue
        mu_a, mu_b, sig_a, sig_b = fit
        da = (Xd[:, ca] - mu_a) / (ns * sig_a)
        db = (Xd[:, cb] - mu_b) / (ns * sig_b)
        keep &= (Xd[:, ca] > _BEAM_LO) & (Xd[:, cb] > _BEAM_LO) \
            & (da * da + db * db < 1.0)
    return keep


_GATE_FIT = None  # (s0,s1) beam-ellipse fit, set once on the full reservoir


def _set_gate_fit(X):
    """Fit the (s0,s1) beam ellipse ONCE on the full reservoir and cache it, so
    the `beamgate` feature scores membership relative to the whole-reservoir
    beam (not the per-slice data)."""
    global _GATE_FIT
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    ca, cb = 0 - first, 1 - first
    if 0 <= ca < X.shape[1] and 0 <= cb < X.shape[1]:
        Xd = X.astype(np.float64)
        _GATE_FIT = _fit_beam_ellipse(Xd[:, ca], Xd[:, cb])


def _beam_gate_value(X):
    """Binary `beamgate` feature: 1 if the event is INSIDE the cached (s0,s1)
    beam ellipse (entered as beam), 0 otherwise -- the beam-entrance gate used
    as a clustering feature instead of a hard cut. 0 everywhere if the ellipse
    was not fit. Strips map to columns via the guard offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    ca, cb = 0 - first, 1 - first
    Xd = X.astype(np.float64)
    if (_GATE_FIT is None
            or not (0 <= ca < Xd.shape[1] and 0 <= cb < Xd.shape[1])):
        return np.zeros(Xd.shape[0], dtype=np.float64)
    mu_a, mu_b, sig_a, sig_b = _GATE_FIT
    ns = config.BLIND_BEAM_NSIGMA
    da = (Xd[:, ca] - mu_a) / (ns * sig_a)
    db = (Xd[:, cb] - mu_b) / (ns * sig_b)
    inside = ((Xd[:, ca] > _BEAM_LO) & (Xd[:, cb] > _BEAM_LO)
              & (da * da + db * db < 1.0))
    return inside.astype(np.float64)


_PREBEAM_FITS = None  # {strip: (mu_a,mu_b,sig_a,sig_b)} beam ellipse per (s-1,s-2)


def _set_prebeam_fits(X):
    """Fit the beam ellipse in the (s-1, s-2) plane for every candidate trigger
    strip s (from config.BLIND_PREBEAM_MIN_STRIP up) ONCE on the full reservoir,
    cached for the `prebeam` feature/cut -- so membership scores against the
    whole-reservoir beam blob at those pre-trigger strips, not the per-slice
    data. Strips map to columns via the guard offset."""
    global _PREBEAM_FITS
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    Xd = X.astype(np.float64)
    fits = {}
    for s in range(config.BLIND_PREBEAM_MIN_STRIP, 18):
        ca, cb = (s - 1) - first, (s - 2) - first
        if 0 <= cb and ca < Xd.shape[1]:
            fits[s] = _fit_beam_ellipse(Xd[:, ca], Xd[:, cb])
    _PREBEAM_FITS = fits


def _pre_trigger_beam(X, reac):
    """Binary `prebeam` feature: 1 if the event is STILL inside the beam blob at
    the two strips just before its trigger (#DeltaE(reac-1), #DeltaE(reac-2)), 0
    otherwise. A clean reaction is beam-like right up to the trigger; pileup/junk
    has already departed the beam by then. Applies to every trigger at or after
    config.BLIND_PREBEAM_MIN_STRIP (floors at 2, where reac-1, reac-2 = strips 1
    and 0, the guard -- both valid columns); triggers below it or with no valid
    pre-trigger pair return 1 (pass, "still beam"). The per-pair beam ellipse is
    the cached _PREBEAM_FITS; strips map to columns via the guard offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    Xd = X.astype(np.float64)
    reac = np.asarray(reac)
    ns = config.BLIND_BEAM_NSIGMA
    out = np.ones(Xd.shape[0], dtype=np.float64)  # default: pass (still beam)
    if _PREBEAM_FITS is None:
        return out
    for s in _PREBEAM_FITS:
        fit = _PREBEAM_FITS[s]
        sel = reac == s
        if fit is None or not bool(sel.any()):
            continue
        ca, cb = (s - 1) - first, (s - 2) - first
        if not (0 <= cb and ca < Xd.shape[1]):
            continue
        mu_a, mu_b, sig_a, sig_b = fit
        da = (Xd[sel, ca] - mu_a) / (ns * sig_a)
        db = (Xd[sel, cb] - mu_b) / (ns * sig_b)
        inside = ((Xd[sel, ca] > _BEAM_LO) & (Xd[sel, cb] > _BEAM_LO)
                  & (da * da + db * db < 1.0))
        out[sel] = inside.astype(np.float64)
    return out


def _apply_prebeam_cut(X, both, reac):
    """Optional hard cut (config.BLIND_REJECT_PREBEAM): drop events whose
    pre-trigger pair is NOT still beam (prebeam == 0). Needs `reac`, so it runs
    after the reaction strip rather than in step 0. Returns (X, both, reac),
    filtered or unchanged."""
    if not config.BLIND_REJECT_PREBEAM:
        return X, both, reac
    keep = _pre_trigger_beam(X, reac) > 0.5
    print(
        f"  pre-trigger beam cut (triggers >= s{config.BLIND_PREBEAM_MIN_STRIP}"
        f" still beam at reac-1,reac-2): kept {int(keep.sum())} of "
        f"{X.shape[0]}")
    return X[keep], both[keep], reac[keep]


def _is_pileup(X):
    """StripSumScatter IsPileup: True for events with at least
    config.BLIND_PILEUP_MIN_STRIPS long strips (1-16) at or above
    config.BLIND_PILEUP_THRESH (overlapping-beam pileup). Strips map to
    columns via the guard offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    cols = [s - first for s in range(1, 17) if 0 <= s - first < X.shape[1]]
    over = (X.astype(np.float64)[:, cols]
            >= config.BLIND_PILEUP_THRESH).sum(axis=1)
    return over >= config.BLIND_PILEUP_MIN_STRIPS


def _is_noise(X):
    """Hard noise flag: True for events with at least
    config.BLIND_NOISE_MIN_STRIPS strips (all columns 0..n-1, matching the
    long-side trace view) below config.BLIND_NOISE_THRESH. A real deposit is
    never uniformly sub-threshold; flat-noise drops are. Strips map to
    columns via the guard offset."""
    Xd = X.astype(np.float64)
    below = (Xd < config.BLIND_NOISE_THRESH).sum(axis=1)
    return below >= config.BLIND_NOISE_MIN_STRIPS


def _is_offbeam(X, beam_ref=None):
    """Off-beam junk flag, consolidating pileup + noise: True for events with
    at least config.BLIND_OFFBEAM_MIN_STRIPS long strips (1-16) at least
    config.BLIND_OFFBEAM_DIST (absolute) from beam -- elevated (pileup) OR
    below (dropout). A real (a,n) departs over only a few strips (reaction +
    collapse), so MIN_STRIPS set high leaves it at 0. beam_ref per strip (or
    1.0). Strips map to columns via the guard offset."""
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    cols = [s - first for s in range(1, 17) if 0 <= s - first < X.shape[1]]
    Xd = X.astype(np.float64)[:, cols]
    thr = 1.0 if beam_ref is None else np.asarray(beam_ref,
                                                  dtype=np.float64)[cols]
    off = (np.abs(Xd - thr) >= config.BLIND_OFFBEAM_DIST).sum(axis=1)
    return off >= config.BLIND_OFFBEAM_MIN_STRIPS


def _beam_reference(X):
    """PURE-beam events -- entered AND left as beam (inside the fitted
    (s0,s1) AND (s16,s17) beam ellipses). Returns (baseline, sigma):
    baseline is the per-strip mean trace (the beam reference overlaid on
    every gallery), sigma is the pooled RMS of the baseline-subtracted
    excess -- the noise level the leading-edge threshold is set from.
    (None, None) if no pure-beam events are found."""
    keep = _beam_gate(X, ((0, 1), (16, 17)))
    if not keep.any():
        print("  beam reference: no pure-beam events found")
        return None, None
    Xb = X[keep].astype(np.float64)
    baseline = Xb.mean(axis=0)
    sigma = float((Xb - baseline).std())
    print(f"  beam reference: mean+RMS of {int(keep.sum())} pure-beam events "
          f"(fitted s0,s1 & s16,s17 ellipses); noise sigma={sigma:.4f}")
    return baseline, sigma


def step0_reservoir():
    """The input reservoir. By default EVERY calibrated event with ZERO cuts
    (no pileup/noise rejection, no reaction-shape selection) -- the per-channel
    calibration puts the data beam at ~1 per strip and the long-side trace is
    taken with unit gains (rows with non-positive peak dropped). When
    config.BLIND_GATE_BEAM_01 is on, a minimal strip-0/1 beam-entrance gate is
    applied here. Returns (X (n, long_w), both (n, 18) -- the per-strip
    both-channel fired mask, carried through every cut alongside X)."""
    print("step 0: calibrated events" +
          (" (strip-0/1 beam gate)" if config.
           BLIND_GATE_BEAM_01 else " (ZERO cuts)"))
    totals, both = data.load_experimental_totals(
        max_files=config.BLIND_MAX_FILES,
        max_events_per_file=config.BLIND_MAX_EVENTS_PER_FILE)
    long_w = config.block_widths()[0]
    X, kept = _trace(totals)
    X = X[:, :long_w]
    both = both[kept]
    _set_gate_fit(X)  # cache the (s0,s1) ellipse for the `beamgate` feature
    _set_prebeam_fits(X)  # cache per-(s-1,s-2) ellipses for the `prebeam` test
    if config.BLIND_GATE_BEAM_01:
        keep = _beam_gate(X, ((0, 1), ))
        print(f"  beam-0/1 gate (fitted s0,s1 ellipse, "
              f"{config.BLIND_BEAM_NSIGMA}-sigma): kept {int(keep.sum())} "
              f"of {X.shape[0]}")
        X, both = X[keep], both[keep]
    # Per-cut hard cuts, each honored independently of whether the same
    # condition is also used as a clustering feature (a cut set to both just
    # leaves survivors reading the flag as 0). beam_ref is not built yet here,
    # so offbeam compares against the flat normed beam level (1.0).
    if config.BLIND_REJECT_NOISE:
        keep = ~_is_noise(X)
        print(f"  noise cut (>= {config.BLIND_NOISE_MIN_STRIPS} strips "
              f"< {config.BLIND_NOISE_THRESH}): kept {int(keep.sum())} "
              f"of {X.shape[0]}")
        X, both = X[keep], both[keep]
    if config.BLIND_REJECT_PILEUP:
        keep = ~_is_pileup(X)
        print(f"  pileup cut (>= {config.BLIND_PILEUP_MIN_STRIPS} strips "
              f">= {config.BLIND_PILEUP_THRESH}): kept {int(keep.sum())} "
              f"of {X.shape[0]}")
        X, both = X[keep], both[keep]
    if config.BLIND_REJECT_OFFBEAM:
        keep = ~_is_offbeam(X)
        print(f"  offbeam cut (>= {config.BLIND_OFFBEAM_MIN_STRIPS} strips "
              f">= {config.BLIND_OFFBEAM_DIST} from beam): kept "
              f"{int(keep.sum())} of {X.shape[0]}")
        X, both = X[keep], both[keep]
    mult = both[:, 1:17].sum(axis=1)
    print(f"  reservoir: {X.shape[0]} events ({long_w} long-side strips)")
    print(f"  both-channel mult: {int((mult > 0).sum())} of {X.shape[0]} "
          f"events nonzero (max {int(mult.max()) if mult.size else 0})")
    return X, both


# Column indices into the FULL (step-2) feature matrix, used for the 2-D
# scatters: energy, peak3, plateau, tail, mult, ... in that build order.
_COL_ENERGY, _COL_PLATEAU, _COL_TAIL = 0, 2, 3


def _near_reaction_both(both, reac):
    """Binary per event: 1 if a both-channel firing occurred within +-1 strip
    of the reaction strip (reac-1 .. reac+1, clamped to the split strips
    1-16). On-axis deposits light only the long anode (no both-firing); a
    both-channel hit AT the reaction is an off-axis / topology signature."""
    n = both.shape[0]
    rows = np.arange(n)
    out = np.zeros(n, dtype=bool)
    for d in (-1, 0, 1):
        s = reac + d
        valid = (s >= 1) & (s <= 16)
        out |= both[rows, np.clip(s, 1, 16)] & valid
    return out


def _peak3(X, reac):
    """3-strip #DeltaE sum centered on the TRIGGER (reaction) strip (matches
    StripSumScatter peak3): #DeltaE(reac-1)+#DeltaE(reac)+#DeltaE(reac+1). 0
    where there is no trigger (reac < 0) -- flat beam scores 0, the
    discriminator. Strips map to columns via the guard offset."""
    Xd = np.pad(X.astype(np.float64), ((0, 0), (1, 1)))  # zero-pad both edges
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    ncol = X.shape[1]
    rows = np.arange(X.shape[0])
    has = reac >= 0
    base_col = reac - first  # column of the trigger strip in the unpadded view
    out = np.zeros(X.shape[0], dtype=np.float64)
    for d in (-1, 0, 1):
        c = base_col + d
        in_range = has & (c >= 0) & (c < ncol)
        out += np.where(in_range, Xd[rows, np.clip(c, 0, ncol - 1) + 1], 0.0)
    return out


def _trig_tail_dev(X, reac, beam_ref):
    """|#DeltaE - beam| at the TRIGGER (reaction) strip plus |#DeltaE - beam|
    at the end strip s17 (StripSumScatter trigtaildev): an (a,n) rises at the
    trigger AND collapses at s17, so both magnitudes add; flat beam scores ~0.
    The trigger term is 0 where there is no trigger (reac < 0). Strips map to
    columns via the guard offset."""
    Xd = X.astype(np.float64)
    br = np.asarray(beam_ref, dtype=np.float64)
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    rows = np.arange(Xd.shape[0])
    has = reac >= 0
    trig_col = np.clip(reac - first, 0, Xd.shape[1] - 1)
    tail_col = min(17 - first, Xd.shape[1] - 1)
    trig_dev = np.where(has, np.abs(Xd[rows, trig_col] - br[trig_col]), 0.0)
    tail_dev = np.abs(Xd[:, tail_col] - br[tail_col])
    return trig_dev + tail_dev


def _beam_deviation(X, beam_ref):
    """Back-half beam deviation (StripSumScatter beamdev): RMS of
    (#DeltaE - beam) over strips 8-17. LOW = beam-like (flat at beam level),
    high for a reaction plateau/collapse or pileup elevation. Amplitude-aware
    (NOT max-normalized) and sawtooth-free (beam subtraction). Strips map to
    columns via the guard offset."""
    Xd = X.astype(np.float64)
    br = np.asarray(beam_ref, dtype=np.float64)
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    cols = [s - first for s in range(8, 18) if 0 <= s - first < Xd.shape[1]]
    d = Xd[:, cols] - br[cols]
    return np.sqrt((d * d).mean(axis=1))


def _cluster_matrix(X,
                    both,
                    reac,
                    drop_flags=(),
                    beam_ref=None,
                    include=None,
                    X_sg=None,
                    reac_sg=None):
    """The clustering feature matrix (n, d) and column names. The continuous
    features (energy above beam #Sigma(#DeltaE-1), trigger 3-sum `peak3`
    centered on the reaction strip -- 0 for no-trigger events, plateau excess
    #Sigma_{2-10}(#DeltaE-1), raw end strip #DeltaE(s17), both-channel
    multiplicity -- the row sum of the per-strip mask `both`); plus a BINARY
    flag per junk cut whose own config.BLIND_*_FEATURE is on (noise / pileup /
    offbeam) so the clustering can group the garbage rather than (or as well as)
    it being hard-cut upstream; plus, when config.BLIND_NEAR_MULT_FEATURE, the
    binary near-reaction both-channel multiplicity (a firing within +-1 of the
    reaction strip); plus, when `beam_ref` is given and
    config.BLIND_TRIGTAIL_FEATURE, the trigger+tail deviation |#DeltaE-beam| at
    the reaction strip plus at s17.

    `reac` (the no-fallback trigger strip, -1 when none) is REQUIRED -- peak3
    depends on it. `include`, when given, is an AUTHORITATIVE whitelist of
    feature names (step 1 uses config.BLIND_STEP1_FEATURES, step 2
    config.BLIND_STEP2_FEATURES): a listed feature is taken regardless of its
    per-feature config.BLIND_*_FEATURE toggle (those toggles only gate the
    include=None default/blacklist path). drop_flags removes names either way.
    The only hard veto is a missing data dependency -- offbeam / trigtaildev /
    beamdev need beam_ref. Excluded features are never computed.

    When `X_sg` and `reac_sg` are both provided, feature computation uses the
    SG-filtered trace and its reaction strip (matching the C++ SG pass in
    ClusterVarHists). The `both` mask is always derived from raw ADC and is
    unchanged. The beam reference is always raw."""
    # Use SG trace/reac when both are provided; otherwise fall back to raw.
    Xt = X_sg if X_sg is not None and reac_sg is not None else X
    ract = reac_sg if reac_sg is not None else reac
    Xd = Xt.astype(np.float64)
    cols, names = [], []
    _toggle = {
        "pileup": config.BLIND_PILEUP_FEATURE,
        "noise": config.BLIND_NOISE_FEATURE,
        "offbeam": config.BLIND_OFFBEAM_FEATURE,
        "beamgate": config.BLIND_BEAM_GATE_FEATURE,
        "prebeam": config.BLIND_PREBEAM_FEATURE,
        "near_mult": config.BLIND_NEAR_MULT_FEATURE,
        "beamdev": config.BLIND_BEAMDEV_FEATURE,
    }
    _needs_beam = ("offbeam", "beamdev")

    def want(name):
        if name in drop_flags:
            return False
        if name in _needs_beam and beam_ref is None:
            return False  # hard data dependency, can't compute
        if include is not None:
            return name in include  # whitelist is authoritative
        return _toggle.get(name, True)  # blacklist path honors the toggles

    if want("energy"):
        cols.append(Xd.sum(axis=1) - Xd.shape[1])
        names.append("energy")
    if want("peak3"):
        cols.append(_peak3(Xt, ract))
        names.append("peak3")
    if want("plateau"):
        cols.append(_plateau_excess(Xt, ract))
        names.append("plateau")
    if want("tail"):
        cols.append(_tail_value(Xt))
        names.append("tail")
    if want("mult"):
        cols.append(both[:, 1:17].sum(axis=1).astype(np.float64))
        names.append("mult")
    if want("trig"):  # 1 if the event has a reaction trigger (ract >= 0)
        cols.append((np.asarray(ract) >= 0).astype(np.float64))
        names.append("trig")
    if want("pileup"):
        cols.append(_is_pileup(Xt).astype(np.float64))
        names.append("pileup")
    if want("noise"):
        cols.append(_is_noise(Xt).astype(np.float64))
        names.append("noise")
    if want("offbeam"):
        cols.append(_is_offbeam(Xt, beam_ref).astype(np.float64))
        names.append("offbeam")
    if want("beamgate"):
        cols.append(_beam_gate_value(Xt))
        names.append("beamgate")
    if want("prebeam"):
        cols.append(_pre_trigger_beam(Xt, ract))
        names.append("prebeam")
    if want("near_mult"):
        cols.append(_near_reaction_both(both, ract).astype(np.float64))
        names.append("near_mult")
    if want("beamdev"):
        cols.append(_beam_deviation(Xt, beam_ref))
        names.append("beamdev")
    return np.column_stack(cols), names


def _print_cluster_means(means, labels, names, indent="    "):
    """Per-cluster mean of every named feature (flag means read as the
    fraction of the cluster carrying that flag)."""
    for j in range(means.shape[0]):
        vals = ", ".join(f"{names[c]}={float(means[j, c]):.2f}"
                         for c in range(len(names)))
        print(
            f"{indent}cluster {j}: {int((labels == j).sum())} events; {vals}")


def _print_importance(feats, labels, names):
    """Per-feature clustering importance = the fraction of each (standardized)
    feature's variance that is BETWEEN clusters (ANOVA eta^2): 1 = clusters
    fully separated along that feature, 0 = identical. The unsupervised analog
    of feature importance; marginal, so it can miss a feature that only
    matters through correlations the full-cov GMM uses."""
    F = feats.astype(np.float64)
    Z = (F - F.mean(axis=0)) / (F.std(axis=0) + 1.0e-9)
    imp = np.zeros(Z.shape[1])
    for lab in np.unique(labels):
        m = labels == lab
        imp += float(m.mean()) * Z[m].mean(axis=0)**2
    print("  feature importance (between-cluster var fraction): " +
          " ".join(f"{n}={v:.2f}" for n, v in zip(names, imp)))


_AXIS_TITLE = {
    "energy": "total energy #minus beam  #Sigma(#DeltaE#minus1) [a.u.]",
    "peak3": "trigger 3-sum  #Sigma_{trig#pm1}#DeltaE [a.u.]",
    "plateau":
    "post-trigger plateau  #Sigma_{trig+1..+POST}(#DeltaE#minus1) [a.u.]",
    "tail": "raw end strip  #DeltaE(s17) [a.u.]",
    "mult": "both-channel multiplicity (strips 1-16)",
    "trigtaildev":
    "trigger+tail dev  |#DeltaE#minusbeam|_{trig} + _{s17} [a.u.]",
    "beamdev": "back-half beam dev  RMS(#DeltaE#minusbeam)_{8-17} [a.u.]",
}

# Curated continuous feature pairs for the 2-D density hists (drawn for each
# pair whose BOTH names are present in the slice's feature matrix). Ordered so
# the (a,n) signature pair (plateau up, end-strip collapse) comes first.
_HIST2D_PAIRS = (
    ("plateau", "tail"),
    ("energy", "plateau"),
    ("peak3", "tail"),
    ("trigtaildev", "beamdev"),
    ("plateau", "beamdev"),
    ("energy", "tail"),
    ("mult", "plateau"),
    ("beamdev", "mult"),  # the continuous step-1 pair (flags aside)
)


def _draw_feature_hist2d(feats, names, labels, k, tag, subdir):
    """Per-cluster 2-D feature-density hists (one TH2 BOX per cluster, coloured
    by cluster index like the trace overlays; noise label -1 in grey) for each
    curated continuous pair present in `names`. Skipped when config.BLIND_FEAT2D
    is off. Lands in a feat2d/ subdir of the stage's plot dir."""
    if not config.BLIND_FEAT2D:
        return
    idx = {n: i for i, n in enumerate(names)}
    sub = f"{subdir}/feat2d"
    for xn, yn in _HIST2D_PAIRS:
        if xn not in idx or yn not in idx:
            continue
        xc = feats[:, idx[xn]].astype(np.float64)
        yc = feats[:, idx[yn]].astype(np.float64)
        pts = []
        noise = labels == -1
        if noise.any():
            pts.append(
                (f"noise (n={int(noise.sum())})", xc[noise], yc[noise], -1))
        for j in range(k):
            m = labels == j
            if not m.any():
                continue
            pts.append((f"cluster {j} (n={int(m.sum())})", xc[m], yc[m], j))
        if pts:
            mlplots.feature_hist2d(pts,
                                   f"blind_{tag}_h2d_{xn}_{yn}",
                                   _AXIS_TITLE.get(xn, xn),
                                   _AXIS_TITLE.get(yn, yn),
                                   subdir=sub)


def step1_reaction_vs_background(X,
                                 both,
                                 reac,
                                 beam_ref=None,
                                 subdir=config.PLOT_SUBDIR,
                                 tag="step1",
                                 X_sg=None,
                                 reac_sg=None):
    """Step 1: ONE global split on config.BLIND_STEP1_FEATURES into
    config.BLIND_STEP1_K classes (None = auto-k). KEEP a cluster if the mean of
    ANY config.BLIND_STEP1_KEEP_ALWAYS feature is > 0.5 (hard keep, EXEMPT from
    the veto), OR if ANY config.BLIND_STEP1_KEEP_IF feature mean is > 0.5 AND
    EVERY
    config.BLIND_STEP1_DROP_IF feature mean is < 0.5 (conditional keep + junk
    veto). Flag-based -- no `energy` assumed. Returns the keep boolean mask.

    When `X_sg` and `reac_sg` are provided, feature computation uses the
    SG-filtered trace and its reaction strip."""
    feats, names = _cluster_matrix(X,
                                   both,
                                   reac,
                                   beam_ref=beam_ref,
                                   include=config.BLIND_STEP1_FEATURES,
                                   X_sg=X_sg,
                                   reac_sg=reac_sg)
    labels, means, k, _ = cluster_auto(feats,
                                       k=config.BLIND_STEP1_K,
                                       tag=tag,
                                       subdir=subdir,
                                       backend=config.BLIND_STEP1_BACKEND)

    def _majority(flags):  # clusters where the mean of ANY listed flag > 0.5
        out = np.zeros(k, dtype=bool)
        for f in flags:
            if f in names:
                out |= means[:, names.index(f)] > 0.5
        return out

    hard = _majority(config.BLIND_STEP1_KEEP_ALWAYS)  # exempt from the veto
    present = [f for f in config.BLIND_STEP1_KEEP_IF if f in names]
    cond = _majority(config.BLIND_STEP1_KEEP_IF) if present \
        else np.ones(k, dtype=bool)
    for f in config.BLIND_STEP1_DROP_IF:  # veto drops junk-dominated clusters
        if f in names:
            cond &= means[:, names.index(f)] < 0.5
    keep_ok = hard | cond
    keep_c = np.nonzero(keep_ok)[0]
    if keep_c.size == 0:
        print("  WARNING: step 1 keep rule matched no cluster; keeping all")
        keep_c = np.arange(k)
    keep = np.isin(labels, keep_c)
    print(f"step 1: k={k} ({', '.join(names)}) -> keep {int(keep.sum())} in "
          f"clusters {sorted(int(c) for c in keep_c)} (always "
          f"{config.BLIND_STEP1_KEEP_ALWAYS}, keep-if "
          f"{config.BLIND_STEP1_KEEP_IF}, veto {config.BLIND_STEP1_DROP_IF}), "
          f"background {int((~keep).sum())}")
    _print_cluster_means(means, labels, names, indent="  ")
    _print_importance(feats, labels, names)
    _draw_means_only(X, labels, k, "step1", beam_ref, subdir)
    # Sample noise traces (label -1), drawn same style as cluster overlay.
    _draw_noise_traces(X, labels, k, "step1", beam_ref, subdir)
    _draw_feature_hist2d(feats, names, labels, k, "step1", subdir)
    return keep


def cluster_per_reaction_strip(X,
                               both,
                               reac,
                               beam_ref,
                               strips,
                               tag,
                               drop_flags=("trig", "pileup", "offbeam"),
                               subdir=config.PLOT_SUBDIR,
                               force_k=config.BLIND_STEP2_K,
                               backend=config.BLIND_STEP2_BACKEND,
                               include=None,
                               X_sg=None,
                               reac_sg=None):
    """Partition X by the (no-fallback) leading-edge reaction strip `reac`
    (-1 = no trigger, dropped here since the slices are 3/4/5), then cluster
    WITHIN each strip in `strips`. `include`, when given (the dual / per-strip /
    beamnorm pipelines pass config.BLIND_STEP2_FEATURES), is the authoritative
    feature whitelist; when None the feature set is the blacklist path
    (everything minus drop_flags, gated by the BLIND_*_FEATURE toggles -- what
    blind_combined uses to also group junk). force_k forces the cluster count
    (None = auto-k).

    The per-cluster mean +- sigma overlay of ALL clusters is drawn; then, when
    config.BLIND_TEMPLATE_PRUNE, EVERY cluster is template-pruned to its own
    self-consistent shape (the (a,n) matching -- pruned members get label -1)
    and a `_final` overlay is drawn. Per-strip working plots (initial means +
    the template-prune rounds) go into a per-strip subdir (<subdir>/reac<strip>);
    only the `_final` overlays stay in the main <subdir> (next to step 1). The
    sliding `plateau` needs POST strips after the trigger, so a strip whose
    window would run past s17 is SKIPPED. Shared by blind_an (dual-step) and
    blind_combined (single-step).

    When `X_sg` and `reac_sg` are provided, feature computation uses the
    SG-filtered trace and its reaction strip. Partitioning by strip still uses
    the raw `reac` (plotting also uses raw `X`)."""
    last_strip = (0 if config.INCLUDE_GUARD_STRIPS else 1) + X.shape[1] - 1
    for strip in strips:
        at = reac == strip
        n_at = int(at.sum())
        print(f"reaction strip {strip}: {n_at} events")
        if n_at < 10:
            print("  too few; skipping")
            continue
        # Standalone mean of ALL events triggered at this strip (no clustering),
        # to the main dir -- the raw triggered shape per strip at a glance.
        _draw_trig_mean(X[at], f"{tag}_reac{strip}", beam_ref, subdir)
        if strip + config.BLIND_PLATEAU_POST > last_strip:
            print(
                f"  sliding plateau needs {config.BLIND_PLATEAU_POST} strips "
                f"past the trigger (> s{last_strip}); skipping")
            continue
        feats, names = _cluster_matrix(
            X[at],
            both[at],
            reac[at],
            drop_flags=drop_flags,
            beam_ref=beam_ref,
            include=include,
            X_sg=X_sg[at] if X_sg is not None else None,
            reac_sg=reac_sg[at] if reac_sg is not None else None)
        strip_subdir = f"{subdir}/reac{strip}"
        # Step 2 isolates the RARE (a,n); GMM noise tagging would clip the
        # lowest-density tail -- which can BE the (a,n) -- so it is forced off
        # (noise_pctl=None). The template prune / co-assign + an_select do the
        # (a,n) isolation instead.
        labels, means, k, bics = cluster_auto(feats,
                                              k=force_k,
                                              tag=tag,
                                              subdir=strip_subdir,
                                              backend=backend,
                                              noise_pctl=None)
        mode = "auto-k" if force_k is None else f"forced k={force_k}"
        bic_str = [round(b) if math.isfinite(b) else "nan" for b in bics]
        print(f"  {mode} on ({', '.join(names)}) -> k={k} (BIC {bic_str})")
        _print_cluster_means(means, labels, names)
        _print_importance(feats, labels, names)
        X_at = X[at]
        # Per-strip working plots (initial means + every template-prune round)
        # go into a per-strip subdir; only the `_final` overlay lands in the
        # main dir alongside step 1.
        # Combined goes many-cluster (auto-k); drop the sigma shading at k>=4
        # so the overlapping bands do not wash the canvas out.
        bands = not (tag == "combined" and k >= 4)
        _draw_means_only(X_at,
                         labels,
                         k,
                         f"{tag}_reac{strip}",
                         beam_ref,
                         strip_subdir,
                         bands=bands)
        _draw_feature_hist2d(feats, names, labels, k, f"{tag}_reac{strip}",
                             strip_subdir)
        if not config.BLIND_TEMPLATE_PRUNE:
            continue
        # Co-assign by nearest cluster mean (shape residual), or fall back to
        # the old per-cluster template_prune loop.  Pruned members get label -1
        # and drop from the _final view.
        final = labels.copy()
        if config.BLIND_TEMPLATE_CO_ASSIGN:
            final = _co_assign(X_at, final, k, tag, beam_ref, strip_subdir,
                               strip)
        else:
            for j in range(k):
                sel = labels == j
                if int(sel.sum()) < 10:
                    continue
                keep_j = template_prune(X_at[sel], f"{tag}_reac{strip}_c{j}",
                                        beam_ref, strip_subdir)
                idx = np.nonzero(sel)[0]
                final[idx[~keep_j]] = -1
        _draw_means_only(X_at,
                         final,
                         k,
                         f"{tag}_reac{strip}_final",
                         beam_ref,
                         subdir,
                         bands=bands)
        # Noise traces for the final per-strip clustering.
        _draw_noise_traces(X_at, final, k, f"{tag}_reac{strip}", beam_ref,
                           subdir)


def _draw_trig_mean(X_at, name, beam_ref, subdir):
    """One standalone plot of the MEAN trace over all events triggered at this
    strip (no clustering), with a +-1 sigma band and the beam reference dashed
    (blind_<name>_trigmean). The raw per-strip triggered shape on its own."""
    t = X_at.astype(np.float64)
    if t.shape[0] == 0:
        return
    groups = [(f"triggered (n={t.shape[0]})", t.mean(axis=0), t.std(axis=0), 1)
              ]
    mlplots.overlay_cluster_means(groups,
                                  f"blind_{name}_trigmean",
                                  name,
                                  subdir=subdir,
                                  reference=beam_ref)


def _draw_savgol_sample(X, X_sg, name, subdir):
    """Draw a sample of raw vs SG-filtered traces for comparison with C++.

    Each sample event gets its own canvas: the raw trace (dashed) and the
    SG-filtered trace (solid) overlaid, plus the beam reference (dotted).
    A fixed number of events are drawn (config.BLIND_OVERLAY_N), selected
    randomly with a fixed seed. The plot lands in the given <subdir> under
    the name blind_<name>_savgol_sample_<idx>."""
    R = mlplots._root()
    long_w = config.block_widths()[0]
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    strip_x = np.arange(first, first + long_w, dtype=np.float64)
    n = X.shape[0]
    if n == 0:
        return
    rng = np.random.default_rng(config.SEED)
    n_draw = 5
    order = np.sort(rng.choice(n, n_draw, replace=False))
    X_raw = X[order].astype(np.float64)
    X_sg_arr = X_sg[order].astype(np.float64)
    for i in range(n_draw):
        c, frame, _ = mlplots._frame_and_canvas(R, f"savgol_sample_{i}",
                                                "#DeltaE [a.u.]")
        # Raw trace (dashed).
        g_raw = mlplots._graph(R, strip_x, X_raw[i, :long_w])
        g_raw.SetLineColor(R.kGray + 2)
        g_raw.SetLineWidth(2 * R.PlottingUtils.GetLineWidth())
        g_raw.SetLineStyle(2)
        g_raw.Draw("L SAME")
        # SG-filtered trace (solid).
        g_sg = mlplots._graph(R, strip_x, X_sg_arr[i, :long_w])
        g_sg.SetLineColor(R.kAzure + 2)
        g_sg.SetLineWidth(2 * R.PlottingUtils.GetLineWidth())
        g_sg.Draw("L SAME")
        # Beam reference (dotted).
        if frame.GetYaxis().GetTitleOffset() > 0:
            pass  # frame already set up
        leg = R.PlottingUtils.AddLegend()
        leg.AddEntry(g_raw, "raw", "l")
        leg.AddEntry(g_sg, "savgol", "l")
        leg.Draw()
        label = R.PlottingUtils.AddText(f"event {order[i]}", 0.42, 0.85)
        label.Draw()
        R.PlottingUtils.SaveFigure(c, f"blind_{name}_savgol_sample_{i}",
                                   subdir, R.PlotSaveOptions.kLINEAR)
    print(f"  SG sample traces: drew {n_draw} of {n} events "
          f"(blind_{name}_savgol_sample_0..{n_draw - 1})")


def _draw_means_only(X_at, labels, k, name, beam_ref, subdir, bands=True):
    """The reduced plot set: ONLY the per-cluster mean +- sigma overlay of all
    clusters (blind_<name>_means), coloured by cluster index. No scatters,
    trace overlays or (a,n) plots. bands=False drops the sigma shading (mean
    lines only) for many-cluster plots that the bands would wash out."""
    groups = []
    for j in range(k):
        tj = X_at[labels == j].astype(np.float64)
        if tj.shape[0] == 0:
            continue
        groups.append((f"cluster {j} (n={tj.shape[0]})", tj.mean(axis=0),
                       tj.std(axis=0), j))
    if groups:
        mlplots.overlay_cluster_means(groups,
                                      f"blind_{name}_means",
                                      name,
                                      subdir=subdir,
                                      reference=beam_ref,
                                      bands=bands)


def _draw_noise_traces(X_at, labels, k, tag, beam_ref, subdir):
    """Draw sample traces of noise events (label -1) in the same style as
    draw_cluster_overlay. Uses a grey line so noise traces are visually
    distinct from cluster traces."""
    noise_sel = labels == -1
    if not noise_sel.any():
        return
    X_noise = X_at[noise_sel].astype(np.float64)
    n_avail = X_noise.shape[0]
    rng = np.random.default_rng(config.SEED)
    n_draw = min(config.BLIND_OVERLAY_N, n_avail)
    order = np.sort(rng.choice(n_avail, n_draw, replace=False))
    traces = [(f"noise (n={n_avail})", X_noise[order], 2)]
    if traces:
        print(f"  {tag} noise trace overlay: "
              f"{n_draw} of {n_avail} noise events")
        mlplots.overlay_cluster_traces(traces,
                                       f"blind_{tag}_noise_traces",
                                       f"{tag} noise traces",
                                       subdir=subdir,
                                       reference=beam_ref)


def _draw_cluster_plots(feats, X_at, labels, k, name, beam_ref, subdir):
    """The standard per-strip cluster views: the (energy, plateau) and
    (plateau, tail) scatters plus the per-cluster mean-trace overlay. Labels
    of -1 are excluded (pruned), so this draws both the initial clusters and
    the post-prune `_final` clusters identically."""
    _cluster_scatter(
        feats[:, _COL_ENERGY], feats[:, _COL_PLATEAU], labels, k,
        f"blind_{name}_energy_plateau",
        "total energy #minus beam  #Sigma(#DeltaE#minus1) [a.u.]",
        "post-trigger plateau  #Sigma_{trig+1..+POST}(#DeltaE#minus1) [a.u.]",
        subdir)
    _cluster_scatter(
        feats[:, _COL_PLATEAU], feats[:, _COL_TAIL], labels, k,
        f"blind_{name}_plateau_tail",
        "post-trigger plateau  #Sigma_{trig+1..+POST}(#DeltaE#minus1) [a.u.]",
        "end strip  #DeltaE(s17) [a.u.]", subdir)
    draw_cluster_overlay(X_at, labels, k, name, beam_ref, subdir)


def _first_valley(values, nbins=120, smooth=3):
    """Threshold at the VALLEY between the dominant (largest, near-0) residual
    peak and the NEXT peak to its right -- the boundary between the clean (a,n)
    lump and the first contamination bump. Robust to a long sparse tail, which
    a global Otsu split chases instead (it lands past every peak). Returns the
    data max (cut nothing) when there is no second peak (unimodal)."""
    v = np.asarray(values, dtype=np.float64)
    lo, hi = float(v.min()), float(v.max())
    if not hi > lo:
        return hi
    hist, edges = np.histogram(v, bins=nbins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])
    h = hist.astype(np.float64)
    if smooth > 1:
        h = np.convolve(h, np.ones(smooth) / smooth, mode="same")
    n = h.shape[0]
    p1 = int(np.argmax(h))  # dominant peak (largest; near 0)
    i = p1 + 1
    while i < n and h[i] <= h[i - 1]:  # descend off the main peak
        i += 1
    if i >= n:  # monotonic fall -> no second peak
        return hi
    while i + 1 < n and h[i + 1] >= h[i]:  # climb to the next peak
        i += 1
    valley = p1 + int(np.argmin(h[p1:i + 1]))  # lowest bin between the peaks
    return float(centers[valley])


def _cut_threshold(rk):
    """Upper residual cut + a short description, per config.BLIND_TEMPLATE_CUT:
    'valley' = the first valley between the clean lump and the next bump
    (single pass); 'mad' = median + N*MAD (robust sigma-clip, iterated)."""
    if config.BLIND_TEMPLATE_CUT == "valley":
        return _first_valley(rk), "first valley"
    med = float(np.median(rk))
    mad = float(np.median(np.abs(rk - med))) * 1.4826
    return (med + config.BLIND_TEMPLATE_NMAD * mad,
            f"med {med:.3f} + {config.BLIND_TEMPLATE_NMAD:g}*MAD {mad:.3f}")


def template_prune(X_an, tag, beam_ref, subdir):
    """Prune the off-shape tail of a cluster by template match (run on EVERY
    step-2 cluster; for the (a,n) cluster the tail is pileup). The template is
    the max-normalized survivor MEAN; each member is max-normalized too (peak
    -> 1, amplitude divided out so only SHAPE is compared) and scored by its L2
    residual to the template. The cut is config.BLIND_TEMPLATE_CUT: 'valley'
    cuts ONCE at the Otsu valley of the residual histogram (re-averaging then
    re-cutting an already-unimodal distribution would just nibble); 'mad'
    iterates median + N*MAD, re-averaging survivors, until a round prunes
    nothing or BLIND_TEMPLATE_MAX_ROUNDS. Each round draws the residual
    histogram (with the cut line) and a kept-vs-pruned mean overlay. Returns
    the boolean keep mask over X_an (the survivors)."""
    valley = config.BLIND_TEMPLATE_CUT == "valley"
    keep = np.ones(X_an.shape[0], dtype=bool)
    converged = False
    r = 0
    while r < config.BLIND_TEMPLATE_MAX_ROUNDS:
        r += 1
        template = X_an[keep].astype(np.float64).mean(axis=0)
        resid = _shape_residual(X_an, template)
        rk = resid[keep]
        thresh, desc = _cut_threshold(rk)
        new_keep = keep & (resid <= thresh)
        n_pruned = int(keep.sum() - new_keep.sum())
        print(
            f"  template prune {tag} round {r}: thresh={thresh:.3f} ({desc});"
            f" kept {int(new_keep.sum())} of {int(keep.sum())} "
            f"(pruned {n_pruned})")
        mlplots.value_hist(rk,
                           f"blind_{tag}_template_resid_r{r}",
                           "shape residual to (a,n) template [a.u.]",
                           subdir=subdir,
                           vline=thresh,
                           label=f"{tag} round {r}")
        groups = []
        pruned = keep & ~new_keep
        if new_keep.any(
        ):  # green = kept (palette[3]), red = pruned (palette[2])
            kt = X_an[new_keep].astype(np.float64)
            groups.append((f"kept (n={int(new_keep.sum())})", kt.mean(axis=0),
                           kt.std(axis=0), 3))
        if pruned.any():
            pt = X_an[pruned].astype(np.float64)
            groups.append((f"pruned (n={int(pruned.sum())})", pt.mean(axis=0),
                           pt.std(axis=0), 2))
        if groups:
            mlplots.overlay_cluster_means(groups,
                                          f"blind_{tag}_template_r{r}",
                                          f"{tag} template round {r}",
                                          subdir=subdir,
                                          reference=beam_ref)
        keep = new_keep
        if valley or n_pruned == 0:  # valley = single pass; mad = until stable
            converged = valley or n_pruned == 0
            break
    if converged:
        print(
            f"  template prune {tag}: {'valley cut' if valley else 'converged'}"
            f" after {r} round(s), {int(keep.sum())} (a,n) survivors")
    else:
        print(f"  template prune {tag}: hit max "
              f"{config.BLIND_TEMPLATE_MAX_ROUNDS} rounds, "
              f"{int(keep.sum())} (a,n) survivors")
    return keep


def _co_assign(X_all, labels, k, tag, beam_ref, strip_subdir, strip):
    """Reassign events by nearest cluster mean (shape residual) across ALL
    clusters simultaneously, iterating until assignments stabilise.

    This replaces the old per-cluster ``template_prune`` loop: instead of each
    cluster pruning in isolation (discarding outliers without re-assignment),
    every event is scored against *every* cluster's max-normalized mean and
    moved to the nearest one.  Events whose minimum residual exceeds the cut
    are given label ``-1`` (pruned) permanently.

    Returns the updated label array."""
    valley = config.BLIND_TEMPLATE_CUT == "valley"
    converged = False
    r = 0
    while r < config.BLIND_CO_MAX_ROUNDS:
        r += 1
        # Build per-cluster templates from current survivors.
        templates = []
        for j in range(k):
            tj = X_all[labels == j]
            if tj.shape[0] == 0:
                # Empty cluster -- pad with zeros (same width) so the residual
                # matrix stays aligned; it will never win a nearest-neighbour.
                templates.append(np.zeros(X_all.shape[1], dtype=np.float64))
            else:
                templates.append(tj.astype(np.float64).mean(axis=0))
        tstack = np.asarray(templates)  # (k, n_cols)
        # Compute shape residuals for every event against every cluster mean.
        if tstack.size == 0:
            # No clusters (k=0) -- skip co-assign, keep labels as-is.
            print(f"  co-assign {tag}: k=0, skipping co-assign")
            break
        resid_all = _shape_residual(X_all, tstack)  # (n_events, k)
        best_j = resid_all.argmin(axis=1)  # (n_events,)
        min_resid = resid_all.min(axis=1)  # (n_events,)
        # Apply the cut globally on the *minimum* residual.  If all assigned
        # events were pruned in the previous round, just reassign them to
        # nearest mean and move on.
        assign_prev = min_resid[labels >= 0]
        if len(assign_prev) == 0:
            thresh = float(min_resid.max())
            desc = "all-pruned prev round; use max"
        else:
            thresh, desc = _cut_threshold(assign_prev)
        new_labels = np.where(min_resid <= thresh, best_j, -1).astype(int)
        n_changed = int((new_labels != labels).sum())
        print(
            f"  co-assign {tag} round {r}: thresh={thresh:.3f} ({desc});"
            f" {int((new_labels >= 0).sum())} assigned, {int((new_labels == -1).sum())} pruned "
            f"(changed {n_changed} labels)")
        # Residual histogram (global, over all currently-assigned events).
        assign_mask = new_labels >= 0
        if assign_mask.any():
            mlplots.value_hist(min_resid[assign_mask],
                               f"blind_{tag}_reac{strip}_co_resid_r{r}",
                               "min shape residual across clusters [a.u.]",
                               subdir=strip_subdir,
                               vline=thresh,
                               label=f"{tag} reac{strip} co-assign round {r}")
            # Per-cluster mean overlay: colour assigned survivors by new label,
            # pruned ones in red for visibility.
            groups = []
            for j in range(k):
                sj = new_labels == j
                if not sj.any():
                    continue
                tj = X_all[sj].astype(np.float64)
                groups.append((f"cluster {j} (n={int(sj.sum())})",
                               tj.mean(axis=0), tj.std(axis=0), j))
            # Pruned events.
            pm = new_labels == -1
            if pm.any():
                pt = X_all[pm].astype(np.float64)
                groups.append((f"pruned (n={int(pm.sum())})", pt.mean(axis=0),
                               pt.std(axis=0), 2))
            mlplots.overlay_cluster_means(
                groups,
                f"blind_{tag}_reac{strip}_co_r{r}",
                f"{tag} reac{strip} co-assign round {r}",
                subdir=strip_subdir,
                reference=beam_ref)
        labels = new_labels
        if valley or n_changed == 0:
            converged = valley or n_changed == 0
            break
    if converged:
        print(
            f"  co-assign {tag} reac{strip}: "
            f"{'valley cut' if valley else 'converged'} after {r} round(s), "
            f"{int((labels >= 0).sum())} assigned, {int((labels == -1).sum())} pruned"
        )
    else:
        print(f"  co-assign {tag} reac{strip}: hit max "
              f"{config.BLIND_CO_MAX_ROUNDS} rounds, "
              f"{int((labels >= 0).sum())} assigned, "
              f"{int((labels == -1).sum())} pruned")
    return labels


def _max_norm(traces):
    """Each trace divided by its own max so the peak sits at 1 -- SHAPE only,
    amplitude removed. Rows with a non-positive max are left unscaled."""
    t = traces.astype(np.float64)
    m = t.max(axis=1, keepdims=True)
    m[m <= 0.0] = 1.0
    return t / m


def _shape_residual(traces, template):
    """L2 distance between each max-normalized trace and the max-normalized
    template -- the shape mismatch with amplitude divided out.

    When *template* is 1-D ``(ncols,)`` returns residuals of shape ``(n_events,)``.
    When *template* is 2-D ``(k, ncols)`` returns a distance matrix of shape
    ``(n_events, k)`` (each event scored against every cluster mean)."""
    traces = np.asarray(traces, dtype=np.float64)
    if traces.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    tn = np.asarray(template, dtype=np.float64)
    if tn.size == 0:
        # No templates (k=0) -- return large residuals so nothing wins.
        return np.full(traces.shape[0], 1e9, dtype=np.float64)
    tnorm = _max_norm(traces)  # (n, m) or (1, m)
    if tn.ndim == 1:  # single template row
        tn = tn / (tn.max() if tn.max() > 0.0 else 1.0)  # (m,)
        d = tnorm - tn[np.newaxis, :]  # (n, m)
        return np.sqrt((d * d).sum(axis=1))  # (n,)
    # 2-D stack of templates (k, m) -- normalize each row independently.
    rows_max = tn.max(axis=1, keepdims=True)  # (k, 1)
    rows_max[rows_max <= 0.0] = 1.0
    tn_normed = tn / rows_max  # (k, m)
    d = tnorm[:, np.newaxis, :] - tn_normed[np.newaxis, :, :]  # (n, k, m)
    return np.sqrt((d * d).sum(axis=2))  # (n, k)


def _cluster_scatter(x,
                     y,
                     labels,
                     k,
                     name,
                     x_title,
                     y_title,
                     subdir=config.PLOT_SUBDIR):
    """Scatter of (x, y) colored by cluster label (one series per cluster).
    The 4th tuple element is the cluster index `j`, which keys the colour so a
    cluster keeps it across every plot."""
    mlplots.cluster_scatter([(f"cluster {j} (n={int((labels == j).sum())})",
                              x[labels == j], y[labels == j], j)
                             for j in range(k)],
                            name,
                            x_title,
                            y_title,
                            subdir=subdir)


def draw_cluster_overlay(X,
                         labels,
                         k,
                         tag,
                         beam_ref=None,
                         subdir=config.PLOT_SUBDIR):
    """Per-cluster trace plots (StripSumScatter DrawRegionTraces style, fixed
    0.8-1.3). Always: each cluster's MEAN trace (over ALL its events) with a
    +-1 sigma band, one bold line per cluster, beam reference dashed
    (blind_<tag>_means). When k <= config.BLIND_OVERLAY_MAX_K, ALSO the
    StripSumScatter-style overlay of up to config.BLIND_OVERLAY_N individual
    traces per cluster (blind_<tag>_traces). Each group carries its cluster
    index `j` as the colour key, so a cluster keeps its colour everywhere."""
    means = []
    for j in range(k):
        tj = X[labels == j].astype(np.float64)
        if tj.shape[0] == 0:
            continue
        means.append((f"cluster {j} (n={tj.shape[0]})", tj.mean(axis=0),
                      tj.std(axis=0), j))
    if not means:
        return
    print(f"  {tag} cluster means (+-1 sigma): {len(means)} clusters")
    mlplots.overlay_cluster_means(means,
                                  f"blind_{tag}_means",
                                  tag,
                                  subdir=subdir,
                                  reference=beam_ref)
    if k <= config.BLIND_OVERLAY_MAX_K or config.BLIND_OVERLAY_MAX_K == -1:
        rng = np.random.default_rng(config.SEED)
        traces = []
        for j in range(k):
            tj = X[labels == j]
            n_avail = tj.shape[0]
            if n_avail == 0:
                continue
            n_draw = min(config.BLIND_OVERLAY_N, n_avail)
            order = np.sort(rng.choice(n_avail, n_draw, replace=False))
            traces.append((f"cluster {j} (n={n_avail})", tj[order], j))
        if traces:
            print(f"  {tag} trace overlay (k<={config.BLIND_OVERLAY_MAX_K}): "
                  f"up to {config.BLIND_OVERLAY_N} traces each")
            mlplots.overlay_cluster_traces(traces,
                                           f"blind_{tag}_traces",
                                           tag,
                                           subdir=subdir,
                                           reference=beam_ref)


def main():
    random.seed(config.SEED)
    np.random.seed(config.SEED)
    print(f"music-ml: dataset {config.DATASET}, BLIND (a,n) dual-step "
          "(global k=2 -> per reaction-strip)")
    X, both = step0_reservoir()
    beam_ref, beam_sigma = _beam_reference(X)
    if beam_ref is None:
        raise RuntimeError("no pure-beam events; cannot set the "
                           "reaction-strip threshold")
    reac = _reaction_strip_all(X, beam_ref, beam_sigma)
    X, both, reac = _apply_prebeam_cut(X, both, reac)
    X_sg, reac_sg = None, None
    if config.BLIND_SAVITZKY_GOLAY:
        X_sg = data.savgol_filter_trace(X)
        reac_sg = _reaction_strip_all(X_sg, beam_ref, beam_sigma)
        print(f"  SG filter: reac detected {int((reac_sg >= 0).sum())} "
              f"triggered (vs {int((reac >= 0).sum())} raw)")
    subdir = f"{config.PLOT_SUBDIR}/an"
    if config.BLIND_SAVITZKY_GOLAY:
        _draw_savgol_sample(X, X_sg, "an_savgol", subdir)
    # Step 1: drop flat background / noise with one global k=2 cut.
    keep = step1_reaction_vs_background(X,
                                        both,
                                        reac,
                                        beam_ref,
                                        subdir,
                                        X_sg=X_sg,
                                        reac_sg=reac_sg)
    # Step 2: per reaction-strip clustering on the kept (real-deposit) half.
    cluster_per_reaction_strip(
        X[keep],
        both[keep],
        reac[keep],
        beam_ref,
        config.BLIND_COMBINED_STRIPS,
        "dual",
        subdir=subdir,
        include=config.BLIND_STEP2_FEATURES,
        X_sg=X_sg[keep] if X_sg is not None else None,
        reac_sg=reac_sg[keep] if reac_sg is not None else None)
    print("done")


if __name__ == "__main__":
    main()
