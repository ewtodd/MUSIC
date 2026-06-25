import os
from pathlib import Path


def _require_env(name):
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(
            f"{name} is not set. Enter a dataset dev shell first "
            "(e.g. `nix develop .#87Rb`).")
    return value


DATASET = _require_env("MUSIC_DATASET")
DATASET_DIR = Path(_require_env("MUSIC_DATASET_DIR"))
RESULTS_DIR = Path(os.environ.get("MUSIC_RESULTS_DIR", str(DATASET_DIR)))

CACHE_DIR = RESULTS_DIR / "ml_cache"
PLOTS_DIR = RESULTS_DIR / "plots"
ROOT_FILES_DIR = RESULTS_DIR / "root_files"
PLOT_SUBDIR = "ml"

N_STRIPS = 18

# Trace view (mirrors EnergyView::Decode); match C++ Constants.
IGNORE_SHORT_STRIPS = False  # split strips 1-16: long end only
INCLUDE_GUARD_STRIPS = True  # keep guards 0/17
INCLUDE_DERIVATIVE = False
SEED = 42


def block_widths():
    """Long-side block width: 18 columns, or 16 with guards dropped."""
    long_w = N_STRIPS if INCLUDE_GUARD_STRIPS else N_STRIPS - 2
    return (long_w, )


BLIND_MAX_FILES = None  # None = every events file
BLIND_MAX_EVENTS_PER_FILE = None  # None = every event per file

# Each junk metric is both a hard cut (BLIND_REJECT_*, applied in step 0) and a
# 0/1 clustering feature (BLIND_*_FEATURE). Beam ellipse is fit, not assumed.
BLIND_GATE_BEAM_01 = True  # hard (s0,s1)-ellipse beam gate in step 0
BLIND_BEAM_NSIGMA = 5.0  # ellipse half-width (gate AND beamgate feature)
BLIND_BEAM_GATE_FEATURE = False  # `beamgate`: inside (s0,s1) ellipse

BLIND_REJECT_PILEUP = True  # `pileup`: >= MIN_STRIPS strips (1-16) >= THRESH
BLIND_PILEUP_THRESH = 1.3
BLIND_PILEUP_MIN_STRIPS = 1
BLIND_PILEUP_FEATURE = False

BLIND_REJECT_NOISE = True  # `noise`: >= MIN_STRIPS strips < THRESH
BLIND_NOISE_THRESH = 0.4
BLIND_NOISE_MIN_STRIPS = 1
BLIND_NOISE_FEATURE = False

# `offbeam`: >= MIN_STRIPS strips (1-16) >= DIST from beam (above OR below).
# MIN_STRIPS high so a real (a,n) -- departs over a few strips -- reads 0.
BLIND_REJECT_OFFBEAM = True
BLIND_OFFBEAM_FEATURE = False
BLIND_OFFBEAM_DIST = 0.3
BLIND_OFFBEAM_MIN_STRIPS = 4

# `prebeam`: still inside the beam blob at (reac-1, reac-2)? 1 = still beam.
BLIND_REJECT_PREBEAM = False
BLIND_PREBEAM_FEATURE = True
BLIND_PREBEAM_MIN_STRIP = 2  # floor: pre-trigger pair = strips 1 and 0 (guard)

# Savitzky-Golay smoothing: 5-point, cubic, with edge renormalization
# (matches StripSumScatter::SavitzkyGolay in C++). When enabled, clustering
# features and the reaction strip for clustering are computed from the
# SG-filtered trace rather than the raw trace. The beam reference remains
# raw (computed from raw beam-flat events).
BLIND_SAVITZKY_GOLAY = True

# Shape/topology clustering features.
BLIND_PLATEAU_POST = 3  # `plateau`: Sigma(dE-1) over reac+1..reac+POST
BLIND_MULT_THRESH = 0.0  # `mult`: count strips 1-16 where both ends fire
BLIND_NEAR_MULT_FEATURE = True  # `near_mult`: both-fire within +-1 of trigger
BLIND_BEAMDEV_FEATURE = False  # `beamdev`: RMS(dE-beam) over strips 8-17

# Step 1 (global real-vs-background split). Keep a cluster if ANY KEEP_ALWAYS
# mean > 0.5 (exempt from veto), OR ANY KEEP_IF mean > 0.5 AND EVERY DROP_IF
# mean < 0.5 (conditional keep + junk veto).
BLIND_STEP1_FEATURES = ("beamgate", "offbeam", "pileup", "beamdev", "mult")
BLIND_STEP1_K = 3  # None = auto-k (gmm only)
BLIND_STEP1_KEEP_ALWAYS = ()
BLIND_STEP1_KEEP_IF = ("beamgate", )
BLIND_STEP1_DROP_IF = ("offbeam", "pileup")

# Step 2 (per-strip shape clustering): curated feature whitelist. AUTHORITATIVE
# -- listed features are used regardless of their BLIND_*_FEATURE toggle (those
# toggles only gate the default blacklist path). None = fall back to that
# blacklist (everything minus drop_flags). The (a,n) signature is plateau-up +
# end-strip collapse, so the default is the shape axes only -- no amplitude
# (energy/peak3) or junk flags to dilute the split. blind_combined ignores this
# (keeps its all-features blacklist by design).
BLIND_STEP2_FEATURES = ("plateau", "tail", "beamdev", "near_mult")

# Cluster count per stage: None = auto-k (k minimizing GMM BIC), int = force.
# BLIND_MAX_K bounds the auto-k search (-1 = unbounded 1..19). Honored by the
# "gmm" backend only -- torch_dbscan returns its own component count.
BLIND_MAX_K = -1
BLIND_STEP2_K = 4
BLIND_COMBINED_K = None

# Clustering backend (cluster_auto):
#   "gmm"/"none"    -- sklearn GaussianMixture, BIC model selection, native
#                      predict; optional noise tail via BLIND_GMM_NOISE_PCTL.
#   "torch_dbscan"  -- GPU DBSCAN (radius-graph cores -> connected components),
#                      fit on a capped subsample, every row assigned to nearest
#                      core within eps (min dist > eps -> noise, label -1).
BLIND_NOISE_CLUSTERING = "torch_dbscan"
BLIND_NOISE_FIT_CAP = int(1e4)  # rows to FIT on; None = all

# Per-stage backend override (None = follow BLIND_NOISE_CLUSTERING). Step 2
# resolves an (a,n)/(a,a') CONTINUUM -- gmm carves it by BIC; torch_dbscan,
# needing a density gap, returns one blob (k=1).
BLIND_STEP1_BACKEND = "gmm"
BLIND_STEP2_BACKEND = "gmm"

# GMM noise tail (gmm backend only): bottom-percentile of per-point logL ->
# noise (-1). None = off.
BLIND_GMM_NOISE_PCTL = 5.0

# torch_dbscan params (z-scored space).
BLIND_DBSCAN_EPS = None  # None = k-distance elbow; float = manual
BLIND_DBSCAN_MIN_SAMPLES = 100  # min points (incl. self) for a core
BLIND_TORCH_DEVICE = "auto"  # "auto"/"cuda"/"cpu"
BLIND_TORCH_TILE = 16384  # rows per GPU distance tile

# Free GPU before torch work: unload any resident llama-swap model (it reloads
# on its next request). URL "" or BLIND_FREE_GPU = False disables.
BLIND_FREE_GPU = True
BLIND_LLAMA_SWAP_URL = "http://127.0.0.1:8080"

# Reaction-strip onset (constant-fraction discriminator).
BLIND_REAC_ONSET_NSIGMA = 3.0  # trigger-exists gate: peak excess > N*beam RMS
BLIND_REAC_ONSET_FRAC = 0.30  # onset = first strip reaching this frac of peak
BLIND_COMBINED_STRIPS = tuple(range(2, 14))  # candidate reaction strips

# Step-2 shape cleanup: co-assign reassigns events to the nearest cluster
# mean-trace (max-norm shape residual), pruning shape outliers to -1; the
# per-cluster template_prune loop is the fallback when CO_ASSIGN is off.
BLIND_TEMPLATE_PRUNE = True
BLIND_TEMPLATE_CUT = "mad"  # "valley" (first valley) | "mad" (median+N*MAD)
BLIND_TEMPLATE_NMAD = 1.5  # "mad" only
BLIND_TEMPLATE_MAX_ROUNDS = 10
BLIND_TEMPLATE_CO_ASSIGN = True
BLIND_CO_MAX_ROUNDS = 10

# Plotting.
BLIND_FIT_CAP = 50_000  # rows to FIT on (predict on all); None = all
BLIND_OVERLAY_MAX_K = 3  # also draw the N-trace overlay when k <= this
BLIND_OVERLAY_N = 40
BLIND_FEAT2D = True  # per-stage 2D feature-density hists (feat2d/ subdir)
BLIND_FEAT2D_MAXPTS = 1_000_000  # per-cluster subsample before filling
