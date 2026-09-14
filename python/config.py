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

# MUSIC_DATASET_DIR is <git root>/analysis/<dataset>, so two levels up is the
# repository root -- where the gitignored secrets sit.
REPO_ROOT = DATASET_DIR.parent.parent

CACHE_DIR = RESULTS_DIR / "ml_cache"
PLOTS_DIR = RESULTS_DIR / "plots"
ROOT_FILES_DIR = RESULTS_DIR / "root_files"
SIM_ROOT_FILES_DIR = RESULTS_DIR / "sim_root_files"
PLOT_SUBDIR = "ml"

N_STRIPS = 18

# Trace view (mirrors EnergyView::Decode); match C++ Constants.
IGNORE_SHORT_STRIPS = False  # split strips 1-16: long end only
INCLUDE_UNSEGMENTED_STRIPS = True  # keep the unsegmented strips 0/17
INCLUDE_DERIVATIVE = False
SEED = 42


def block_widths():
    """Long-side block width: 18 columns, or 16 with strips 0/17 dropped."""
    long_w = N_STRIPS if INCLUDE_UNSEGMENTED_STRIPS else N_STRIPS - 2
    return (long_w, )


BLIND_MAX_FILES = None  # None = every events file
BLIND_MAX_EVENTS_PER_FILE = None  # None = every event per file

# Each junk metric is both a hard cut (BLIND_REJECT_*, applied in step 0) and a
# 0/1 clustering feature (BLIND_*_FEATURE). Beam ellipse is fit, not assumed.
BLIND_GATE_BEAM_01 = True  # hard (s0,s1)-ellipse beam gate in step 0
BLIND_BEAM_NSIGMA = 3.0  # ellipse half-width (gate AND beamgate feature)
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
BLIND_PREBEAM_MIN_STRIP = 2  # floor: pre-trigger pair = strips 1 and 0

# Savitzky-Golay smoothing: 5-point, cubic, edge-renormalized (matches
# StripSumScatter::SavitzkyGolay); clustering uses the SG trace, beam ref stays raw.
BLIND_SAVITZKY_GOLAY = True

# Shape/topology clustering features.
BLIND_PLATEAU_POST = 3  # `plateau`: Sigma(dE-1) over reac+1..reac+POST
BLIND_MULT_THRESH = 0.0  # `mult`: count strips 1-16 where both ends fire
BLIND_NEAR_MULT_FEATURE = True  # `near_mult`: both-fire within +-1 of trigger
BLIND_BEAMDEV_FEATURE = False  # `beamdev`: RMS(dE-beam) over strips 8-17

# Step 1 (real-vs-background split): keep a cluster if ANY KEEP_ALWAYS mean
# > 0.5, OR ANY KEEP_IF mean > 0.5 AND EVERY DROP_IF mean < 0.5.
BLIND_STEP1_FEATURES = ("beamgate", "offbeam", "pileup", "beamdev", "mult")
BLIND_STEP1_K = 3  # None = auto-k (GMM only)
BLIND_STEP1_KEEP_ALWAYS = ()
BLIND_STEP1_KEEP_IF = ("beamgate", )
BLIND_STEP1_DROP_IF = ("offbeam", "pileup")

# Step 2 (per-strip shape clustering): curated whitelist, AUTHORITATIVE (see
# _cluster_matrix). (a,n) = plateau-up + end-strip collapse -> shape axes only.
BLIND_STEP2_FEATURES = ("plateau", "tail", "beamdev", "near_mult")

# Cluster count per stage: None = auto-k (k minimizing GMM BIC), int = force.
# BLIND_MAX_K bounds the auto-k search (-1 = unbounded 1..19).
BLIND_MAX_K = -1
BLIND_STEP2_K = 4
BLIND_COMBINED_K = 4

# Clustering backend (cluster_auto): "gmm" (sklearn GaussianMixture, BIC model
# selection, native predict) or "none"; optional noise tail (BLIND_GMM_NOISE_PCTL).
BLIND_NOISE_CLUSTERING = "gmm"
BLIND_STEP1_BACKEND = None  # override BLIND_NOISE_CLUSTERING for step 1
BLIND_STEP2_BACKEND = None  # override BLIND_NOISE_CLUSTERING for step 2

# GMM noise tail: bottom-percentile of per-point logL -> noise (-1). None = off.
BLIND_GMM_NOISE_PCTL = 5.0

# Reaction-strip onset (constant-fraction discriminator).
BLIND_REAC_ONSET_NSIGMA = 5.0  # trigger-exists gate: peak excess > N*beam RMS
BLIND_REAC_ONSET_FRAC = 0.30  # onset = first strip reaching this frac of peak
BLIND_COMBINED_STRIPS = tuple(range(2, 14))  # candidate reaction strips

# Step-2 shape cleanup: co-assign reassigns events to the nearest cluster
# mean-trace (shape residual), pruning outliers to -1; template_prune is the fallback.
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

# --- VLM event classification (vlm.py, vlm_eval.py) ---
# Zero-shot per-event classification: each event's trace is rasterized to a
# small RGB image in memory and pushed through a vision-language model, and
# the class posteriors are read from the logits at the last prompt position.
# The model NEVER generates -- see vlm.EventClassifier for why that matters.

# The scaling ladder, smallest first. The question this run asks is not
# only "does it work" but "how much model does it take" -- if E2B matches
# 12B on the sim confusion matrix, the task is being solved by the vision
# encoder and shallow pattern-matching and scale buys nothing; if it is bad
# at EVERY rung, the rendering is wrong and no rung will save it.
#
# The E-series counts are effective, not total: E2B and E4B carry per-layer
# embeddings, so E4B is 4.5B effective but ~8B of weights on disk and in
# VRAM. Size the batch off the total, not the effective count.
#
# The 12B id is the one rung not confirmed against the Hub in this session.
VLM_MODEL_LADDER = (
    "google/gemma-4-E2B-it",
    "google/gemma-4-E4B-it",
    "google/gemma-4-12B-it",
)
VLM_MODEL = VLM_MODEL_LADDER[2]

# E2B result, so the ladder is not being climbed on a hunch: on 200 events
# matched at plateau +0.056 and differing ONLY in end level (0.79 vs 1.02 --
# the (a,n) collapse), the separation in p(an) was +0.022, +0.022, -0.001,
# -0.045 across four render/budget configurations. No trend, no
# discrimination. The patch grid was 30x24 even at the lowest budget (1.33
# columns per strip), so resolution was never the constraint.

# 11.95B in bf16 is ~24 GB and will not fit a 24 GB card alongside
# activations. "8bit" is ~12 GB and the safer choice for a discrimination
# that was already marginal; "4bit" is ~6.5 GB if that still does not fit.
# None loads bf16, for a card with room.
VLM_LOAD_IN = "8bit"
# Gemma is a gated repo. vlm.py reads the token from $HF_TOKEN when that is
# set, and otherwise from this file -- one line, gitignored.
VLM_HF_TOKEN_FILE = REPO_ROOT / "hftoken"
VLM_DTYPE = "bfloat16"
VLM_DEVICE = "cuda"
VLM_ATTN = "sdpa"
# EVENTS per forward pass -- with a reference figure that is two images per
# event, and Gemma 4's vision tower allocates a one_hot over pre-pool patch
# positions, so memory grows with the square of the patch grid rather than
# with the token count. classify() halves this on OOM and carries on, so it
# is a starting point rather than a value to get right.
VLM_BATCH = 32

# Subfiles the VLM run reads. Its OWN knob, not BLIND_MAX_FILES: capping the
# VLM run must not quietly change the blind pipeline's reservoir. One subfile
# is ~205k gated events -- a real run producing a real table, in about an
# hour. None means all 48, which is order a day per seed.
VLM_MAX_FILES = 1

# Visual token budget. A Gemma 4 image processor emits a FIXED number of soft
# tokens per image, selectable from {70, 140, 280, 560, 1120} (default 280);
# Google's guidance is that classification wants the low end. Prefill cost is
# linear in this, so it is the one knob that sets the run time. vlm.py
# MEASURES what the processor actually produced and refuses to run quietly if
# the request did not take.
# 70, and the reference figure IS legible at it: p(an) on a random gated draw
# was 0.06 with the figure at this budget, against 0.73 without the figure.
# The tiers are 70/140/280/560/1120, but memory in Gemma 4's vision tower
# grows with the SQUARE of the patch grid -- 280 needs ~3 GB of one_hot per
# image, which on a 24 GB card with two images per event only runs at batch
# 1. Move up a tier only with evidence that 70 is costing accuracy.
VLM_TOKEN_BUDGET = 70
# How the budget is handed to the processor. transformers 5.5 wants
# per-call processor arguments inside a `processor_kwargs` dict rather than
# as bare **kwargs, and the argument's NAME is not documented on the model
# card -- setting attributes on the image processor does nothing (the
# self-test reports "set via no known attr" and measures 256). Run
# `vlm_selftest.py 2`, which dumps the image processor's config keys and its
# __call__ signature, and put the right name here. None leaves the default.
VLM_BUDGET_KWARG = None

# Rasterization (vlm.render_traces): numpy only, a whole batch at a time. No
# matplotlib (~10-30 ms/figure would dominate everything else) and no disk.
VLM_IMG_H = 112
VLM_IMG_W = 112
VLM_LINE_HALFWIDTH = 1.0  # rows drawn either side of the polyline
# Absolute dE window in calibrated a.u. (beam == 1). These are the values
# StripSumScatter::DrawRegionTraces frames region_traces_reac<r> with -- the
# figure that shows the three classes separating -- so the model is handed
# the same view a person reads the classification off. Both features that
# distinguish (a,n) live in it: the raised plateau near 1.1 across the middle
# strips, and the collapse to 0.7-0.9 over the last few. A tighter window
# clips the collapse off the canvas and leaves (a,n) and (a,a') looking
# alike. Every event is drawn on THIS window, never rescaled per trace.
VLM_DE_MIN = 0.6
VLM_DE_MAX = 1.6
VLM_DRAW_BEAM_REF = True  # beam reference drawn under the trace
VLM_SPLIT_LR = False  # draw L and R as two curves instead of their sum

VLM_BG_RGB = (255, 255, 255)
VLM_TRACE_RGB = (0, 0, 0)
VLM_BEAM_RGB = (200, 60, 60)
VLM_SECOND_RGB = (40, 90, 200)  # the R curve when VLM_SPLIT_LR is on

# Classes, the single token that stands for each, and the one-line
# description that goes in the prompt. Every token MUST tokenize to exactly
# one token (vlm.py asserts this): the answer is a softmax over these four
# token ids alone, so the model never has to obey a formatting instruction
# and a small model that would otherwise ramble still scores cleanly.
#
# The letters are assigned to classes at run time rather than baked into the
# prompt, because which letter means which class is one of the arbitrary
# choices the seed ensemble below perturbs.
VLM_CLASSES = ("an", "aa", "beam", "other")
VLM_CLASS_TOKENS = ("A", "B", "C", "D")
VLM_CLASS_DESC = {
    "an":
    "(alpha,n) -- steps UP at the reaction strip to about 1.1, holds that "
    "raised plateau across the middle strips, then collapses to 0.7-0.9 over "
    "the last few strips. BOTH features together are the signature.",
    "aa":
    "(alpha,alpha\') -- makes a small but SUSTAINED step up at the reaction "
    "strip, holding a level slightly above 1.0 for several strips, then "
    "returning to 1.0. No collapse at the end. Choose this over beam only "
    "when a step is actually visible, not for ordinary noise.",
    "beam":
    "beam -- no reaction. The trace wanders around 1.0 with small "
    "strip-to-strip noise and never makes a sustained step to a new level. "
    "MOST events are this one.",
    "other":
    "none of the above -- noise, pileup, or an incomplete track.",
}

# A worked example, shown before the event: the region-trace figure the C++
# already draws, where the three classes are overlaid in colour at one
# reaction strip. Telling the model which strip it is for and then asking
# about other strips is the generalization test. Set to None to prompt from
# the text descriptions alone -- worth comparing, because the reference costs
# a second image on EVERY request (see vlm.PromptBuilder).
# Measured, not assumed: with this figure p(an) on a random gated draw (~95%
# unreacted beam) is 0.06; without it, 0.73. The example is what tells the
# model that most events are beam -- the text descriptions alone do not.
VLM_REFERENCE_IMAGE = (PLOTS_DIR / "strip_sum_scatter" /
                       "region_traces_reac2.png")
VLM_REFERENCE_STRIP = 2
VLM_REFERENCE_DESC = (
    "The first image is a reference, not the event to classify. It overlays "
    "many measured events at reaction strip {strip}, on the same axes as the "
    "second image: grey = beam, blue = (alpha,alpha\'), red = (alpha,n). "
    "Learn the three shapes from it.")

VLM_PROMPT_HEADER = (
    "These are events from a MUSIC active-target ionization chamber: energy "
    "loss per strip along the beam axis, strip 0 at the left, plotted on a "
    "fixed vertical scale of 0.6 to 1.6 where unreacted beam sits at 1.0.")
VLM_PROMPT_QUESTION = (
    "Classify the single event in the last image. Its red line marks the "
    "unreacted beam level.")
VLM_PROMPT_FOOTER = "Answer with one letter."

# --- seed ensemble -------------------------------------------------------
# The forward pass is deterministic and the label is an argmax, so re-running
# with a different seed changes NOTHING on its own -- the spread would be
# exactly zero and the systematic meaningless. A seed here therefore perturbs
# the choices that are genuinely arbitrary in the method, and the spread in
# the resulting cross section across seeds is what gets quoted:
#
#   - which letter stands for which class (VLMs carry documented label-token
#     and option-order biases, so this is a real axis, not noise);
#   - sub-pixel rasterization jitter, which asks whether the answer survives
#     choices made with no physics behind them.
#
# Each seed writes its own tag-efficiency store, so the C++ cross section is
# run once per seed by pointing CROSS_SECTION_CONFIG.TAG_EFFICIENCY_FILE at
# each in turn, and the spread is taken over its per-strip output.
VLM_SEEDS = (42, 43, 44)
VLM_SEED_PERMUTE_LABELS = True
VLM_SEED_JITTER = True
VLM_JITTER_DY = 1.0  # max |row offset| in pixels
VLM_JITTER_DHALF = 0.5  # max change in the line half-width

# Where the records land. `{seed}` is substituted per seed; the basename is
# what CROSS_SECTION_CONFIG.TAG_EFFICIENCY_FILE has to be set to, and living
# under its own name leaves the bootstrap tag_efficiency.root untouched.
VLM_CHANNEL = "an"
VLM_TAG_STORE_FMT = "tag_efficiency_vlm_seed{seed}.root"
VLM_METHOD_FMT = "vlm:{model}:seed{seed}"

# The efficiency written into the tag store. There is no labelled truth set
# here by design -- the check is against the PUBLISHED analysis, not against
# labels invented locally -- so nothing measures this yet and 1.0 means the
# cross section that comes out is an uncorrected count. Good enough to ask
# whether the shape and scale reproduce; not a publishable absolute.
VLM_ASSUMED_EFF = 1.0

# The decision is (a,n) vs NOT, on p(an) alone. (a,a') is not a result we
# want, only a place for (a,a')-like events to go so they are not pushed
# into the (a,n) column -- which is why the prompt keeps four classes while
# nothing downstream distinguishes the other three. A threshold on p(an) is
# used rather than the argmax: an event can be the argmax of a flat
# distribution at p(an)=0.3, and the whole reason for reading posteriors
# instead of generating a label is that this line can be swept afterwards
# from the saved table without re-running the model.
VLM_AN_THRESHOLD = 0.5

# Cache the shared prompt prefix's KV once and reuse it across every event.
# With the reference-then-text-then-event ordering, everything but the
# event's own ~62 image tokens is identical on every request, so this is most
# of the forward pass. It is verified against the uncached path at startup
# and switched off automatically if the logits disagree -- an optimization
# that changes answers is a bug, not an optimization.
# OFF: it does not currently reproduce the uncached logits (0.81 against a
# 0.05 tolerance), and the verification refuses to use it. The likely cause
# is that Gemma 4 interleaves local sliding-window attention with global
# attention, and those layers need the model's own hybrid cache class rather
# than a plain cache continued across two forwards -- so this is a different
# cache object, not a slicing fix. Left in place because the ~2x it would buy
# only matters for a full 48-subfile run.
VLM_PREFIX_CACHE = True

# Print a render / processor / forward breakdown every N batches; 0 = off.
# Worth leaving on: the batch size going 16 -> 32 moved throughput 31 -> 32
# events/s, which says the bottleneck is not the GPU, and guessing which
# serial stage it is has been wrong more than once. Timing the forward
# requires a CUDA synchronise, so this costs a little of what it measures.
VLM_TIMING_EVERY = 20
# NOT an accept/reject threshold -- the cache is accepted or rejected on
# whether any real event changes its (a,n) tag, which is the only thing that
# reaches a yield. This is the level above which the measured |dp(an)| is
# called out as a systematic worth comparing against the seed spread. Three
# earlier attempts to gate on a constant (raw logits, then all four
# posteriors, then this) all rejected a cache that flipped nothing.
VLM_PREFIX_CACHE_TOL = 0.02

# Self-test slice: subfiles to read, and the reaction strips to show a row
# of on the contact sheet. One file is 200k events, plenty to find triggered
# events at every strip without waiting on the full run.
VLM_SELFTEST_FILES = 1
VLM_SHEET_PER_ROW = 8
VLM_SHEET_ROWS = 4

VLM_PLOT_SUBDIR = "vlm"
