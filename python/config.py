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

SIM_ROOT_DIR = DATASET_DIR / "sim_root_files"
MODELS_DIR = RESULTS_DIR / "models"
CACHE_DIR = RESULTS_DIR / "ml_cache"
PLOTS_DIR = RESULTS_DIR / "plots"
ROOT_FILES_DIR = RESULTS_DIR / "root_files"
PLOT_SUBDIR = "ml"

N_STRIPS = 18

REACTION_STRIP_MIN = 3
REACTION_STRIP_MAX = 15

INCLUDE_SHORT_STRIPS = False
INCLUDE_GUARD_STRIPS = True
INCLUDE_DERIVATIVE = True
SUBTRACT_BEAM_LEVEL = True


def block_widths():
    """Column blocks of the totals view, in order.

    Always starts with the long-side trace (18 columns, or 16 with the
    guard strips dropped). With INCLUDE_SHORT_STRIPS a second block of the
    16 short-side ends (strips 1-16) follows as separate columns -- the
    ends are never summed, so the L/R sharing survives as a feature.
    """
    long_w = N_STRIPS if INCLUDE_GUARD_STRIPS else N_STRIPS - 2
    if INCLUDE_SHORT_STRIPS:
        return (long_w, N_STRIPS - 2)
    return (long_w, )


def n_strip_inputs():
    """Strip columns in the input vector (all blocks)."""
    return sum(block_widths())


def ae_input_dim():
    """Autoencoder input width: trace plus optional appended derivative.

    The derivative is taken within each block (w - 1 diffs per block), so
    no junk column spans the long/short boundary.
    """
    n = n_strip_inputs()
    if INCLUDE_DERIVATIVE:
        n += sum(w - 1 for w in block_widths())
    return n


SEED = 42

# autoencoder for discriminating beam vs. not beam
AE_LATENT_DIM = 5
AE_HIDDEN_DIMS = (48, 24)
AE_BATCH_SIZE = 512
AE_MAX_EPOCHS = 1000
AE_PATIENCE = 100
AE_LEARNING_RATE = 4.0e-3
AE_VAL_FRACTION = 0.5
# Anomaly threshold: this percentile of the beam validation reconstruction
# error (i.e. the accepted beam false-positive rate is 100 - percentile).
AE_THRESHOLD_PERCENTILE = 99.5

# penalize the model when a sim reacted trace reconstructs too well. Per negative trace the hinge
# max(0, log(R * <beam batch MSE>) - log(MSE_neg)) pushes its error above
# R times the beam level.
AE_OUTLIER_EXPOSURE = True
AE_OE_POPULATIONS = ("aa", )  # exposure negatives; (a,n) is easy already
AE_OE_MARGIN_RATIO = 5.0  # R: target err_neg >= R * mean beam err
AE_OE_WEIGHT = 1.0  # hinge weight relative to the beam MSE term
# Negatives split: this fraction is held out of OE training, then divided
# evenly into an early-stopping half and a final-AUC test half.
AE_OE_HOLDOUT_FRACTION = 0.5

TRACE_N_SAMPLES = 5
TRACE_STRIP = 4  # reaction strip for the sim (a,a') / (a,n) samples
TRACE_DRAW_DATA = True  # also sample experimental events

# OAT sweeps + randomized search for the feature configuration currently
# set above.
TUNE_SEEDS = (42, 123, 256)
# Reduced training budget for the study (retrain the winner with the
# production AE_MAX_EPOCHS/AE_PATIENCE afterwards).
TUNE_MAX_EPOCHS = 150
TUNE_PATIENCE = 20
# Beam traces per fit and reacted traces per population for the AUC; caps
# keep one fit ~O(1 min) so the full grid stays overnight-sized.
TUNE_TRAIN_CAP = 40000
TUNE_EVAL_CAP_PER_POP = 20000
# One-at-a-time sweep values (other params held at the AE_* values above).
TUNE_LATENT_VALUES = (1, 2, 3, 4, 6, 8)
TUNE_HIDDEN_VALUES = ((12, 6), (16, 8), (24, 12), (32, 16), (48, 24))
TUNE_LR_VALUES = (3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)
TUNE_BATCH_VALUES = (128, 256, 512, 1024)
# Outlier-exposure sweep values (only swept when AE_OUTLIER_EXPOSURE).
TUNE_OE_RATIO_VALUES = (3.0, 10.0, 30.0, 100.0)
TUNE_OE_WEIGHT_VALUES = (0.3, 1.0, 3.0, 10.0)
# Randomized-search iterations around the OAT recommendations.
TUNE_RS_ITER = 25
# KernelSHAP pass: background beam traces / reacted traces to explain, and
# how many top-ranked combos get the (expensive) SHAP treatment.
TUNE_SHAP_BG = 100
TUNE_SHAP_EXPLAIN = 200
TUNE_TOP_COMBOS = 3
