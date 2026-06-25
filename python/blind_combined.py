"""Single-step blind clustering on all variables combined.

A second strategy alongside blind_an.py (which it imports and reuses, never
modifies): partitions the (optionally beam-gated) reservoir by the no-fallback
leading-edge reaction strip FIRST (blind_an.reaction_strip), then clusters
WITHIN each slice in config.BLIND_COMBINED_STRIPS -- auto-k, or forced to
config.BLIND_COMBINED_K classes (its own knob, default auto-k) -- on the same
feature vector blind_an step 2 uses

    (energy above beam #Sigma(#DeltaE-1), trigger 3-sum peak3, plateau excess
     #Sigma_{2-10}(#DeltaE-1), raw end strip #DeltaE(s17), both-channel
     multiplicity, near-reaction both-channel multiplicity, trigtaildev, and
     the binary noise / pileup / offbeam flags whose config.BLIND_*_FEATURE is
     on)

Unlike blind_an's dual step there is NO step-1 split, so this drops only
`trig` (constant within a slice) and KEEPS the junk flags (noise / pileup /
offbeam, which step 1 would otherwise have grouped) so the single clustering
groups that garbage itself. Only the per-cluster mean +- sigma overlay is drawn
per slice (the reduced plot set). The reaction strip is the slice dimension, NOT
a clustering feature.

Run inside the dataset dev shell, from python/:

    python -m music_ml.blind_combined
"""

import random

import numpy as np

import config
import data
import blind_an


def main():
    random.seed(config.SEED)
    np.random.seed(config.SEED)
    print(f"music-ml: dataset {config.DATASET}, BLIND combined clustering "
          "(one step, all variables)")

    X, both = blind_an.step0_reservoir()
    beam_ref, beam_sigma = blind_an._beam_reference(X)
    if beam_ref is None:
        raise RuntimeError("no pure-beam events; cannot set the "
                           "reaction-strip threshold")
    reac = blind_an._reaction_strip_all(X, beam_ref, beam_sigma)
    X, both, reac = blind_an._apply_prebeam_cut(X, both, reac)
    X_sg, reac_sg = None, None
    if config.BLIND_SAVITZKY_GOLAY:
        X_sg = data.savgol_filter_trace(X)
        reac_sg = blind_an._reaction_strip_all(X_sg, beam_ref, beam_sigma)
        print(f"  SG filter: reac detected {int((reac_sg >= 0).sum())} "
              f"triggered (vs {int((reac >= 0).sum())} raw)")
    subdir = f"{config.PLOT_SUBDIR}/combined"
    if config.BLIND_SAVITZKY_GOLAY:
        blind_an._draw_savgol_sample(X, X_sg, "combined_savgol", subdir)
    # Single step: partition by reaction strip, cluster within each slice.
    # Drop only `trig` (constant in a slice); KEEP `noise` (no step 1 here).
    blind_an.cluster_per_reaction_strip(X, both, reac, beam_ref,
                                        config.BLIND_COMBINED_STRIPS,
                                        "combined", drop_flags=("trig",),
                                        subdir=subdir,
                                        force_k=config.BLIND_COMBINED_K,
                                        X_sg=X_sg, reac_sg=reac_sg)
    print("done")


if __name__ == "__main__":
    main()
