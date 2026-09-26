"""Classify the experimental events and write what the comparison needs.

The reservoir is EVERY calibrated event that passes the strip-0 vs strip-1
beam-entrance gate, and nothing else: no noise cut, no pileup cut, no
off-beam cut, no pre-beam cut. The gate only asks that the particle entered
as beam, which every candidate reaction must be true of. Pruning further
would answer the question for the model.

That gate is blind_an's own -- _fit_beam_ellipse / _beam_gate at
config.BLIND_BEAM_NSIGMA -- called directly rather than reimplemented, so the
two pipelines cannot drift apart at the one cut they share.

Classification is strip-agnostic. Every gated event is classified and
written out with its SeedTs and its four posteriors; binning by reaction
strip is a separate, later step on that table.

It is deliberately NOT done here with blind_an.reaction_strip. That CFD
gates on a peak excess above BLIND_REAC_ONSET_NSIGMA times the beam RMS,
which on run 16 is 5 * 0.0346 = 0.173 -- while the (a,n) plateau in
region_traces_reac<r>.png sits about 0.10-0.15 above beam. The threshold is
above the signal, so it keeps the upper tail (pileup and large
fluctuations) and drops most real reactions: 297 of 205614 events on one
subfile. Selecting on it would decide the answer before the model saw
anything.

What comes out is the per-event table the two published-comparison checks
both start from: SeedTs is what the upstream comparison joins on, so "did
this pick out the same events the published build tagged" is a join rather
than an argument, and a per-strip yield for the cross section is a groupby
on whatever strip assignment turns out to be trustworthy.

Run inside the dataset dev shell, from python/:

    python vlm_apply.py
"""

import numpy as np

import blind_an
import config
import data
import vlm


def beam_gated_reservoir(max_files=None):
    """Calibrated traces passing the (s0,s1) beam gate, and nothing else.

    Mirrors the first half of blind_an.step0_reservoir -- same loader, same
    unit-gain trace view, same fitted ellipse -- and then stops, where
    step0_reservoir would go on to apply the noise / pileup / offbeam cuts.
    Returns (X, both, seed_ts) with seed_ts in the same row order.
    """
    print("vlm: reservoir = calibrated events, (s0,s1) beam gate ONLY")
    totals, both, seed_ts = data.load_experimental_totals(
        max_files=config.VLM_MAX_FILES if max_files is None else max_files,
        max_events_per_file=config.BLIND_MAX_EVENTS_PER_FILE,
        with_seed_ts=True)
    long_w = config.block_widths()[0]
    unit = np.ones(totals.shape[1], dtype=np.float64)
    X, kept = data.max_normalize(totals, unit)
    X = X[:, :long_w]
    both, seed_ts = both[kept], seed_ts[kept]
    blind_an._set_gate_fit(X)
    keep = blind_an._beam_gate(X, ((0, 1), ))
    print(f"  beam-0/1 gate (fitted s0,s1 ellipse, "
          f"{config.BLIND_BEAM_NSIGMA}-sigma): kept {int(keep.sum())} "
          f"of {X.shape[0]}")
    X, both, seed_ts = X[keep], both[keep], seed_ts[keep]
    print(f"  reservoir: {X.shape[0]} events ({long_w} long-side strips) "
          "-- no other cut applied")
    return X, both, seed_ts


def main():
    np.random.seed(config.SEED)
    print(f"music-ml: dataset {config.DATASET}, VLM applied to data "
          f"({config.VLM_MODEL})")

    X, _both, seed_ts = beam_gated_reservoir()
    beam_ref, _sigma = blind_an._beam_reference(X)
    if beam_ref is None:
        raise RuntimeError("no pure-beam events; cannot draw the in-image "
                           "beam reference the prompt refers to")

    for seed in config.VLM_SEEDS:
        order, dy, dhalf = vlm.seed_perturbation(seed)
        letters = ", ".join(f"{config.VLM_CLASS_TOKENS[i]}="
                            f"{config.VLM_CLASSES[c]}"
                            for i, c in enumerate(order))
        print(f"\nseed {seed}: {letters}; dy={dy:+.2f}px dhalf={dhalf:+.2f}px")
        clf = vlm.EventClassifier(letter_order=order)
        probs, labels = clf.classify(X, beam_ref, dy=dy, dhalf=dhalf)
        del clf

        for i, name in enumerate(config.VLM_CLASSES):
            n = int((labels == i).sum())
            print(f"  {name:>6}: {n:>9}  ({100.0 * n / labels.size:.2f}%)")
        # The number that matters: the tag is a threshold on p(an), and the
        # saved posteriors let it be moved afterwards without re-running.
        pan = vlm.an_probability(probs)
        tagged = vlm.is_an(probs)
        print(f"  tagged (a,n) at p>={config.VLM_AN_THRESHOLD}: "
              f"{int(tagged.sum())} of {tagged.size} "
              f"({100.0 * tagged.mean():.3f}%)")
        for thr in (0.3, 0.5, 0.7, 0.9):
            print(f"    p>={thr}: {int((pan >= thr).sum())}")

        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        out = config.CACHE_DIR / f"vlm_events_seed{seed}.npz"
        np.savez(out,
                 seed_ts=seed_ts,
                 probs=probs.astype(np.float32),
                 labels=labels.astype(np.int8),
                 classes=np.array(config.VLM_CLASSES))
        print(f"  per-event table -> {out}")

    print("\nEach table is (seed_ts, four posteriors, label) per gated event.")
    print("Join on seed_ts against the published tags for the event-identity")
    print("check; group by whatever strip assignment is trusted for the yield.")
    print("done")


if __name__ == "__main__":
    main()
