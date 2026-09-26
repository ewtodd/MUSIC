"""Relate TabFM scores on experimental tagged events to trace topology."""

import numpy as np

import config
import vlm_domain


def main():
    score_path = config.CACHE_DIR / "tabfm_compute_regions.npz"
    input_path = config.CACHE_DIR / "tabfm_compute_regions_input.npz"
    scores = np.load(score_path)
    inputs = np.load(input_path)
    pan = scores["probabilities"][:, list(config.VLM_CLASSES).index("an")]
    reaction_strip = inputs["reac"].astype(np.int32)
    tagged = reaction_strip >= 0
    features = vlm_domain.trace_features(inputs["X"][tagged],
                                         reaction_strip[tagged])
    names = ("pre", "plateau", "tail", "rise", "collapse")
    print("Pearson correlation with p(an), existing tagged events")
    for i, name in enumerate(names):
        correlation = np.corrcoef(features[:, i], pan[tagged])[0, 1]
        print(f"  {name:<8}: {correlation:+.3f}")

    tagged_indices = np.flatnonzero(tagged)
    order = tagged_indices[np.argsort(pan[tagged])]
    for label, indices in (("lowest", order[:20]), ("highest", order[-20:])):
        topology = vlm_domain.trace_features(inputs["X"][indices],
                                             reaction_strip[indices])
        print(f"{label} 20 tagged scores: p(an) "
              f"{np.median(pan[indices]):.4f}")
        print("  " + vlm_domain.summarize(label, topology))
        print("  reaction strips "
              + " ".join(str(int(value)) for value in reaction_strip[indices]))


if __name__ == "__main__":
    main()
