"""Zero-shot TabFM evaluation on the simulator leave-one-strip split.

`fit` supplies labelled simulator rows as in-context examples; it does not
optimize TabFM's pretrained weights. This is therefore the tabular analogue of
the original no-fine-tuning Gemma test, evaluated on exactly the same traces as
the supervised baselines.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

import config
import vlm_baseline
import vlm_finetune


def load_experimental_subfile():
    """Load a pre-exported beam-gated subfile without importing ROOT here."""
    path = config.CACHE_DIR / "tabfm_experimental_input.npz"
    if not path.is_file():
        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        code = (
            "import numpy as np; import blind_an, vlm_apply; "
            "X,both,ts=vlm_apply.beam_gated_reservoir(max_files=1); "
            "beam,_=blind_an._beam_reference(X); live=beam>0; "
            f"np.savez({str(path)!r}, X=X[:,live]/beam[None,live], "
            "seed_ts=ts, live=live)"
        )
        subprocess.run([sys.executable, "-c", code], check=True,
                       cwd=Path(__file__).parent)
    arrays = np.load(path)
    return arrays["X"], arrays["seed_ts"]


def load_compute_regions_reservoir():
    """Load the small experimental tag reservoir before TabFM imports."""
    path = config.CACHE_DIR / "tabfm_compute_regions_input.npz"
    if not path.is_file():
        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        code = (
            "import numpy as np; import vlm_reservoir; "
            "X,m=vlm_reservoir.load_reservoir(); "
            "flat=m['beam_flat']; beam=X[flat].mean(axis=0); live=beam>0; "
            f"np.savez({str(path)!r}, X=X[:,live]/beam[None,live], "
            "seed_ts=m['seed_ts'], reac=m['reac'], beam_flat=flat)"
        )
        subprocess.run([sys.executable, "-c", code], check=True,
                       cwd=Path(__file__).parent)
    arrays = np.load(path)
    return (arrays["X"], arrays["seed_ts"], arrays["reac"],
            arrays["beam_flat"])


def predict_in_chunks(classifier, X, chunk, output, seed_ts, extra=None):
    """Resumable bounded prediction, writing after every independent chunk."""
    completed = 0
    probabilities = []
    if output.is_file():
        saved = np.load(output)
        old = saved["probabilities"]
        if old.shape[0] <= X.shape[0]:
            completed = old.shape[0]
            probabilities.append(old)
            print(f"  resuming TabFM at {completed}/{X.shape[0]}")
    for start in range(completed, X.shape[0], chunk):
        stop = min(start + chunk, X.shape[0])
        probabilities.append(classifier.predict_proba(X[start:stop]))
        current = np.concatenate(probabilities)
        payload = {
            "seed_ts": seed_ts[:current.shape[0]],
            "probabilities": current.astype(np.float32),
            "classes": np.array(config.VLM_CLASSES),
        }
        if extra:
            for name, values in extra.items():
                payload[name] = values[:current.shape[0]]
        np.savez(output, **payload)
        print(f"  TabFM {stop}/{X.shape[0]}")
    return np.concatenate(probabilities)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-per-class", type=int, default=500)
    parser.add_argument("--holdout-strip", type=int, default=8)
    parser.add_argument("--context-rows", type=int, default=100)
    parser.add_argument("--estimators", type=int, default=8)
    parser.add_argument("--batch", type=int, default=32)
    parser.add_argument("--sim-negative-scale", type=float, default=1.0,
                        help="scale simulated below-beam deviations toward beam")
    parser.add_argument("--experimental", action="store_true",
                        help="also classify a beam-gated experimental sample")
    parser.add_argument("--reservoir", action="store_true",
                        help="classify the compute-regions tagged/beam-flat cache")
    parser.add_argument("--experimental-max", type=int, default=1000,
                        help="maximum experimental rows; raise only after a safe probe")
    parser.add_argument("--experimental-chunk", type=int, default=32,
                        help="rows per independent TabFM predict_proba call")
    parser.add_argument("--experimental-offset", type=int, default=0,
                        help="first experimental row to classify")
    args = parser.parse_args()
    config.VLM_FINETUNE_VAL_STRIPS = (args.holdout_strip, )

    X, y, holdout, beam_ref, _strips = vlm_finetune.load_examples(
        args.max_per_class)
    if args.sim_negative_scale != 1.0:
        delta = X - beam_ref[np.newaxis, :]
        X = beam_ref[np.newaxis, :] + np.where(
            delta < 0.0, delta * args.sim_negative_scale, delta)
        print(f"scaled simulated below-beam deviations by "
              f"{args.sim_negative_scale:.3f}")
    train_rows = np.flatnonzero(~holdout)
    validation_rows = np.flatnonzero(holdout)

    experimental = None
    experimental_seed_ts = None
    reservoir = None
    if args.experimental:
        experimental, experimental_seed_ts = load_experimental_subfile()
        if experimental.shape[1] != X.shape[1]:
            # Strip 17 is disabled in this dataset and absent from experimental
            # traces. Remove the same dead feature from simulator context.
            live = np.any(experimental != 0.0, axis=0)
            if live.size == X.shape[1] and int(live.sum()) == experimental.shape[1]:
                X = X[:, live]
            elif X.shape[1] - experimental.shape[1] == 1:
                X = X[:, :-1]
            else:
                raise RuntimeError(f"simulator has {X.shape[1]} features but "
                                   f"experiment has {experimental.shape[1]}")
    if args.reservoir:
        reservoir = load_compute_regions_reservoir()
        if reservoir[0].shape[1] != X.shape[1]:
            if X.shape[1] - reservoir[0].shape[1] == 1:
                X = X[:, :-1]
            else:
                raise RuntimeError(f"simulator has {X.shape[1]} features but "
                                   f"reservoir has {reservoir[0].shape[1]}")

    from tabfm import TabFMClassifier
    from tabfm import tabfm_v1_0_0_pytorch as tabfm_v1_0_0
    model = tabfm_v1_0_0.load(model_type="classification")
    classifier = TabFMClassifier(
        model=model,
        max_num_rows=args.context_rows,
        n_estimators=args.estimators,
        batch_size=args.batch,
        random_state=config.VLM_FINETUNE_SEED)
    classifier.fit(X[train_rows], y[train_rows])
    probabilities = classifier.predict_proba(X[validation_rows])
    prediction = probabilities.argmax(axis=1)
    matrix = vlm_baseline.confusion_matrix(y[validation_rows], prediction)
    print(f"TabFM zero-shot: {train_rows.size} context-pool rows, "
          f"{args.context_rows} rows per context, {args.estimators} estimators")
    vlm_baseline.print_report(matrix, probabilities, y[validation_rows])

    if args.experimental:
        if args.experimental_max < 1 or args.experimental_chunk < 1:
            parser.error("experimental limits must be positive")
        if args.experimental_offset < 0:
            parser.error("--experimental-offset cannot be negative")
        begin = args.experimental_offset
        end = min(begin + args.experimental_max, experimental.shape[0])
        experimental = experimental[begin:end]
        experimental_seed_ts = experimental_seed_ts[begin:end]
        output = config.CACHE_DIR / f"tabfm_experimental_{begin}_{end}.npz"
        experimental_probabilities = predict_in_chunks(
            classifier, experimental, args.experimental_chunk, output,
            experimental_seed_ts)
        pan = experimental_probabilities[:, list(config.VLM_CLASSES).index("an")]
        threshold, _an_recall, background_fpr, _beam_fpr = \
            vlm_finetune.threshold_metrics(probabilities, y[validation_rows])
        selected = pan >= threshold
        print(f"experimental subfile: {selected.sum()} of {selected.size} events "
              f"({selected.mean():.4%}) at simulated p(an)>={threshold:.4f}; "
              f"validation background leakage {background_fpr:.3f}")
        np.savez(output, seed_ts=experimental_seed_ts,
                 probabilities=experimental_probabilities.astype(np.float32),
                 classes=np.array(config.VLM_CLASSES), threshold=threshold)
        print(f"  per-event table -> {output}")

    if args.reservoir:
        reservoir_X, reservoir_seed_ts, reac, beam_flat = reservoir
        scale_tag = str(args.sim_negative_scale).replace(".", "p")
        output = config.CACHE_DIR / f"tabfm_compute_regions_neg{scale_tag}.npz"
        reservoir_probabilities = predict_in_chunks(
            classifier, reservoir_X, args.experimental_chunk, output,
            reservoir_seed_ts, {"reac": reac, "beam_flat": beam_flat})
        pan = reservoir_probabilities[:, list(config.VLM_CLASSES).index("an")]
        threshold, _an_recall, background_fpr, _beam_fpr = \
            vlm_finetune.threshold_metrics(probabilities, y[validation_rows])
        tagged = reac >= 0
        selected = pan >= threshold
        print(f"compute-regions at p(an)>={threshold:.4f}: tagged acceptance "
              f"{selected[tagged].mean():.3f}, beam-flat false-positive "
              f"{selected[beam_flat].mean():.3f}, simulated background "
              f"leakage {background_fpr:.3f}")
        for strip in sorted(set(int(value) for value in reac[tagged])):
            rows = reac == strip
            print(f"  strip {strip}: {selected[rows].sum()}/{rows.sum()} "
                  f"({selected[rows].mean():.3f}), median p(an) "
                  f"{np.median(pan[rows]):.4f}")


if __name__ == "__main__":
    main()
