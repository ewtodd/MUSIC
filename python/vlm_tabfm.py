"""Zero-shot TabFM evaluation on the simulator leave-one-strip split.

`fit` supplies labelled simulator rows as in-context examples; it does not
optimize TabFM's pretrained weights. This is therefore the tabular analogue of
the original no-fine-tuning Gemma test, evaluated on exactly the same traces as
the supervised baselines.
"""

import argparse

import numpy as np

import config
import vlm_baseline
import vlm_finetune


def load_experimental_subfile():
    """Load a pre-exported beam-gated subfile without importing ROOT here."""
    path = config.CACHE_DIR / "tabfm_experimental_input.npz"
    if not path.is_file():
        raise FileNotFoundError(
            f"no experimental input at {path}; export it with the established "
            "experimental loader before importing TabFM. The rebased Python "
            "3.14 environment currently segfaults in that ROOT/NumPy path.")
    arrays = np.load(path)
    return arrays["X"], arrays["seed_ts"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-per-class", type=int, default=500)
    parser.add_argument("--holdout-strip", type=int, default=8)
    parser.add_argument("--context-rows", type=int, default=100)
    parser.add_argument("--estimators", type=int, default=8)
    parser.add_argument("--batch", type=int, default=32)
    parser.add_argument("--experimental", action="store_true",
                        help="also classify one beam-gated experimental subfile")
    args = parser.parse_args()
    config.VLM_FINETUNE_VAL_STRIPS = (args.holdout_strip, )

    X, y, holdout, _beam_ref, _strips = vlm_finetune.load_examples(
        args.max_per_class)
    train_rows = np.flatnonzero(~holdout)
    validation_rows = np.flatnonzero(holdout)

    experimental = None
    experimental_seed_ts = None
    if args.experimental:
        experimental, experimental_seed_ts = load_experimental_subfile()

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
        experimental_probabilities = classifier.predict_proba(experimental)
        pan = experimental_probabilities[:, list(config.VLM_CLASSES).index("an")]
        threshold, _an_recall, background_fpr, _beam_fpr = \
            vlm_finetune.threshold_metrics(probabilities, y[validation_rows])
        selected = pan >= threshold
        print(f"experimental subfile: {selected.sum()} of {selected.size} events "
              f"({selected.mean():.4%}) at simulated p(an)>={threshold:.4f}; "
              f"validation background leakage {background_fpr:.3f}")
        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        output = config.CACHE_DIR / "tabfm_experimental_subfile.npz"
        np.savez(output, seed_ts=experimental_seed_ts,
                 probabilities=experimental_probabilities.astype(np.float32),
                 classes=np.array(config.VLM_CLASSES), threshold=threshold)
        print(f"  per-event table -> {output}")


if __name__ == "__main__":
    main()
