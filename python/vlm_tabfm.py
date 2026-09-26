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


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-per-class", type=int, default=500)
    parser.add_argument("--holdout-strip", type=int, default=8)
    parser.add_argument("--context-rows", type=int, default=100)
    parser.add_argument("--estimators", type=int, default=8)
    parser.add_argument("--batch", type=int, default=32)
    parser.add_argument("--experimental", action="store_true",
                        help="also classify a beam-gated experimental sample")
    parser.add_argument("--experimental-max", type=int, default=1000,
                        help="maximum experimental rows; raise only after a safe probe")
    parser.add_argument("--experimental-chunk", type=int, default=32,
                        help="rows per independent TabFM predict_proba call")
    parser.add_argument("--experimental-offset", type=int, default=0,
                        help="first experimental row to classify")
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
        completed = 0
        chunks = []
        if output.is_file():
            saved = np.load(output)
            old_probabilities = saved["probabilities"]
            if old_probabilities.shape[0] <= experimental.shape[0]:
                completed = old_probabilities.shape[0]
                chunks.append(old_probabilities)
                print(f"  resuming TabFM experimental at {completed}/{experimental.shape[0]}")
        for start in range(completed, experimental.shape[0],
                           args.experimental_chunk):
            stop = min(start + args.experimental_chunk, experimental.shape[0])
            chunks.append(classifier.predict_proba(experimental[start:stop]))
            current = np.concatenate(chunks)
            np.savez(output, seed_ts=experimental_seed_ts[:current.shape[0]],
                     probabilities=current.astype(np.float32),
                     classes=np.array(config.VLM_CLASSES))
            print(f"  TabFM experimental {stop}/{experimental.shape[0]}")
        experimental_probabilities = np.concatenate(chunks)
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


if __name__ == "__main__":
    main()
