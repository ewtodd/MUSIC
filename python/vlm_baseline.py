"""Train an XGBoost baseline on the VLM simulator split.

This deliberately imports vlm_finetune's loader so Gemma and XGBoost see the
same normalized traces, balanced classes, random sample, and held-out reaction
strip. It is the cheap test of whether a large vision-language model adds
anything to an 18-number classification problem.
"""

import argparse

import numpy as np

import config
import vlm_finetune


def confusion_matrix(truth, prediction):
    n_classes = len(config.VLM_CLASSES)
    matrix = np.zeros((n_classes, n_classes), dtype=np.int64)
    np.add.at(matrix, (truth, prediction), 1)
    return matrix


def print_report(matrix, probabilities, truth):
    print("  rows = truth, columns = prediction")
    print("           " + " ".join(f"{name:>7}" for name in config.VLM_CLASSES))
    for i, name in enumerate(config.VLM_CLASSES):
        count = matrix[i].sum()
        recall = float(matrix[i, i] / count) if count else float("nan")
        cells = " ".join(f"{n:7d}" for n in matrix[i])
        print(f"  {name:>7} {cells}   recall {recall:.3f}")
    threshold, an_recall, background_fpr, beam_fpr = \
        vlm_finetune.threshold_metrics(probabilities, truth)
    print(f"  selected p(an)>={threshold:.4f}: an recall {an_recall:.3f}, "
          f"background false-positive {background_fpr:.3f}, "
          f"beam false-positive {beam_fpr:.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-per-class", type=int,
                        default=config.VLM_FINETUNE_MAX_PER_CLASS)
    parser.add_argument("--holdout-strip", type=int, default=None)
    args = parser.parse_args()
    if args.holdout_strip is not None:
        config.VLM_FINETUNE_VAL_STRIPS = (args.holdout_strip, )

    X, y, holdout, _beam_ref, _strips = vlm_finetune.load_examples(
        args.max_per_class)
    train_rows = np.flatnonzero(~holdout)
    validation_rows = np.flatnonzero(holdout)

    from xgboost import XGBClassifier
    classifier = XGBClassifier(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="multi:softprob",
        num_class=len(config.VLM_CLASSES),
        eval_metric="mlogloss",
        random_state=config.VLM_FINETUNE_SEED,
        n_jobs=-1)
    classifier.fit(X[train_rows], y[train_rows])
    probabilities = classifier.predict_proba(X[validation_rows])
    prediction = probabilities.argmax(axis=1)
    matrix = confusion_matrix(y[validation_rows], prediction)
    print_report(matrix, probabilities, y[validation_rows])


if __name__ == "__main__":
    main()
