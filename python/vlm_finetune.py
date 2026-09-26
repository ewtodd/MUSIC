"""Fine-tune the small Gemma VLM on labelled Remix-MUSIC-Sim traces.

The simulator supplies the labels: files matching the per-class globs in
config.VLM_FINETUNE_FILE_GLOBS become (alpha,n), (alpha,alpha'), beam, and
other examples. The default `other` source is simulated (alpha,p), a real
reaction topology outside the requested classes. The MC tree's reaction_strip
makes the validation split physical:
all reactions at VLM_FINETUNE_VAL_STRIPS are held out, so validation measures
transfer to a vertex location not used for optimization. Beam events are split
deterministically by event index.

Generate high-statistics, detector-noise-matched files first, for example by
raising `run.n_events` in the controls under analysis/37Cl/sim_control/ and
running them with ~/Remix-MUSIC-Sim. The output paths must be under this
dataset's sim_root_files directory and match the configured file globs.

Run from python/ inside the dataset shell:

    python vlm_finetune.py
    python vlm_finetune.py --epochs 5 --max-per-class 50000

The output directory contains only the LoRA adapter and training metadata. Use
--adapter with vlm_finetune.py --evaluate to score it later; applying an adapter
to the production VLM path is deliberately a separate change from training.
"""

import argparse
import inspect
import json
from pathlib import Path

import numpy as np

import config
import vlm


def _sim_files():
    """Existing simulator files grouped by canonical VLM class name."""
    found = {}
    base = config.VLM_FINETUNE_SIM_DIR
    for name, globs in config.VLM_FINETUNE_FILE_GLOBS.items():
        paths = []
        for pattern in globs:
            paths.extend(base.glob(pattern))
        found[name] = sorted(set(paths))
        if not found[name]:
            raise FileNotFoundError(f"no {name} simulation files in {base} "
                                    f"matching {globs}")
    return found


def _load_sim_file(path):
    """Return raw (n, 18) MeV totals and MC reaction strips from one file."""
    import ROOT

    events = ROOT.RDataFrame("events_MeV", str(path)).AsNumpy(
        ["LeftdE", "RightdE", "Strip0dE", "Strip17dE"])
    truth = ROOT.RDataFrame("MC", str(path)).AsNumpy(["reaction_strip"])
    n = len(events["Strip0dE"])
    if len(truth["reaction_strip"]) != n:
        raise RuntimeError(f"{path}: events_MeV has {n} rows but MC has "
                           f"{len(truth['reaction_strip'])}")
    left = np.stack(events["LeftdE"]).astype(np.float32)
    right = np.stack(events["RightdE"]).astype(np.float32)
    if left.shape != (n, 16) or right.shape != (n, 16):
        raise RuntimeError(f"{path}: expected 16 segmented strips, got "
                           f"LeftdE {left.shape}, RightdE {right.shape}")
    totals = np.empty((n, config.N_STRIPS), dtype=np.float32)
    totals[:, 0] = events["Strip0dE"]
    totals[:, 1:17] = left + right
    totals[:, 17] = events["Strip17dE"]
    return totals, np.asarray(truth["reaction_strip"], dtype=np.int32)


def _sample_indices(n, maximum, rng):
    if maximum is None or n <= maximum:
        return np.arange(n, dtype=np.int64)
    return np.sort(rng.choice(n, maximum, replace=False))


def load_examples(max_per_class=None):
    """Read balanced, beam-normalized simulator examples and their split.

    The beam reference is computed exclusively from simulated unreacted-beam
    files. This preserves amplitude information in reaction traces while using
    the identical beam-equals-one convention as vlm_apply.py.
    """
    rng = np.random.default_rng(config.VLM_FINETUNE_SEED)
    files = _sim_files()
    beam_raw = []
    for path in files["beam"]:
        raw, strip = _load_sim_file(path)
        if not np.all(strip == -1):
            raise RuntimeError(f"{path}: beam training file has reacted events")
        beam_raw.append(raw)
    beam_ref = np.concatenate(beam_raw).mean(axis=0, dtype=np.float64)
    if np.any(beam_ref <= 0.0):
        raise RuntimeError("simulated beam reference has a non-positive strip")

    per_class = {}
    for name in config.VLM_CLASSES:
        rows = []
        row_strips = []
        for path in files[name]:
            raw, reaction_strip = _load_sim_file(path)
            expected = -1 if name == "beam" else None
            if expected is not None and not np.all(reaction_strip == expected):
                raise RuntimeError(f"{path}: expected beam truth strip -1")
            if expected is None and np.any(reaction_strip < 0):
                raise RuntimeError(f"{path}: reaction-class file has beam events")
            rows.append(raw / beam_ref[np.newaxis, :])
            row_strips.append(reaction_strip)
        per_class[name] = (np.concatenate(rows).astype(np.float32),
                           np.concatenate(row_strips))
    n_each = min(X.shape[0] for X, _strip in per_class.values())
    if max_per_class is not None:
        n_each = min(n_each, max_per_class)
    if n_each == 0:
        raise RuntimeError("at least one configured simulation class has no events")
    traces = []
    labels = []
    strips = []
    for name in config.VLM_CLASSES:
        X, reaction_strip = per_class[name]
        take = _sample_indices(X.shape[0], n_each, rng)
        X, reaction_strip = X[take], reaction_strip[take]
        traces.append(X)
        labels.append(np.full(X.shape[0], list(config.VLM_CLASSES).index(name),
                              dtype=np.int64))
        strips.append(reaction_strip)
        print(f"  {name:>5}: {X.shape[0]} balanced examples from "
              f"{len(files[name])} file(s)")

    X = np.concatenate(traces)
    y = np.concatenate(labels)
    reaction_strip = np.concatenate(strips)
    holdout = np.isin(reaction_strip, config.VLM_FINETUNE_VAL_STRIPS)
    beam_rows = reaction_strip == -1
    beam_holdout = rng.random(int(beam_rows.sum())) < 0.2
    holdout[beam_rows] = beam_holdout
    if not holdout.any() or holdout.all():
        raise RuntimeError("empty train or validation split; adjust the simulated "
                           "strip coverage or VLM_FINETUNE_VAL_STRIPS")
    print(f"  split: {int((~holdout).sum())} train, {int(holdout.sum())} validation "
          f"(held-out reaction strips {config.VLM_FINETUNE_VAL_STRIPS})")
    return X, y, holdout, beam_ref.astype(np.float32), reaction_strip


def _load_base_model(model_id, load_in):
    """Load the E-series or unified Gemma auto class, optionally quantized."""
    import torch
    import transformers

    kw = {"device_map": "auto", "token": vlm._hf_token()}
    if load_in:
        from transformers import BitsAndBytesConfig
        if load_in != "4bit":
            raise ValueError("fine-tuning currently supports VLM_FINETUNE_LOAD_IN="
                             "'4bit' or None")
        kw["quantization_config"] = BitsAndBytesConfig(
            load_in_4bit=True, bnb_4bit_quant_type="nf4",
            bnb_4bit_compute_dtype=torch.bfloat16)
    else:
        kw["dtype"] = torch.bfloat16
    last = None
    for name in ("AutoModelForImageTextToText", "AutoModelForMultimodalLM"):
        model_class = getattr(transformers, name, None)
        if model_class is None:
            continue
        try:
            model = model_class.from_pretrained(model_id, **kw)
            print(f"  loaded via {name}")
            return model
        except (KeyError, TypeError, ValueError) as exc:
            last = exc
    raise RuntimeError(f"no auto class could load {model_id}: {last}")


def _load_trainable_model(model_id, load_in):
    """Load Gemma with a trainable LoRA adapter, including 4-bit QLoRA."""
    from peft import LoraConfig, get_peft_model, prepare_model_for_kbit_training

    model = _load_base_model(model_id, load_in)
    if load_in:
        model = prepare_model_for_kbit_training(model,
                                                use_gradient_checkpointing=True)
    else:
        model.gradient_checkpointing_enable()
        model.enable_input_require_grads()
    # Bare suffixes also match Gemma4ClippableLinear wrappers in the vision and
    # audio towers, which PEFT cannot wrap. Restrict LoRA to the quantized
    # language-model projections; the visual representation still trains
    # through how those tokens are interpreted by the LM.
    targets = (r"model\.language_model\.layers\.\d+\."
               r"(self_attn\.(q_proj|k_proj|v_proj|o_proj)|"
               r"mlp\.(gate_proj|up_proj|down_proj))")
    adapter = LoraConfig(r=config.VLM_FINETUNE_LORA_RANK,
                         lora_alpha=config.VLM_FINETUNE_LORA_ALPHA,
                         lora_dropout=config.VLM_FINETUNE_LORA_DROPOUT,
                         target_modules=targets,
                         bias="none")
    model = get_peft_model(model, adapter)
    model.print_trainable_parameters()
    return model


def _model_device(model):
    """Input device for an accelerate-dispatched model, or its first parameter."""
    if hasattr(model, "hf_device_map"):
        for device in model.hf_device_map.values():
            if device not in ("cpu", "disk"):
                return device
    return next(model.parameters()).device


def _model_inputs(builder, images, model, dtype):
    inputs = builder.build_inputs(images).to(_model_device(model))
    if "pixel_values" in inputs:
        inputs["pixel_values"] = inputs["pixel_values"].to(dtype)
    return inputs


def _forward_kwargs(model):
    """Keep only final-position logits when this transformers version supports it."""
    parameters = inspect.signature(model.forward).parameters
    for name in ("logits_to_keep", "num_logits_to_keep"):
        if name in parameters:
            return {name: 1}
    return {}


def translate_from_beam(X, beam_ref, shifts):
    """Translate each trace's deviation from beam by an integer strip count."""
    X = np.asarray(X, dtype=np.float32)
    beam_ref = np.asarray(beam_ref, dtype=np.float32)
    out = np.broadcast_to(beam_ref, X.shape).copy()
    delta = X - beam_ref[np.newaxis, :]
    for i, shift in enumerate(shifts):
        if shift > 0:
            out[i, shift:] += delta[i, :-shift]
        elif shift < 0:
            out[i, :shift] += delta[i, -shift:]
        else:
            out[i] = X[i]
    return out


def _batch_loss(model, builder, class_ids, X, y, beam_ref, rows, train,
                rng=None):
    import torch

    traces = X[rows]
    if train and config.VLM_FINETUNE_SHIFT_STRIPS > 0:
        shifts = rng.integers(-config.VLM_FINETUNE_SHIFT_STRIPS,
                             config.VLM_FINETUNE_SHIFT_STRIPS + 1,
                             size=rows.size)
        traces = translate_from_beam(traces, beam_ref, shifts)
    images = vlm.render_traces(traces, beam_ref)
    inputs = _model_inputs(builder, images, model, torch.bfloat16)
    context = torch.enable_grad() if train else torch.inference_mode()
    with context:
        logits = model(**inputs, use_cache=False,
                       **_forward_kwargs(model)).logits[:, -1, :].float()
        chosen = logits[:, class_ids]
        targets = torch.as_tensor(y[rows], device=chosen.device)
        loss = torch.nn.functional.cross_entropy(chosen, targets)
    probabilities = torch.softmax(chosen, dim=1).detach().cpu().numpy()
    return loss, probabilities.argmax(axis=1), probabilities


def evaluate(model, builder, class_ids, X, y, beam_ref, rows, batch):
    """Return loss, confusion matrix, and per-class recall on held-out traces."""
    n_cls = len(config.VLM_CLASSES)
    confusion = np.zeros((n_cls, n_cls), dtype=np.int64)
    probabilities = np.empty((rows.size, n_cls), dtype=np.float32)
    loss_sum = 0.0
    model.eval()
    cursor = 0
    for start in range(0, rows.size, batch):
        subset = rows[start:start + batch]
        loss, predicted, probs = _batch_loss(model, builder, class_ids, X, y,
                                             beam_ref, subset, train=False)
        loss_sum += float(loss) * subset.size
        probabilities[cursor:cursor + subset.size] = probs
        cursor += subset.size
        np.add.at(confusion, (y[subset], predicted), 1)
    recall = np.divide(np.diag(confusion), confusion.sum(axis=1),
                       out=np.full(n_cls, np.nan),
                       where=confusion.sum(axis=1) > 0)
    if rows.size == 0:
        raise RuntimeError("cannot evaluate an empty validation selection")
    return loss_sum / rows.size, confusion, recall, probabilities


def threshold_metrics(probabilities, truth):
    """Best p(an) cut subject to a cap on all simulated background leakage."""
    an_index = list(config.VLM_CLASSES).index("an")
    beam_index = list(config.VLM_CLASSES).index("beam")
    pan = probabilities[:, an_index]
    signal = truth == an_index
    background = ~signal
    beam = truth == beam_index
    candidates = np.unique(np.concatenate(([0.0, 1.0], pan)))
    best = None
    for threshold in candidates:
        accepted = pan >= threshold
        recall = float(accepted[signal].mean())
        background_fpr = float(accepted[background].mean())
        beam_fpr = float(accepted[beam].mean())
        if background_fpr <= config.VLM_FINETUNE_MAX_BACKGROUND_FPR:
            candidate = (recall, -background_fpr, float(threshold), beam_fpr)
            if best is None or candidate[:3] > best[:3]:
                best = candidate
    if best is None:
        return 1.0, 0.0, 0.0, 0.0
    return best[2], best[0], -best[1], best[3]


def _print_evaluation(loss, confusion, recall, probabilities, truth):
    threshold, an_recall, background_fpr, beam_fpr = threshold_metrics(
        probabilities, truth)
    print(f"  validation loss: {loss:.4f}")
    print("  rows = truth, columns = prediction")
    print("           " + " ".join(f"{name:>7}" for name in config.VLM_CLASSES))
    for i, name in enumerate(config.VLM_CLASSES):
        cells = " ".join(f"{n:7d}" for n in confusion[i])
        print(f"  {name:>7} {cells}   recall {recall[i]:.3f}")
    print(f"  selected p(an)>={threshold:.4f}: an recall {an_recall:.3f}, "
          f"background false-positive {background_fpr:.3f}, "
          f"beam false-positive {beam_fpr:.3f}")
    return threshold, an_recall, background_fpr


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--epochs", type=int, default=config.VLM_FINETUNE_EPOCHS)
    ap.add_argument("--max-per-class", type=int,
                    default=config.VLM_FINETUNE_MAX_PER_CLASS)
    ap.add_argument("--adapter", type=Path, default=None,
                    help="existing LoRA adapter to evaluate instead of training")
    ap.add_argument("--evaluate", action="store_true",
                    help="evaluate --adapter on the held-out simulator traces")
    args = ap.parse_args()
    if args.evaluate and args.adapter is None:
        ap.error("--evaluate requires --adapter")
    if args.adapter is not None and not args.adapter.is_dir():
        ap.error(f"adapter directory does not exist: {args.adapter}")

    X, y, holdout, beam_ref, _strips = load_examples(args.max_per_class)
    train_rows = np.flatnonzero(~holdout)
    validation_rows = np.flatnonzero(holdout)
    for i, name in enumerate(config.VLM_CLASSES):
        n_train = int((y[train_rows] == i).sum())
        n_validation = int((y[validation_rows] == i).sum())
        if n_train == 0 or n_validation == 0:
            raise RuntimeError(f"{name}: train={n_train}, validation={n_validation}; "
                               "each class needs examples in both splits")
    builder = vlm.PromptBuilder(model_id=config.VLM_FINETUNE_MODEL)
    class_ids = builder.class_token_ids

    if args.adapter is not None:
        from peft import PeftModel
        base = _load_base_model(config.VLM_FINETUNE_MODEL, load_in=None)
        model = PeftModel.from_pretrained(base, str(args.adapter))
        loss, confusion, recall, probabilities = evaluate(
            model, builder, class_ids, X, y, beam_ref, validation_rows,
            config.VLM_FINETUNE_BATCH)
        _print_evaluation(loss, confusion, recall, probabilities,
                          y[validation_rows])
        if args.evaluate:
            return
        raise RuntimeError("an adapter was supplied without --evaluate; refusing "
                           "to overwrite it")
    if args.evaluate:
        ap.error("--evaluate requires --adapter")

    if args.epochs < 1:
        ap.error("--epochs must be at least 1")
    if args.max_per_class is not None and args.max_per_class < 1:
        ap.error("--max-per-class must be at least 1")

    import torch
    model = _load_trainable_model(config.VLM_FINETUNE_MODEL,
                                  config.VLM_FINETUNE_LOAD_IN)
    trainable = [p for p in model.parameters() if p.requires_grad]
    if not trainable:
        raise RuntimeError("LoRA attached no trainable parameters; check the "
                           "Gemma module names against the installed model")
    optimizer = torch.optim.AdamW(trainable, lr=config.VLM_FINETUNE_LR)
    steps_per_epoch = int(np.ceil(train_rows.size / config.VLM_FINETUNE_BATCH))
    total_steps = steps_per_epoch * args.epochs
    warmup_steps = int(total_steps * config.VLM_FINETUNE_WARMUP_FRAC)
    scheduler = None
    if warmup_steps:
        scheduler = torch.optim.lr_scheduler.LinearLR(
            optimizer, start_factor=0.1, end_factor=1.0,
            total_iters=warmup_steps)
    rng = np.random.default_rng(config.VLM_FINETUNE_SEED)
    best_score = None
    best_epoch = -1
    best_threshold = config.VLM_AN_THRESHOLD
    for epoch in range(args.epochs):
        order = rng.permutation(train_rows)
        model.train()
        optimizer.zero_grad(set_to_none=True)
        running = 0.0
        for step, start in enumerate(range(0, order.size,
                                           config.VLM_FINETUNE_BATCH), 1):
            rows = order[start:start + config.VLM_FINETUNE_BATCH]
            loss, _pred, _probabilities = _batch_loss(
                model, builder, class_ids, X, y, beam_ref, rows, train=True,
                rng=rng)
            (loss / config.VLM_FINETUNE_GRAD_ACCUM).backward()
            running += float(loss.detach())
            if step % config.VLM_FINETUNE_GRAD_ACCUM == 0 or step == steps_per_epoch:
                optimizer.step()
                optimizer.zero_grad(set_to_none=True)
                if scheduler is not None and scheduler.last_epoch < warmup_steps:
                    scheduler.step()
        val_loss, confusion, recall, probabilities = evaluate(
            model, builder, class_ids, X, y, beam_ref, validation_rows,
            config.VLM_FINETUNE_BATCH)
        print(f"epoch {epoch + 1}/{args.epochs}: train loss "
              f"{running / steps_per_epoch:.4f}")
        threshold, an_recall, background_fpr = _print_evaluation(
            val_loss, confusion, recall, probabilities, y[validation_rows])
        checkpoint = config.VLM_FINETUNE_OUTPUT_DIR / f"epoch-{epoch + 1}"
        checkpoint.mkdir(parents=True, exist_ok=True)
        model.save_pretrained(checkpoint)
        score = (an_recall, -background_fpr, -val_loss)
        if epoch == 0 or score > best_score:
            best_score = score
            best_epoch = epoch + 1
            best_threshold = threshold

    output = config.VLM_FINETUNE_OUTPUT_DIR
    best = output / f"epoch-{best_epoch}"
    print(f"selected epoch {best_epoch}: p(an)>={best_threshold:.4f}, "
          f"score {best_score}")
    metadata = {
        "base_model": config.VLM_FINETUNE_MODEL,
        "classes": list(config.VLM_CLASSES),
        "class_tokens": list(config.VLM_CLASS_TOKENS),
        "validation_strips": list(config.VLM_FINETUNE_VAL_STRIPS),
        "train_examples": int(train_rows.size),
        "validation_examples": int(validation_rows.size),
        "translation_strips": config.VLM_FINETUNE_SHIFT_STRIPS,
        "selected_epoch": best_epoch,
        "selected_adapter": str(best),
        "selected_threshold": best_threshold,
        "maximum_background_false_positive":
        config.VLM_FINETUNE_MAX_BACKGROUND_FPR,
    }
    (output / "music_finetune_metadata.json").write_text(
        json.dumps(metadata, indent=2))
    print(f"saved selected LoRA adapter -> {best}")


if __name__ == "__main__":
    main()
