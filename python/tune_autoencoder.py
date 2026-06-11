"""Autoencoder hyperparameter study for the active feature configuration.

For the feature knobs currently set in config.py this runs:

1. one-at-a-time sweeps (latent dim, hidden widths, learning rate, batch
size, and -- with AE_OUTLIER_EXPOSURE -- the OE margin ratio and
hinge weight) around the AE_* config values, multi-seed, scored by
AUC of reacted (aa+an) vs beam-validation reconstruction error;
2. a randomized search around the OAT-recommended values to explore 
interactions;
3. a KernelSHAP pass on the reconstruction-error score of the winner,
plotting mean |SHAP| vs strip index to expose which strips drive the
anomaly score.

Training mirrors train_autoencoder: with AE_OUTLIER_EXPOSURE the fits see
the pooled OE negatives, and those populations contribute only their
held-out slice to the AUC evaluation pool.

Every (configuration, seed) result is cached under ml_cache/ae_tune/ keyed
by the feature/normalization/OE knobs, so interrupted or repeated runs
only compute what is missing and switching knobs never reuses stale
results.

Run inside the dataset dev shell, from python/:

    python -m music_ml.tune_autoencoder
"""

import copy
import json
import math

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

import config, data, mlplots, util
from models import BeamAutoencoder

TUNE_CACHE = config.CACHE_DIR / "ae_tune"

# Tune only the feature configuration currently set in config.py; the knob
# combos are explored by hand, the harness tunes hyperparameters for the
# active one.
COMBOS = [{
    "short": config.INCLUDE_SHORT_STRIPS,
    "guard": config.INCLUDE_GUARD_STRIPS,
    "deriv": config.INCLUDE_DERIVATIVE,
}]

SWEEPS = [
    ("latent_dim", config.TUNE_LATENT_VALUES, "Latent dimension", False),
    ("hidden_dims", config.TUNE_HIDDEN_VALUES, "Hidden widths", False),
    ("learning_rate", config.TUNE_LR_VALUES, "Learning rate", True),
    ("batch_size", config.TUNE_BATCH_VALUES, "Batch size", True),
]
if config.AE_OUTLIER_EXPOSURE:
    SWEEPS.append(("oe_margin_ratio", config.TUNE_OE_RATIO_VALUES,
                   "OE margin ratio", True))
    SWEEPS.append(
        ("oe_weight", config.TUNE_OE_WEIGHT_VALUES, "OE hinge weight", True))


def combo_tag(combo):
    # The global normalization/OE knobs are folded into the tag so cached
    # results never silently cross settings (the OE holdout changes the
    # evaluation pool, hence it is part of the key; the swept OE ratio /
    # weight live in hp_key instead).
    tag = (f"sh{int(combo['short'])}"
           f"_gd{int(combo['guard'])}"
           f"_dv{int(combo['deriv'])}"
           f"_sb{int(config.SUBTRACT_BEAM_LEVEL)}")
    if config.AE_OUTLIER_EXPOSURE:
        tag += (f"_oe{'+'.join(config.AE_OE_POPULATIONS)}"
                f"_ho{config.AE_OE_HOLDOUT_FRACTION:g}")
    return tag


def combo_label(combo):
    parts = []
    if combo["short"]:
        parts.append("short")
    if combo["guard"]:
        parts.append("guard")
    if combo["deriv"]:
        parts.append("deriv")
    return "+".join(parts) if parts else "long-only"


def apply_combo(combo):
    config.INCLUDE_SHORT_STRIPS = combo["short"]
    config.INCLUDE_GUARD_STRIPS = combo["guard"]
    config.INCLUDE_DERIVATIVE = combo["deriv"]


def hp_key(hp):
    hidden = "x".join(str(w) for w in hp["hidden_dims"])
    key = (f"l{hp['latent_dim']}_h{hidden}"
           f"_lr{hp['learning_rate']:.2e}_b{hp['batch_size']}")
    if "oe_margin_ratio" in hp:
        key += f"_R{hp['oe_margin_ratio']:g}_w{hp['oe_weight']:g}"
    return key


def base_hp():
    hp = {
        "latent_dim": config.AE_LATENT_DIM,
        "hidden_dims": tuple(config.AE_HIDDEN_DIMS),
        "learning_rate": config.AE_LEARNING_RATE,
        "batch_size": config.AE_BATCH_SIZE,
    }
    if config.AE_OUTLIER_EXPOSURE:
        hp["oe_margin_ratio"] = config.AE_OE_MARGIN_RATIO
        hp["oe_weight"] = config.AE_OE_WEIGHT
    return hp


def load_combo_data():
    """(X_beam, X_oe, X_reacted) AE feature arrays for the active combo.

    X_oe pools the outlier-exposure training negatives (None with
    AE_OUTLIER_EXPOSURE off). X_reacted concatenates every (a,a') and
    (a,n) population for the AUC, each capped at TUNE_EVAL_CAP_PER_POP
    traces (deterministic subsample) so evaluation stays cheap and
    population-balanced; OE populations contribute only their held-out
    slice, mirroring train_autoencoder.
    """
    specs = data.list_sim_specs(eres_only=True)
    gain = data.beam_gains(specs)
    beam_spec = data.find_beam_spec(specs)
    X_beam = data.ae_features(data.load_population(beam_spec, gain))
    rng = np.random.default_rng(config.SEED)
    oe_parts = []
    reacted = []
    for spec in specs:
        if spec.strip < 0 or spec.base not in ("aa", "an"):
            continue
        X = data.ae_features(data.load_population(spec, gain))
        if (config.AE_OUTLIER_EXPOSURE
                and spec.base in config.AE_OE_POPULATIONS):
            tr_idx, held_idx = util.split_indices(
                X.shape[0], config.AE_OE_HOLDOUT_FRACTION, rng)
            oe_parts.append(X[tr_idx])
            X = X[held_idx]
        if X.shape[0] > config.TUNE_EVAL_CAP_PER_POP:
            sel = rng.choice(X.shape[0],
                             config.TUNE_EVAL_CAP_PER_POP,
                             replace=False)
            X = X[sel]
        reacted.append(X)
    if not reacted:
        raise RuntimeError("no reacted (aa/an) sim populations found")
    if config.AE_OUTLIER_EXPOSURE and not oe_parts:
        raise RuntimeError(f"AE_OUTLIER_EXPOSURE: no sim populations match "
                           f"{config.AE_OE_POPULATIONS}")
    X_oe = np.concatenate(oe_parts) if oe_parts else None
    return X_beam, X_oe, np.concatenate(reacted)


def train_ae(X_train, X_val, hp, seed, device, X_oe=None):
    """One reduced-budget fit; returns the best-val-loss model.

    With X_oe (pooled outlier-exposure negatives) the training loss and
    the early-stopping objective both gain the hinge term, exactly like
    train_autoencoder; a tenth of the pool is carved off as the
    early-stopping slice.
    """
    util.set_seed(seed)
    model = BeamAutoencoder(n_inputs=X_train.shape[1],
                            hidden_dims=hp["hidden_dims"],
                            latent_dim=hp["latent_dim"]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=hp["learning_rate"])
    loader = DataLoader(TensorDataset(torch.from_numpy(X_train)),
                        batch_size=hp["batch_size"],
                        shuffle=True)
    x_val_dev = torch.from_numpy(X_val).to(device)

    oe_iter = None
    x_oe_stop_dev = None
    if X_oe is not None:
        oe_rng = np.random.default_rng(seed)
        oe_tr_idx, oe_stop_idx = util.split_indices(X_oe.shape[0], 0.1, oe_rng)
        oe_iter = util.cycle_batches(
            DataLoader(TensorDataset(torch.from_numpy(X_oe[oe_tr_idx])),
                       batch_size=hp["batch_size"],
                       shuffle=True))
        x_oe_stop_dev = torch.from_numpy(X_oe[oe_stop_idx]).to(device)

    best_val = float("inf")
    best_state = None
    stale = 0
    for _epoch in range(1, config.TUNE_MAX_EPOCHS + 1):
        model.train()
        for (xb, ) in loader:
            xb = xb.to(device)
            optimizer.zero_grad()
            err_beam_batch = util.per_trace_mse(model, xb)
            loss = err_beam_batch.mean()
            if oe_iter is not None:
                (nb, ) = next(oe_iter)
                nb = nb.to(device)
                loss = loss + hp["oe_weight"] * util.oe_hinge(
                    util.per_trace_mse(model, nb),
                    err_beam_batch.mean().detach(), hp["oe_margin_ratio"])
            loss.backward()
            optimizer.step()
        model.eval()
        with torch.no_grad():
            val_beam_mse = util.per_trace_mse(model, x_val_dev).mean()
            val_loss = val_beam_mse.item()
            if oe_iter is not None:
                val_loss += hp["oe_weight"] * util.oe_hinge(
                    util.per_trace_mse(model, x_oe_stop_dev), val_beam_mse,
                    hp["oe_margin_ratio"]).item()
        if val_loss < best_val:
            best_val = val_loss
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if stale >= config.TUNE_PATIENCE:
            break
    model.load_state_dict(best_state)
    return model


def reconstruction_errors(model, X, device, chunk=65536):
    model.eval()
    out = np.empty(X.shape[0], dtype=np.float64)
    for lo in range(0, X.shape[0], chunk):
        xb = torch.from_numpy(X[lo:lo + chunk]).to(device)
        out[lo:lo + chunk] = model.reconstruction_error(xb).cpu().numpy()
    return out


def fit_and_auc(X_beam, X_oe, X_reacted, hp, seed, device):
    """Train at one seed, return AUC(reacted vs beam-val recon error)."""
    rng = np.random.default_rng(seed)
    train_idx, val_idx = util.split_indices(X_beam.shape[0],
                                            config.AE_VAL_FRACTION, rng)
    if train_idx.size > config.TUNE_TRAIN_CAP:
        train_idx = train_idx[:config.TUNE_TRAIN_CAP]
    if X_oe is not None and X_oe.shape[0] > config.TUNE_TRAIN_CAP:
        X_oe = X_oe[rng.choice(X_oe.shape[0],
                               config.TUNE_TRAIN_CAP,
                               replace=False)]
    model = train_ae(X_beam[train_idx],
                     X_beam[val_idx],
                     hp,
                     seed,
                     device,
                     X_oe=X_oe)
    err_beam = reconstruction_errors(model, X_beam[val_idx], device)
    err_reacted = reconstruction_errors(model, X_reacted, device)
    return util.auc(err_reacted, err_beam)


def evaluate_hp(tag, X_beam, X_oe, X_reacted, hp, device):
    """Multi-seed AUC for one configuration, npz-cached per (combo, hp)."""
    seeds = list(config.TUNE_SEEDS)
    cache = TUNE_CACHE / f"{tag}_{hp_key(hp)}.npz"
    if cache.exists():
        d = np.load(cache)
        if list(d["seeds"]) == seeds:
            return float(d["aucs"].mean()), float(d["aucs"].std()), d["aucs"]
    aucs = np.array([
        fit_and_auc(X_beam, X_oe, X_reacted, hp, seed, device)
        for seed in seeds
    ])
    TUNE_CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(cache, aucs=aucs, seeds=np.array(seeds))
    return float(aucs.mean()), float(aucs.std()), aucs


def run_oat(tag, X_beam, X_oe, X_reacted, device):
    """OAT sweeps; returns the recommended hp dict (argmax mean AUC each)."""
    recommended = base_hp()
    for param, values, x_title, log_x in SWEEPS:
        means = []
        stds = []
        for value in values:
            hp = dict(base_hp())
            hp[param] = value
            m, s, _ = evaluate_hp(tag, X_beam, X_oe, X_reacted, hp, device)
            means.append(m)
            stds.append(s)
            print(f"  [{tag}] {param}={value}  AUC = {m:.4f} +/- {s:.4f}")
        best = int(np.argmax(means))
        recommended[param] = values[best]
        print(f"  [{tag}] >> recommended {param} = {values[best]} "
              f"(AUC = {means[best]:.4f})")
        if param == "hidden_dims":
            labels = ["x".join(str(w) for w in v) for v in values]
            mlplots.labeled_points(labels, means, stds, "ROC AUC",
                                   f"ae_tune_{tag}_oat_{param}")
        else:
            mlplots.sweep_plot(list(values),
                               means,
                               stds,
                               x_title,
                               f"ae_tune_{tag}_oat_{param}",
                               log_x=log_x)
    return recommended


def sample_hp(recommended, rng):
    """One randomized-search configuration around the OAT recommendation."""
    rec_latent = recommended["latent_dim"]
    latent = int(rng.integers(max(1, rec_latent - 2), rec_latent + 3))
    w0 = recommended["hidden_dims"][0]
    w = int(rng.integers(max(6, w0 // 2), 2 * w0 + 1))
    hidden = [w]
    for _ in recommended["hidden_dims"][1:]:
        w = max(4, w // 2)
        hidden.append(w)
    rec_lr = recommended["learning_rate"]
    lr = float(rec_lr * math.exp(rng.uniform(math.log(0.2), math.log(5.0))))
    batch = int(rng.choice(np.array(config.TUNE_BATCH_VALUES)))
    hp = {
        "latent_dim": latent,
        "hidden_dims": tuple(hidden),
        "learning_rate": lr,
        "batch_size": batch,
    }
    if "oe_margin_ratio" in recommended:
        hp["oe_margin_ratio"] = float(
            recommended["oe_margin_ratio"] *
            math.exp(rng.uniform(math.log(0.3), math.log(3.0))))
        hp["oe_weight"] = float(
            recommended["oe_weight"] *
            math.exp(rng.uniform(math.log(0.3), math.log(3.0))))
    return hp


def run_randomized_search(tag, X_beam, X_oe, X_reacted, recommended, device):
    """Random configs around the recommendation; returns (best_hp, m, s)."""
    rng = np.random.default_rng(config.SEED)
    candidates = [recommended]
    for _ in range(config.TUNE_RS_ITER):
        candidates.append(sample_hp(recommended, rng))
    best = None
    for hp in candidates:
        m, s, _ = evaluate_hp(tag, X_beam, X_oe, X_reacted, hp, device)
        print(f"  [{tag}] rs {hp_key(hp)}  AUC = {m:.4f} +/- {s:.4f}")
        if best is None or m > best[1]:
            best = (hp, m, s)
    return best


def write_combo_table(rows):
    """LaTeX ranking table"""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Best beam-autoencoder AUC (reacted vs beam) "
                 r"per input-feature combination, after OAT sweeps and "
                 r"randomized search. Errors are the seed spread.}")
    lines.append(r"  \label{tab:ae_tune}")
    lines.append(r"  \begin{tabular}{lccc}")
    lines.append(r"    \hline \hline")
    lines.append(r"    Features & Inputs & Best configuration & "
                 r"AUC \\")
    lines.append(r"    \hline")
    for row in rows:
        hp = row["best_hp"]
        hidden = r"$\times$".join(str(w) for w in hp["hidden_dims"])
        cfg = (f"latent {hp['latent_dim']}, {hidden}, "
               f"lr {hp['learning_rate']:.1e}, batch {hp['batch_size']}")
        if "oe_margin_ratio" in hp:
            cfg += (f", OE R {hp['oe_margin_ratio']:g} "
                    f"w {hp['oe_weight']:g}")
        lines.append(f"    {row['label']} & {row['n_inputs']} & {cfg} & "
                     f"{row['auc']:.4f} $\\pm$ {row['err']:.4f} \\\\")
    lines.append(r"    \hline \hline")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")
    table = "\n".join(lines)
    out = TUNE_CACHE / "combo_table.txt"
    with open(out, "w") as f:
        f.write(table + "\n")
    print(f"LaTeX combo table saved to {out}")
    print(table)


def shap_for_combo(tag, combo, hp, device):
    """Seed-averaged mean |SHAP| of the recon-error score, npy-cached.

    Background = beam sample, explain = reacted sample, KernelSHAP on the
    scalar reconstruction error so importances are per input feature.
    Returns (mean_abs, std_abs, avg_reacted_trace).
    """
    import shap

    apply_combo(combo)
    X_beam, X_oe, X_reacted = load_combo_data()
    rng = np.random.default_rng(config.SEED)
    bg = X_beam[rng.choice(X_beam.shape[0], config.TUNE_SHAP_BG,
                           replace=False)]
    explain = X_reacted[rng.choice(X_reacted.shape[0],
                                   config.TUNE_SHAP_EXPLAIN,
                                   replace=False)]

    per_seed = []
    for seed in config.TUNE_SEEDS:
        cache = TUNE_CACHE / f"{tag}_shap_{hp_key(hp)}_s{seed}.npy"
        if cache.exists():
            per_seed.append(np.load(cache))
            continue
        print(f"  [{tag}] SHAP seed {seed}: training + KernelSHAP...")
        srng = np.random.default_rng(seed)
        train_idx, val_idx = util.split_indices(X_beam.shape[0],
                                                config.AE_VAL_FRACTION, srng)
        if train_idx.size > config.TUNE_TRAIN_CAP:
            train_idx = train_idx[:config.TUNE_TRAIN_CAP]
        X_oe_fit = X_oe
        if X_oe_fit is not None and X_oe_fit.shape[0] > config.TUNE_TRAIN_CAP:
            X_oe_fit = X_oe_fit[srng.choice(X_oe_fit.shape[0],
                                            config.TUNE_TRAIN_CAP,
                                            replace=False)]
        model = train_ae(X_beam[train_idx],
                         X_beam[val_idx],
                         hp,
                         seed,
                         device,
                         X_oe=X_oe_fit)

        def score_fn(X):
            return reconstruction_errors(model, X.astype(np.float32), device)

        explainer = shap.KernelExplainer(score_fn, bg)
        values = explainer.shap_values(explain)
        mean_abs = np.mean(np.abs(values), axis=0)
        np.save(cache, mean_abs)
        per_seed.append(mean_abs)

    stacked = np.stack(per_seed, axis=0)
    n = config.n_strip_inputs()
    avg_reacted = explain[:, :n].mean(axis=0)
    return stacked.mean(axis=0), stacked.std(axis=0), avg_reacted


def shap_pass(rows, device):
    """SHAP-vs-strip plots for the top TUNE_TOP_COMBOS combos.

    One canvas per column block (long-side trace, then "_short" for the
    short-end block when present); the feature vector is laid out as the
    trace blocks followed by the block-local derivative blocks.
    """
    for row in rows[:config.TUNE_TOP_COMBOS]:
        combo = row["combo"]
        tag = row["tag"]
        print(f"SHAP pass: {row['label']} ({tag})")
        mean_abs, std_abs, avg_reacted = shap_for_combo(
            tag, combo, row["best_hp"], device)
        widths = config.block_widths()
        first = 0 if combo["guard"] else 1
        t_lo = 0
        d_lo = config.n_strip_inputs()
        for b in range(len(widths)):
            w = widths[b]
            suffix = "" if b == 0 else "_short"
            if b == 0:
                xs = np.arange(first, first + w, dtype=np.float64)
            else:
                xs = np.arange(1, 17, dtype=np.float64)  # short ends 1-16
            if combo["deriv"]:
                # Derivative feature i sits on the boundary between strips
                # i and i+1 -> plot at the midpoints.
                mlplots.shap_vs_strip(xs,
                                      mean_abs[t_lo:t_lo + w],
                                      std_abs[t_lo:t_lo + w],
                                      avg_reacted[t_lo:t_lo + w],
                                      f"ae_tune_{tag}_shap{suffix}",
                                      deriv_x=xs[:-1] + 0.5,
                                      deriv_shap=mean_abs[d_lo:d_lo + w - 1],
                                      deriv_err=std_abs[d_lo:d_lo + w - 1])
            else:
                mlplots.shap_vs_strip(xs, mean_abs[t_lo:t_lo + w],
                                      std_abs[t_lo:t_lo + w],
                                      avg_reacted[t_lo:t_lo + w],
                                      f"ae_tune_{tag}_shap{suffix}")
            t_lo += w
            d_lo += w - 1


def main():
    util.set_seed()
    device = util.get_device()
    print(f"music-ml: dataset {config.DATASET}, device {device.type}")
    print(f"tuning {len(COMBOS)} feature combos, seeds {config.TUNE_SEEDS}, "
          f"outlier exposure "
          f"{'on' if config.AE_OUTLIER_EXPOSURE else 'off'}")
    TUNE_CACHE.mkdir(parents=True, exist_ok=True)

    rows = []
    for combo in COMBOS:
        apply_combo(combo)
        tag = combo_tag(combo)
        label = combo_label(combo)
        print(f"=== combo {label} ({tag}): "
              f"{config.ae_input_dim()} inputs ===")
        X_beam, X_oe, X_reacted = load_combo_data()
        n_oe = 0 if X_oe is None else X_oe.shape[0]
        print(f"  beam {X_beam.shape[0]}, OE negatives {n_oe}, "
              f"reacted (eval) {X_reacted.shape[0]}")

        recommended = run_oat(tag, X_beam, X_oe, X_reacted, device)
        best_hp, best_auc, best_err = run_randomized_search(
            tag, X_beam, X_oe, X_reacted, recommended, device)
        print(f"  [{tag}] BEST {hp_key(best_hp)}  "
              f"AUC = {best_auc:.4f} +/- {best_err:.4f}")
        rows.append({
            "combo": combo,
            "tag": tag,
            "label": label,
            "n_inputs": config.ae_input_dim(),
            "best_hp": best_hp,
            "auc": best_auc,
            "err": best_err,
        })

    rows.sort(key=lambda r: r["auc"], reverse=True)
    print("=== combo ranking ===")
    for i, row in enumerate(rows):
        print(f"  {i + 1}. {row['label']:>22}  "
              f"AUC = {row['auc']:.4f} +/- {row['err']:.4f}  "
              f"({hp_key(row['best_hp'])})")

    mlplots.labeled_points([r["label"] for r in rows],
                           [r["auc"] for r in rows], [r["err"] for r in rows],
                           "Best ROC AUC", "ae_tune_combo_ranking")
    write_combo_table(rows)
    summary = [{
        "tag": r["tag"],
        "label": r["label"],
        "combo": r["combo"],
        "n_inputs": r["n_inputs"],
        "best_hp": {
            k: (list(v) if k == "hidden_dims" else v)
            for k, v in r["best_hp"].items()
        },
        "auc": r["auc"],
        "err": r["err"],
    } for r in rows]
    with open(TUNE_CACHE / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary saved to {TUNE_CACHE / 'summary.json'}")

    shap_pass(rows, device)


if __name__ == "__main__":
    main()
