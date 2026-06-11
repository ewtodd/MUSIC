"""Stage 1: train the beam autoencoder on simulated beam traces.

Trained on the max-normalized sim beam (_eres) population; the per-trace
mean squared reconstruction error is the anomaly score. Beam-like events
(including pileup, which max-norm collapses onto the flat beam shape)
reconstruct well; reacted events do not. The detection threshold is the
AE_THRESHOLD_PERCENTILE of the beam validation errors.

With AE_OUTLIER_EXPOSURE the sim reacted populations in AE_OE_POPULATIONS
are pooled as exposure negatives: a scale-free hinge term pushes their
reconstruction error above AE_OE_MARGIN_RATIO times the beam level, so the
model is trained not just to know beam but to NOT know reactions. Each
negative population is three-way split -- OE training, an early-stopping
slice, and a held-out test slice that alone feeds the reported AUC /
efficiencies (populations outside AE_OE_POPULATIONS are evaluated in
full, as before).

Run inside the dataset dev shell, from python/:

    python -m music_ml.train_autoencoder
"""

import copy
import json

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

import config, data, mlplots, util
from models import BeamAutoencoder
from plot_traces import draw_population


def _reconstruction_errors(model, X, device, chunk=65536):
    model.eval()
    out = np.empty(X.shape[0], dtype=np.float64)
    for lo in range(0, X.shape[0], chunk):
        xb = torch.from_numpy(X[lo:lo + chunk]).to(device)
        out[lo:lo + chunk] = model.reconstruction_error(xb).cpu().numpy()
    return out


def main():
    util.set_seed()
    device = util.get_device()
    print(f"music-ml: dataset {config.DATASET}, device {device.type}")

    specs = data.list_sim_specs(eres_only=True)
    gain = data.beam_gains(specs)
    beam_spec = data.find_beam_spec(specs)
    beam_trace = data.load_population(beam_spec, gain)
    X_beam = data.ae_features(beam_trace)
    print(f"beam population: {X_beam.shape[0]} traces from {beam_spec.tag} "
          f"({X_beam.shape[1]} inputs)")
    rng_plot = np.random.default_rng(config.SEED)
    draw_population(beam_trace, "sample_traces_beam", "Beam (sim)", rng_plot)

    rng = np.random.default_rng(config.SEED)
    train_idx, val_idx = util.split_indices(X_beam.shape[0],
                                            config.AE_VAL_FRACTION, rng)
    X_train = X_beam[train_idx]
    X_val = X_beam[val_idx]
    loader = DataLoader(TensorDataset(torch.from_numpy(X_train)),
                        batch_size=config.AE_BATCH_SIZE,
                        shuffle=True)
    x_val_dev = torch.from_numpy(X_val).to(device)

    # ---- outlier exposure: pooled sim negatives, three-way split ----
    oe_on = config.AE_OUTLIER_EXPOSURE
    oe_test_idx = {}  # spec.tag -> indices reserved for the final AUC
    oe_iter = None
    x_oe_stop_dev = None
    n_oe_train = 0
    if oe_on:
        oe_train_parts = []
        oe_stop_parts = []
        for spec in specs:
            if spec.strip < 0 or spec.base not in config.AE_OE_POPULATIONS:
                continue
            feats = data.ae_features(data.load_population(spec, gain))
            tr_idx, held_idx = util.split_indices(
                feats.shape[0], config.AE_OE_HOLDOUT_FRACTION, rng)
            half = held_idx.size // 2
            oe_train_parts.append(feats[tr_idx])
            oe_stop_parts.append(feats[held_idx[:half]])
            oe_test_idx[spec.tag] = held_idx[half:]
        if not oe_train_parts:
            raise RuntimeError(
                f"AE_OUTLIER_EXPOSURE: no sim populations match "
                f"{config.AE_OE_POPULATIONS}")
        X_oe_train = np.concatenate(oe_train_parts)
        X_oe_stop = np.concatenate(oe_stop_parts)
        n_oe_train = int(X_oe_train.shape[0])
        oe_iter = util.cycle_batches(
            DataLoader(TensorDataset(torch.from_numpy(X_oe_train)),
                       batch_size=config.AE_BATCH_SIZE,
                       shuffle=True))
        x_oe_stop_dev = torch.from_numpy(X_oe_stop).to(device)
        print(f"outlier exposure: {n_oe_train} negatives "
              f"({'/'.join(config.AE_OE_POPULATIONS)}) for training, "
              f"{X_oe_stop.shape[0]} for early stopping, "
              f"ratio {config.AE_OE_MARGIN_RATIO} "
              f"weight {config.AE_OE_WEIGHT}")

    model = BeamAutoencoder().to(device)
    optimizer = torch.optim.Adam(model.parameters(),
                                 lr=config.AE_LEARNING_RATE)

    best_val = float("inf")
    best_state = None
    stale = 0
    train_hist = []
    val_hist = []
    for epoch in range(1, config.AE_MAX_EPOCHS + 1):
        model.train()
        running = 0.0
        for (xb, ) in loader:
            xb = xb.to(device)
            optimizer.zero_grad()
            err_beam_batch = util.per_trace_mse(model, xb)
            loss = err_beam_batch.mean()
            if oe_iter is not None:
                (nb, ) = next(oe_iter)
                nb = nb.to(device)
                loss = loss + config.AE_OE_WEIGHT * util.oe_hinge(
                    util.per_trace_mse(model, nb),
                    err_beam_batch.mean().detach(), config.AE_OE_MARGIN_RATIO)
            loss.backward()
            optimizer.step()
            running += loss.item() * xb.size(0)
        train_loss = running / X_train.shape[0]
        model.eval()
        with torch.no_grad():
            # Early stopping tracks the full objective: beam val MSE plus
            # (with OE) the hinge on the dedicated early-stopping slice --
            # never on the held-out test slice that feeds the final AUC.
            val_beam_mse = util.per_trace_mse(model, x_val_dev).mean()
            val_loss = val_beam_mse.item()
            if oe_iter is not None:
                val_loss += config.AE_OE_WEIGHT * util.oe_hinge(
                    util.per_trace_mse(model, x_oe_stop_dev), val_beam_mse,
                    config.AE_OE_MARGIN_RATIO).item()
        train_hist.append(train_loss)
        val_hist.append(val_loss)
        if val_loss < best_val:
            best_val = val_loss
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if epoch % 10 == 0 or stale >= config.AE_PATIENCE:
            print(f"epoch {epoch:4d}: train {train_loss:.3e} "
                  f"val {val_loss:.3e} (best {best_val:.3e})")
        if stale >= config.AE_PATIENCE:
            print(f"early stop at epoch {epoch}")
            break
    model.load_state_dict(best_state)

    # ---- evaluation: beam (val) vs every reacted population ----
    # OE populations score only their held-out test slice (the trained-on
    # and early-stopping slices would flatter the AUC); everything else is
    # untouched by training and evaluates in full.
    err_beam = _reconstruction_errors(model, X_val, device)
    err_aa = []
    err_an = []
    for spec in specs:
        if spec.strip < 0:
            continue
        trace = data.load_population(spec, gain)
        if spec.base in ("aa", "an") and spec.strip == config.TRACE_STRIP:
            pretty = ("(#alpha,#alpha')"
                      if spec.base == "aa" else "(#alpha,n)")
            draw_population(trace, f"sample_traces_{spec.base}_s{spec.strip}",
                            f"{pretty} s{spec.strip} (sim)", rng_plot)
        feats = data.ae_features(trace)
        if spec.tag in oe_test_idx:
            feats = feats[oe_test_idx[spec.tag]]
        errs = _reconstruction_errors(model, feats, device)
        if spec.base == "aa":
            err_aa.append(errs)
        elif spec.base == "an":
            err_an.append(errs)
    err_aa = np.concatenate(err_aa) if err_aa else np.empty(0)
    err_an = np.concatenate(err_an) if err_an else np.empty(0)
    err_reacted = np.concatenate([err_aa, err_an])

    threshold = float(np.percentile(err_beam, config.AE_THRESHOLD_PERCENTILE))
    auc_value = util.auc(err_reacted, err_beam)
    auc_aa = util.auc(err_aa, err_beam) if err_aa.size else float("nan")
    auc_an = util.auc(err_an, err_beam) if err_an.size else float("nan")
    eff_aa = float(
        (err_aa > threshold).mean()) if err_aa.size else float("nan")
    eff_an = float(
        (err_an > threshold).mean()) if err_an.size else float("nan")
    print(f"AUC (reacted vs beam): {auc_value:.5f}  "
          f"(a,a') {auc_aa:.5f}  (a,n) {auc_an:.5f}")
    print(f"threshold (P{config.AE_THRESHOLD_PERCENTILE} beam val): "
          f"{threshold:.3e}")
    print(f"anomaly efficiency at threshold: (a,a') {eff_aa:.4f}  "
          f"(a,n) {eff_an:.4f}")

    # ---- plots ----
    mlplots.loss_curves(train_hist, val_hist, "ae_loss")
    mlplots.score_hists(
        [("Beam (val)", np.log10(np.maximum(err_beam, 1e-12))),
         ("(#alpha,#alpha')", np.log10(np.maximum(err_aa, 1e-12))),
         ("(#alpha,n)", np.log10(np.maximum(err_an, 1e-12)))],
        "ae_recon_error", "log_{10} reconstruction MSE")
    fpr, tpr = util.roc_points(err_reacted, err_beam)
    mlplots.roc_plot(fpr, tpr, auc_value, "ae_roc", "Anomaly efficiency")

    # ---- persist model + metadata ----
    config.MODELS_DIR.mkdir(parents=True, exist_ok=True)
    model_path = config.MODELS_DIR / "beam_autoencoder.pt"
    torch.save(model.state_dict(), model_path)
    meta = {
        "dataset": config.DATASET,
        "model": "BeamAutoencoder",
        "n_strips": config.N_STRIPS,
        "include_short_strips": config.INCLUDE_SHORT_STRIPS,
        "include_guard_strips": config.INCLUDE_GUARD_STRIPS,
        "include_derivative": config.INCLUDE_DERIVATIVE,
        "subtract_beam_level": config.SUBTRACT_BEAM_LEVEL,
        "n_strip_inputs": config.n_strip_inputs(),
        "n_inputs": config.ae_input_dim(),
        "hidden_dims": list(config.AE_HIDDEN_DIMS),
        "latent_dim": config.AE_LATENT_DIM,
        "seed": config.SEED,
        "sim_beam_file": beam_spec.path,
        "sim_beam_gains": gain.tolist(),
        "n_train": int(X_train.shape[0]),
        "n_val": int(X_val.shape[0]),
        "outlier_exposure": oe_on,
        "oe_populations": list(config.AE_OE_POPULATIONS) if oe_on else [],
        "oe_margin_ratio": config.AE_OE_MARGIN_RATIO if oe_on else None,
        "oe_weight": config.AE_OE_WEIGHT if oe_on else None,
        "oe_holdout_fraction":
        config.AE_OE_HOLDOUT_FRACTION if oe_on else None,
        "n_oe_train": n_oe_train,
        # With OE on this is the combined objective (beam MSE + weighted
        # hinge on the early-stopping slice), not a pure MSE.
        "best_val_objective": best_val,
        "threshold_percentile": config.AE_THRESHOLD_PERCENTILE,
        "threshold_mse": threshold,
        "auc_reacted_vs_beam": auc_value,
        "auc_aa_vs_beam": auc_aa,
        "auc_an_vs_beam": auc_an,
        "efficiency_at_threshold": {
            "aa": eff_aa,
            "an": eff_an
        },
        "torch_version": torch.__version__,
    }
    meta_path = config.MODELS_DIR / "beam_autoencoder.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"saved {model_path} and {meta_path}")


if __name__ == "__main__":
    main()
