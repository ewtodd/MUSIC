"""Quantify simulator-to-experiment trace mismatch by reaction strip.

The experimental side is the compute-regions cache and therefore describes the
existing tag's population, not labelled truth. The comparison is still useful:
it identifies which trace components change between the simulator context and
the data on which TabFM must transfer.
"""

import argparse

import numpy as np

import config
import vlm_finetune


def trace_features(X, reaction_strip):
    """Reaction-relative rise, plateau, tail, and collapse for each trace."""
    result = np.empty((X.shape[0], 5), dtype=np.float64)
    for i, strip in enumerate(reaction_strip):
        before_lo = max(0, strip - 3)
        before = X[i, before_lo:strip]
        plateau = X[i, strip:min(X.shape[1], strip + 4)]
        tail = X[i, max(strip + 1, X.shape[1] - 3):]
        result[i, 0] = before.mean()
        result[i, 1] = plateau.mean()
        result[i, 2] = tail.mean()
        result[i, 3] = result[i, 1] - result[i, 0]
        result[i, 4] = result[i, 1] - result[i, 2]
    return result


def summarize(name, values):
    med = np.median(values, axis=0)
    return (f"{name:<10} n={values.shape[0]:5d} pre={med[0]:.3f} "
            f"plateau={med[1]:.3f} tail={med[2]:.3f} "
            f"rise={med[3]:+.3f} collapse={med[4]:+.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-per-class", type=int, default=2000)
    args = parser.parse_args()

    X_sim, y_sim, _holdout, _beam, strip_sim = vlm_finetune.load_examples(
        args.max_per_class)
    an_index = list(config.VLM_CLASSES).index("an")
    sim_an = y_sim == an_index

    path = config.CACHE_DIR / "tabfm_compute_regions_input.npz"
    if not path.is_file():
        raise FileNotFoundError(f"run vlm_tabfm.py --reservoir first: {path}")
    experimental = np.load(path)
    X_data = experimental["X"]
    strip_data = experimental["reac"].astype(np.int32)

    if X_sim.shape[1] - X_data.shape[1] == 1:
        X_sim = X_sim[:, :-1]
    if X_sim.shape[1] != X_data.shape[1]:
        raise RuntimeError(f"simulator has {X_sim.shape[1]} live strips, data has "
                           f"{X_data.shape[1]}")

    print("reaction-relative medians")
    print("  pre: three strips before reaction; plateau: reaction through +3")
    print("  tail: final three live strips; collapse: plateau - tail")
    rows = []
    for strip in sorted(set(int(value) for value in strip_data if value >= 0)):
        sim_rows = sim_an & (strip_sim == strip)
        data_rows = strip_data == strip
        if not sim_rows.any() or not data_rows.any():
            continue
        sim_features = trace_features(X_sim[sim_rows], strip_sim[sim_rows])
        data_features = trace_features(X_data[data_rows], strip_data[data_rows])
        sim_median = np.median(sim_features, axis=0)
        data_median = np.median(data_features, axis=0)
        trace_rmse = float(np.sqrt(np.mean(
            (np.median(X_sim[sim_rows], axis=0)
             - np.median(X_data[data_rows], axis=0)) ** 2)))
        print(f"strip {strip}")
        print("  " + summarize("sim an", sim_features))
        print("  " + summarize("data tag", data_features))
        print(f"  data-sim rise={data_median[3] - sim_median[3]:+.3f}, "
              f"collapse={data_median[4] - sim_median[4]:+.3f}, "
              f"median-trace RMSE={trace_rmse:.3f}")
        rows.append(np.concatenate(([strip, trace_rmse], sim_median, data_median)))

    output = config.CACHE_DIR / "vlm_domain_by_strip.npz"
    np.savez(output, columns=np.array([
        "strip", "trace_rmse", "sim_pre", "sim_plateau", "sim_tail",
        "sim_rise", "sim_collapse", "data_pre", "data_plateau",
        "data_tail", "data_rise", "data_collapse"]), rows=np.asarray(rows))
    print(f"summary -> {output}")


if __name__ == "__main__":
    main()
