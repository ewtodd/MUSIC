"""Sample-trace visualization.

Draws TRACE_N_SAMPLES randomly sampled traces per population EXACTLY as
the models see them -- the configured feature view (long/short column
blocks, guard-strip drop) applied by the shared data
pipeline -- with companion "<name>_deriv" (and, with the short-side block
enabled, "<name>_short") canvases.

Populations: sim beam, sim (a,a') and (a,n) at TRACE_STRIP, and (when
TRACE_DRAW_DATA is on and events files exist) a random sample of
experimental events.

Run inside the dataset dev shell, from python/:

    python -m music_ml.plot_traces
"""

import numpy as np

import config, data, mlplots


def draw_population(X, name, label, rng):
    """Sample TRACE_N_SAMPLES traces from X (trace columns only, NOT the
    ae_features vector) and draw them (+ derivative canvas if enabled).

    Shared by the trainers so every training run also emits sample plots.
    """
    if X.shape[0] == 0:
        print(f"  {name}: no events, skipped")
        return
    n = min(config.TRACE_N_SAMPLES, X.shape[0])
    sel = rng.choice(X.shape[0], n, replace=False)
    traces = X[sel]
    deriv = data.derivative(traces) if config.INCLUDE_DERIVATIVE else None
    # Show the exact model-input view (the derivative is shift-invariant).
    if config.SUBTRACT_BEAM_LEVEL:
        traces = traces - np.float32(1.0)
    mlplots.sample_traces(traces, name, label, deriv=deriv)
    print(f"  {name}: {n} traces drawn")


def main():
    rng = np.random.default_rng(config.SEED)
    print(f"music-ml: dataset {config.DATASET}, sample traces "
          f"(shorts={config.INCLUDE_SHORT_STRIPS} "
          f"guards={config.INCLUDE_GUARD_STRIPS} "
          f"deriv={config.INCLUDE_DERIVATIVE} ")

    specs = data.list_sim_specs(eres_only=True)
    gain = data.beam_gains(specs)
    strip = config.TRACE_STRIP
    labels = {
        "beam": "Beam (sim)",
        "aa": f"(#alpha,#alpha') s{strip} (sim)",
        "an": f"(#alpha,n) s{strip} (sim)",
    }
    for spec in specs:
        is_beam = spec.base == "beam" and spec.strip < 0
        is_reac = spec.base in ("aa", "an") and spec.strip == strip
        if not (is_beam or is_reac):
            continue
        name = ("sample_traces_beam"
                if is_beam else f"sample_traces_{spec.base}_s{strip}")
        draw_population(data.load_population(spec, gain), name,
                        labels[spec.base], rng)

    if not config.TRACE_DRAW_DATA:
        return
    try:
        totals = data.load_experimental_totals(max_files=1,
                                               max_events_per_file=100000)
    except FileNotFoundError as exc:
        print(f"skipping experimental samples: {exc}")
        return
    # The calibration already puts the data beam at 1 per strip: unit gains.
    unit = np.ones(totals.shape[1], dtype=np.float64)
    X, _ = data.max_normalize(totals, unit)
    draw_population(X, "sample_traces_data", "Data", rng)


if __name__ == "__main__":
    main()
