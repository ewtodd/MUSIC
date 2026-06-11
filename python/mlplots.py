"""ROOT-based diagnostic plots in the house PlottingUtils style."""

import numpy as np

import analysis_utilities as au

import config

_ROOT = None


def _root():
    """Lazy one-time ROOT + style init (batch mode, house dirs)."""
    global _ROOT
    if _ROOT is None:
        _ROOT = au.set_root_preferences(plots_dir=config.PLOTS_DIR,
                                        root_files_dir=config.ROOT_FILES_DIR)
    return _ROOT


def _graph(R, x, y):
    g = R.TGraph(len(x))
    for i in range(len(x)):
        g.SetPoint(i, float(x[i]), float(y[i]))
    return g


def loss_curves(train_losses,
                val_losses,
                name,
                subdir=config.PLOT_SUBDIR,
                train_label="Train",
                val_label="Validation"):
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kTRUE)
    epochs = np.arange(1, len(train_losses) + 1)
    g_train = _graph(R, epochs, train_losses)
    g_val = _graph(R, epochs, val_losses)
    R.PlottingUtils.ConfigureGraph(g_train, colors[0], ";Epoch;Loss")
    g_train.Draw("AL")
    R.PlottingUtils.ConfigureGraph(g_val, colors[1], "")
    g_val.Draw("L SAME")
    leg = R.PlottingUtils.AddLegend()
    leg.AddEntry(g_train, train_label, "l")
    leg.AddEntry(g_val, val_label, "l")
    leg.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLOG)


def score_hists(named_scores,
                name,
                x_title,
                subdir=config.PLOT_SUBDIR,
                n_bins=100):
    """Overlaid unit-area score distributions, one TH1F per population.

    named_scores: list of (label, np.ndarray) in draw/legend order.
    """
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    finite = np.concatenate(
        [s[np.isfinite(s)] for _, s in named_scores if s.size > 0])
    x_lo, x_hi = float(finite.min()), float(finite.max())
    pad = 0.05 * (x_hi - x_lo) if x_hi > x_lo else 1.0
    x_lo -= pad
    x_hi += pad

    c = R.PlottingUtils.GetConfiguredCanvas(R.kTRUE)
    hists = []
    y_max = 0.0
    for i, (label, scores) in enumerate(named_scores):
        h = R.TH1F(f"h_{name}_{i}", f";{x_title};Normalized counts", n_bins,
                   x_lo, x_hi)
        h.SetDirectory(R.nullptr)
        for v in scores:
            h.Fill(float(v))
        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())
        R.PlottingUtils.ConfigureHistogram(h, colors[i % colors.size()],
                                           f";{x_title};Normalized counts")
        h.SetStats(0)
        y_max = max(y_max, h.GetMaximum())
        hists.append((label, h))
    leg = R.PlottingUtils.AddLegend()
    for i, (label, h) in enumerate(hists):
        if i == 0:
            h.SetMaximum(3.0 * y_max)
            h.SetMinimum(1.0e-5)
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")
        leg.AddEntry(h, label, "l")
    leg.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLOG)


def roc_plot(fpr, tpr, auc_value, name, x_title, subdir=config.PLOT_SUBDIR):
    """Signal efficiency vs background rejection, AUC annotated."""
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    g = _graph(R, tpr, 1.0 - fpr)
    R.PlottingUtils.ConfigureGraph(g, colors[0],
                                   f";{x_title};Background rejection")
    g.GetXaxis().SetLimits(0.0, 1.05)
    g.GetHistogram().SetMinimum(0.0)
    g.GetHistogram().SetMaximum(1.05)
    g.Draw("AL")
    text = R.PlottingUtils.AddText(f"AUC = {auc_value:.4f}", 0.45, 0.3)
    text.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def _graph_errors(R, x, y, ey):
    g = R.TGraphErrors(len(x))
    for i in range(len(x)):
        g.SetPoint(i, float(x[i]), float(y[i]))
        g.SetPointError(i, 0.0, float(ey[i]))
    return g


def sweep_plot(x_values,
               means,
               errs,
               x_title,
               name,
               subdir=config.PLOT_SUBDIR,
               log_x=False):
    """AUC vs one swept hyperparameter, seed spread as error bars."""
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    if log_x:
        c.SetLogx(True)
    g = _graph_errors(R, x_values, means, errs)
    R.PlottingUtils.ConfigureGraph(g, colors[0], f";{x_title};ROC AUC")
    g.SetMarkerColor(colors[0])
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.2)
    y_lo = min(m - e for m, e in zip(means, errs))
    y_hi = max(m + e for m, e in zip(means, errs))
    pad = max(0.002, 0.25 * (y_hi - y_lo))
    g.GetHistogram().SetMinimum(y_lo - pad)
    g.GetHistogram().SetMaximum(y_hi + pad)
    g.Draw("APL")
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def labeled_points(labels,
                   means,
                   errs,
                   y_title,
                   name,
                   subdir=config.PLOT_SUBDIR):
    """AUC per labeled configuration (categorical x-axis with bin labels)."""
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    n = len(labels)
    frame = R.TH1F(f"frame_{name}", f";;{y_title}", n, 0.5, n + 0.5)
    frame.SetDirectory(R.nullptr)
    frame.SetStats(0)
    for i, label in enumerate(labels):
        frame.GetXaxis().SetBinLabel(i + 1, label)
    y_lo = min(m - e for m, e in zip(means, errs))
    y_hi = max(m + e for m, e in zip(means, errs))
    pad = max(0.002, 0.25 * (y_hi - y_lo))
    frame.SetMinimum(y_lo - pad)
    frame.SetMaximum(y_hi + pad)
    frame.Draw()
    x = np.arange(1, n + 1, dtype=np.float64)
    g = _graph_errors(R, x, means, errs)
    g.SetMarkerColor(colors[0])
    g.SetLineColor(colors[0])
    g.SetLineWidth(R.PlottingUtils.GetLineWidth())
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.3)
    g.Draw("PE SAME")
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def sample_traces(traces,
                  name,
                  label,
                  deriv=None,
                  y_title="#DeltaE [a.u.]",
                  subdir=config.PLOT_SUBDIR):
    """Overlay sampled per-strip traces.

    traces: (n, n_strip_inputs()) array in model feature units (i.e. after
    the configured flatten / max-norm), laid out per config.block_widths():
    the long-side trace, then (with INCLUDE_SHORT_STRIPS) the short-side
    block, drawn on a companion "<name>_short" canvas. deriv: optional
    block-local adjacent differences (data.derivative layout) for the same
    events, drawn on "<name>_deriv" (+ "_deriv_short") canvases with x at
    the strip midpoints. The x-axis always spans strips 0-17 so plots stay
    comparable across the guard-strip knob; the y-range is auto-scaled in
    flattened a.u..
    """
    R = _root()
    widths = config.block_widths()
    long_w = widths[0]
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    strip_x = np.arange(first, first + long_w, dtype=np.float64)
    short_x = np.arange(1, 17, dtype=np.float64)  # short ends: strips 1-16
    colors = [
        R.kRed + 2, R.kBlue + 2, R.kGreen + 2, R.kOrange + 2, R.kMagenta + 2
    ]

    def _draw_set(values, xs, out_name, ytitle, y_lo, y_hi):
        c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
        graphs = []
        leg = R.PlottingUtils.AddLegend()
        for i in range(values.shape[0]):
            g = _graph(R, xs, values[i].astype(np.float64))
            g.SetLineColor(colors[i % len(colors)])
            g.SetLineWidth(R.PlottingUtils.GetLineWidth())
            if i == 0:
                g.SetTitle("")
                g.GetXaxis().SetTitle("Strip")
                g.GetYaxis().SetTitle(ytitle)
                g.GetXaxis().SetLimits(-0.5, 17.5)
                g.GetHistogram().SetMinimum(y_lo)
                g.GetHistogram().SetMaximum(y_hi)
                g.Draw("AL")
            else:
                g.Draw("L SAME")
            graphs.append(g)
            leg.AddEntry(g, f"Event {i + 1}", "l")
        leg.Draw()
        text = R.PlottingUtils.AddText(label, 0.42, 0.85)
        text.Draw()
        R.PlottingUtils.SaveFigure(c, out_name, subdir,
                                   R.PlotSaveOptions.kLINEAR)

    beam_level = 0.0 if config.SUBTRACT_BEAM_LEVEL else 1.0
    peak = float(np.max(traces)) if traces.size else beam_level
    y_lo = beam_level - 0.2
    y_hi = max(1.2 * peak, beam_level + 0.2)
    _draw_set(traces[:, :long_w], strip_x, name, y_title, y_lo, y_hi)
    if len(widths) > 1:
        _draw_set(traces[:, long_w:], short_x, f"{name}_short",
                  f"Short-end {y_title}", y_lo, y_hi)

    if deriv is not None:
        d_peak = float(np.max(np.abs(deriv))) if deriv.size else 1.0
        d_peak = max(d_peak, 0.1)
        _draw_set(deriv[:, :long_w - 1], strip_x[:-1] + 0.5, f"{name}_deriv",
                  "Adjacent-strip #Delta(#DeltaE) [a.u.]", -1.2 * d_peak,
                  1.2 * d_peak)
        if len(widths) > 1:
            _draw_set(deriv[:, long_w - 1:], short_x[:-1] + 0.5,
                      f"{name}_deriv_short",
                      "Short-end adjacent-strip #Delta(#DeltaE) [a.u.]",
                      -1.2 * d_peak, 1.2 * d_peak)


def shap_vs_strip(strip_x,
                  trace_shap,
                  trace_err,
                  avg_reacted,
                  name,
                  deriv_x=None,
                  deriv_shap=None,
                  deriv_err=None,
                  subdir=config.PLOT_SUBDIR):
    """Mean |SHAP| per input feature vs strip index.

    The trace-block (and optional derivative-block) importances are each
    normalized to their own max; the average reacted trace (max-normalized)
    is overlaid in gray for orientation.
    """
    R = _root()
    colors = R.PlottingUtils.GetDefaultColors()
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)

    g_trace = _graph(R, strip_x, avg_reacted / max(avg_reacted.max(), 1e-12))
    R.PlottingUtils.ConfigureGraph(g_trace, R.kGray + 2,
                                   ";Strip;Normalized value [a.u.]")
    g_trace.GetHistogram().SetMinimum(-0.05)
    g_trace.GetHistogram().SetMaximum(1.3)
    g_trace.Draw("AL")

    peak = max(float(np.max(trace_shap)), 1e-12)
    g_shap = _graph_errors(R, strip_x, trace_shap / peak, trace_err / peak)
    g_shap.SetLineColor(colors[0])
    g_shap.SetMarkerColor(colors[0])
    g_shap.SetLineWidth(R.PlottingUtils.GetLineWidth())
    g_shap.SetMarkerStyle(20)
    g_shap.SetMarkerSize(1.1)
    g_shap.Draw("PL SAME")

    leg = R.PlottingUtils.AddLegend()
    leg.AddEntry(g_trace, "Average reacted trace", "l")
    leg.AddEntry(g_shap, "|SHAP| (trace)", "lp")

    g_dshap = None
    if deriv_shap is not None:
        dpeak = max(float(np.max(deriv_shap)), 1e-12)
        g_dshap = _graph_errors(R, deriv_x, deriv_shap / dpeak,
                                deriv_err / dpeak)
        g_dshap.SetLineColor(colors[1])
        g_dshap.SetMarkerColor(colors[1])
        g_dshap.SetLineWidth(R.PlottingUtils.GetLineWidth())
        g_dshap.SetMarkerStyle(21)
        g_dshap.SetMarkerSize(1.1)
        g_dshap.Draw("PL SAME")
        leg.AddEntry(g_dshap, "|SHAP| (derivative)", "lp")

    leg.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)
