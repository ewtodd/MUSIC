"""ROOT-based plots in the house PlottingUtils style (blind pipeline)."""

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


def cluster_scatter(named_points, name, x_title, y_title,
                    subdir=config.PLOT_SUBDIR, max_points=20000):
    """2-D feature scatter, one marker color per cluster.

    named_points: list of (label, x, y, cidx) in draw/legend order. The colour
    and marker style are keyed to `cidx` (the cluster index), NOT the draw
    order, so a given cluster keeps its colour across every plot even when
    other clusters are empty/pruned. Each cluster is subsampled to max_points
    so the canvas stays light; axes span the data with padding.
    """
    R = _root()
    palette = _trace_palette(R)
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    all_x = np.concatenate([np.asarray(p[1], dtype=np.float64)
                            for p in named_points if len(p[1]) > 0])
    all_y = np.concatenate([np.asarray(p[2], dtype=np.float64)
                            for p in named_points if len(p[2]) > 0])
    x_lo, x_hi = float(all_x.min()), float(all_x.max())
    y_lo, y_hi = float(all_y.min()), float(all_y.max())
    xpad = 0.05 * (x_hi - x_lo) if x_hi > x_lo else 1.0
    ypad = 0.05 * (y_hi - y_lo) if y_hi > y_lo else 1.0
    leg = R.PlottingUtils.AddLegend()
    graphs = []
    rng = np.random.default_rng(config.SEED)
    for i, (label, x, y, cidx) in enumerate(named_points):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if x.size > max_points:
            sel = rng.choice(x.size, max_points, replace=False)
            x, y = x[sel], y[sel]
        g = _graph(R, x, y) if x.size else R.TGraph()
        g.SetMarkerColor(palette[cidx % len(palette)])
        g.SetMarkerStyle(20 + cidx)
        g.SetMarkerSize(0.5)
        if i == 0:
            g.SetTitle("")
            g.GetXaxis().SetTitle(x_title)
            g.GetYaxis().SetTitle(y_title)
            g.GetXaxis().SetLimits(x_lo - xpad, x_hi + xpad)
            g.GetHistogram().SetMinimum(y_lo - ypad)
            g.GetHistogram().SetMaximum(y_hi + ypad)
            g.Draw("AP")
        else:
            g.Draw("P SAME")
        graphs.append(g)
        leg.AddEntry(g, label, "p")
    leg.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def feature_hist2d(named_points, name, x_title, y_title,
                   subdir=config.PLOT_SUBDIR, nbins=60,
                   max_points=config.BLIND_FEAT2D_MAXPTS):
    """Overlaid per-cluster 2-D density: one TH2 per cluster, filled with that
    cluster's (x, y) and drawn as a translucent BOX whose colour is keyed to the
    cluster index `cidx` (NOT draw order) via the same palette as the trace /
    mean overlays -- so a cluster keeps its colour across every plot. Noise
    (cidx < 0) draws grey. BOX self-normalizes each cluster, so a small cluster
    still shows where it sits.

    named_points: list of (label, x, y, cidx). Shared axes span the 0.5-99.5
    percentile of the pooled data (outlier-robust). Each cluster is subsampled
    to max_points before filling. Saved linear."""
    R = _root()
    palette = _trace_palette(R)
    pts = [p for p in named_points if len(p[1]) > 0]
    if not pts:
        return
    all_x = np.concatenate([np.asarray(p[1], dtype=np.float64) for p in pts])
    all_y = np.concatenate([np.asarray(p[2], dtype=np.float64) for p in pts])
    x_lo, x_hi = np.percentile(all_x, [0.5, 99.5])
    y_lo, y_hi = np.percentile(all_y, [0.5, 99.5])
    if not x_hi > x_lo:
        x_hi = x_lo + 1.0
    if not y_hi > y_lo:
        y_hi = y_lo + 1.0
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    frame = R.TH2F(f"f2d_{name}", f";{x_title};{y_title}", nbins, float(x_lo),
                   float(x_hi), nbins, float(y_lo), float(y_hi))
    frame.SetStats(0)
    frame.SetDirectory(0)
    frame.Draw()
    leg = R.PlottingUtils.AddLegend()
    keep = [frame]
    rng = np.random.default_rng(config.SEED)
    # Noise first (drawn behind), then clusters on top.
    order = sorted(range(len(pts)), key=lambda i: pts[i][3] >= 0)
    for i in order:
        label, x, y, cidx = pts[i]
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if max_points is not None and x.size > max_points:
            sel = rng.choice(x.size, max_points, replace=False)
            x, y = x[sel], y[sel]
        h = R.TH2F(f"h2d_{name}_{cidx}", "", nbins, float(x_lo), float(x_hi),
                   nbins, float(y_lo), float(y_hi))
        h.SetDirectory(0)
        h.SetStats(0)
        h.FillN(int(x.size), x, y, np.ones(x.size, dtype=np.float64))
        col = R.kGray + 1 if cidx < 0 else palette[cidx % len(palette)]
        h.SetLineColor(col)
        h.SetLineWidth(1)
        h.SetFillColorAlpha(col, 0.35)
        h.SetMarkerColor(col)
        h.Draw("BOX SAME")
        keep.append(h)
        leg.AddEntry(h, label, "f")
    leg.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def value_hist(values, name, x_title, subdir=config.PLOT_SUBDIR, vline=None,
               label="", nbins=80):
    """1-D histogram of a scalar per-event quantity (e.g. a template-match
    residual), with an optional dashed red vertical line at `vline` (the cut).
    Saved linear."""
    R = _root()
    v = np.asarray(values, dtype=np.float64)
    if v.size == 0:
        return
    lo, hi = float(v.min()), float(v.max())
    if not hi > lo:
        hi = lo + 1.0
    pad = 0.05 * (hi - lo)
    h = R.TH1F(f"h_{name}", f";{x_title};counts", nbins, lo - pad, hi + pad)
    h.SetDirectory(0)
    h.SetStats(0)
    h.FillN(int(v.size), v, np.ones(v.size, dtype=np.float64))
    h.SetLineColor(R.kAzure + 2)
    h.SetLineWidth(2)
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    h.Draw("HIST")
    keep = [h]
    if vline is not None:
        ln = R.TLine(float(vline), 0.0, float(vline), h.GetMaximum())
        ln.SetLineColor(R.kRed + 1)
        ln.SetLineWidth(2)
        ln.SetLineStyle(2)
        ln.Draw("SAME")
        keep.append(ln)
    if label:
        text = R.PlottingUtils.AddText(label, 0.42, 0.85)
        text.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def _trace_palette(R):
    """Per-cluster line colours, leading with the StripSumScatter
    DrawRegionTraces palette (beam grey, (a,a') azure, (a,n) red) and
    extended for k > 3."""
    return [R.kGray + 2, R.kAzure + 2, R.kRed + 1,
            R.kGreen + 2, R.kOrange + 2, R.kMagenta + 2, R.kCyan + 2]


def _frame_and_canvas(R, name, y_title):
    """A StripSumScatter-style framed linear canvas: TH2F frame (18 x-bins
    -0.5..17.5, y fixed STRIP_DE_MIN/MAX_NORMED = 0.8..1.3) drawn on a fresh
    configured canvas. Returns (canvas, frame, strip_x)."""
    long_w = config.block_widths()[0]
    first = 0 if config.INCLUDE_GUARD_STRIPS else 1
    strip_x = np.arange(first, first + long_w, dtype=np.float64)
    frame = R.TH2F(f"frame_{name}", f";Strip;{y_title}", 18, -0.5, 17.5,
                   100, 0.8, 1.3)
    frame.SetStats(0)
    frame.SetDirectory(0)
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    frame.Draw()
    return c, frame, strip_x


def _beam_ref_graph(R, strip_x, reference, long_w):
    """Dashed black beam-reference TGraph drawn on the current canvas, or None
    when no reference is supplied."""
    if reference is None:
        return None
    g_ref = _graph(R, strip_x, np.asarray(reference)[:long_w])
    g_ref.SetLineColor(R.kBlack)
    g_ref.SetLineWidth(2 * R.PlottingUtils.GetLineWidth())
    g_ref.SetLineStyle(2)  # dashed
    g_ref.Draw("L SAME")
    return g_ref


def overlay_cluster_means(groups,
                          name,
                          label,
                          y_title="#DeltaE [a.u.]",
                          subdir=config.PLOT_SUBDIR,
                          reference=None,
                          bands=True):
    """Overlay each cluster's MEAN per-strip trace, with a +-1 sigma band, on
    ONE StripSumScatter-style canvas (fixed 0.8..1.3): a translucent band per
    cluster (drawn behind), then one bold (width-3) mean line coloured by the
    cluster's `cidx` (NOT draw order), the beam reference dashed on top.
    `groups` is a list of (cluster_label, mean_trace (long_w,),
    sigma_trace (long_w,), cidx); cidx keys the colour so a cluster keeps it
    across every plot. bands=False draws the mean lines ONLY (no sigma shading),
    for many-cluster plots where overlapping bands wash the canvas out."""
    R = _root()
    long_w = config.block_widths()[0]
    c, frame, strip_x = _frame_and_canvas(R, name, y_title)
    palette = _trace_palette(R)
    keep = [frame]
    # Bands first so the mean lines sit on top of every band.
    if bands:
        for gi in range(len(groups)):
            mean = np.asarray(groups[gi][1], dtype=np.float64)[:long_w]
            sigma = np.asarray(groups[gi][2], dtype=np.float64)[:long_w]
            cidx = groups[gi][3]
            band = R.TGraph(2 * long_w)
            for i in range(long_w):
                band.SetPoint(i, float(strip_x[i]), float(mean[i] + sigma[i]))
                band.SetPoint(2 * long_w - 1 - i, float(strip_x[i]),
                              float(mean[i] - sigma[i]))
            band.SetFillColorAlpha(palette[cidx % len(palette)], 0.15)
            band.SetLineWidth(0)
            band.Draw("F SAME")
            keep.append(band)
    entries = []
    for gi in range(len(groups)):
        clabel, mean, _, cidx = groups[gi]
        g = _graph(R, strip_x, np.asarray(mean, dtype=np.float64)[:long_w])
        g.SetLineColor(palette[cidx % len(palette)])
        g.SetLineWidth(3)
        g.Draw("L SAME")
        keep.append(g)
        entries.append((g, clabel))
    g_ref = _beam_ref_graph(R, strip_x, reference, long_w)
    leg = R.PlottingUtils.AddLegend(0.725, 0.875, 0.70, 0.86)
    for gi in range(len(entries)):
        leg.AddEntry(entries[gi][0], entries[gi][1], "l")
    if g_ref is not None:
        leg.AddEntry(g_ref, "beam", "l")
    leg.Draw()
    text = R.PlottingUtils.AddText(label, 0.42, 0.85)
    text.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def overlay_cluster_traces(groups,
                           name,
                           label,
                           y_title="#DeltaE [a.u.]",
                           subdir=config.PLOT_SUBDIR,
                           reference=None):
    """Overlay EVERY cluster's individual per-strip traces on ONE canvas, in
    the StripSumScatter DrawRegionTraces style: opaque width-1 lines coloured
    by the cluster's `cidx` (NOT draw order), a width-3 proxy legend, the beam
    reference dashed on top. `groups` is a list of
    (cluster_label, traces (m, long_w), cidx)."""
    R = _root()
    long_w = config.block_widths()[0]
    c, frame, strip_x = _frame_and_canvas(R, name, y_title)
    palette = _trace_palette(R)
    keep = [frame]
    proxies = []
    for gi in range(len(groups)):
        clabel, traces, cidx = groups[gi]
        color = palette[cidx % len(palette)]
        vals = traces[:, :long_w]
        for i in range(vals.shape[0]):
            g = _graph(R, strip_x, vals[i].astype(np.float64))
            g.SetLineColor(color)  # opaque, like DrawRegionTraces
            g.SetLineWidth(1)
            g.Draw("L SAME")
            keep.append(g)
        p = R.TGraph(1)
        p.SetPoint(0, -1.0e9, -1.0e9)
        p.SetLineColor(color)
        p.SetLineWidth(3)
        keep.append(p)
        proxies.append((p, clabel))
    g_ref = _beam_ref_graph(R, strip_x, reference, long_w)
    leg = R.PlottingUtils.AddLegend(0.725, 0.875, 0.70, 0.86)
    for gi in range(len(proxies)):
        leg.AddEntry(proxies[gi][0], proxies[gi][1], "l")
    if g_ref is not None:
        leg.AddEntry(g_ref, "beam", "l")
    leg.Draw()
    text = R.PlottingUtils.AddText(label, 0.42, 0.85)
    text.Draw()
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def scree_plot(ratios, name, subdir=config.PLOT_SUBDIR, label=""):
    """PCA scree: per-PC explained-variance ratio (azure) and the running
    cumulative (red dashed) vs component index. A few PCs reaching ~1.0 means
    the feature set is largely redundant."""
    R = _root()
    ratios = np.asarray(ratios, dtype=np.float64)
    n = ratios.shape[0]
    if n == 0:
        return
    cum = np.cumsum(ratios)
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    g = R.TGraph(n)
    gc = R.TGraph(n)
    for i in range(n):
        g.SetPoint(i, i + 1, float(ratios[i]))
        gc.SetPoint(i, i + 1, float(cum[i]))
    g.SetTitle("")
    g.GetXaxis().SetTitle("principal component")
    g.GetYaxis().SetTitle("explained variance ratio")
    g.GetHistogram().SetMinimum(0.0)
    g.GetHistogram().SetMaximum(1.05)
    g.SetLineColor(R.kAzure + 2)
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(R.kAzure + 2)
    g.Draw("ALP")
    gc.SetLineColor(R.kRed + 1)
    gc.SetLineWidth(3)
    gc.SetLineStyle(2)
    gc.SetMarkerStyle(21)
    gc.SetMarkerColor(R.kRed + 1)
    gc.Draw("LP SAME")
    leg = R.PlottingUtils.AddLegend()
    leg.AddEntry(g, "per-PC", "lp")
    leg.AddEntry(gc, "cumulative", "lp")
    leg.Draw()
    keep = [g, gc, leg]
    if label:
        text = R.PlottingUtils.AddText(label, 0.42, 0.45)
        text.Draw()
        keep.append(text)
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)


def matrix_hist(matrix, row_labels, col_labels, name, x_title, y_title,
                subdir=config.PLOT_SUBDIR, z_title=""):
    """COLZ heatmap of a (rows x cols) matrix with per-bin axis labels -- used
    for PCA loadings (PC x feature) and the feature correlation matrix. Rows map
    to the y-axis (bottom-up), cols to the x-axis."""
    R = _root()
    M = np.asarray(matrix, dtype=np.float64)
    nr, nc = M.shape
    h = R.TH2F(f"m_{name}", f";{x_title};{y_title}", nc, 0.0, float(nc), nr,
               0.0, float(nr))
    h.SetStats(0)
    h.SetDirectory(0)
    for i in range(nr):
        for j in range(nc):
            h.SetBinContent(j + 1, i + 1, float(M[i, j]))
    for j in range(nc):
        h.GetXaxis().SetBinLabel(j + 1, str(col_labels[j]))
    for i in range(nr):
        h.GetYaxis().SetBinLabel(i + 1, str(row_labels[i]))
    if z_title:
        h.GetZaxis().SetTitle(z_title)
    h.SetMarkerSize(1.2)
    c = R.PlottingUtils.GetConfiguredCanvas(R.kFALSE)
    h.Draw("COLZ TEXT")
    R.PlottingUtils.SaveFigure(c, name, subdir, R.PlotSaveOptions.kLINEAR)
