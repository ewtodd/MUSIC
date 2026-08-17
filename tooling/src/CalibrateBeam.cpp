#include "CalibrateBeamInternal.hpp"

const Int_t kMaxChannels = 35;

std::vector<ChannelCal> CalibrateBeam::BuildChannels() {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  c.name = "Strip0";
  c.side = 'S';
  c.strip = 0;
  chans.push_back(c);
  c.name = "Strip17";
  c.side = 'S';
  c.strip = 17;
  chans.push_back(c);
  for (Int_t s = 1; s <= 16; s++) {
    c.name = Form("L%d", s);
    c.side = 'L';
    c.strip = s;
    chans.push_back(c);
    c.name = Form("R%d", s);
    c.side = 'R';
    c.strip = s;
    chans.push_back(c);
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  chans.push_back(c);
  return chans;
}

Int_t ChannelToEresIndex(const ChannelCal &c) {
  if (c.side == 'C')
    return 0;
  if (c.side == 'S' && c.strip == 0)
    return 1;
  if (c.side == 'S' && c.strip == 17)
    return 2;
  if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
    return 3 + (c.strip - 1);
  if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
    return 19 + (c.strip - 1);
  return -1;
}

// Per-channel ADC histogram of the samples that fed each beam anchor. One file
// per channel under <plot_subdir>/beam_peak, named beam_peak_<channel>.
void SaveBeamPeakChannelHistograms(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<Float_t>> &samples,
    const std::vector<TF1 *> &fits, const TString &plot_subdir) {
  TString subdir = plot_subdir + "/beam_peak";
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    const std::vector<Float_t> &v = samples[i];
    if (Long64_t(v.size()) < kMinSamples)
      continue;
    // Same binning recipe the fit used, so the overlaid Gaussian's amplitude
    // matches this histogram exactly.
    Double_t mode = 0.0, sigma = 0.0;
    RobustPeakSeed(v, mode, sigma);
    TH1F *h = MakeBeamPeakHist(Form("h_beam_peak_%s", c.name.Data()),
                               Form(";%s #DeltaE [ADC];Counts", c.name.Data()),
                               v, mode, sigma);
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h, kBlack);
    TF1 *fit = fits[i];
    if (fit) {
      fit->SetLineColor(kViolet + 2);
      fit->SetLineWidth(2);
      fit->Draw("L SAME");
    }
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, Form("beam_peak_%s", c.name.Data()), subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete h;
  }
}

// Writes a one-row `calibration` tree into the open file `dst` (typically a
// per-subfile .cal.root). Layout matches AggregateEresTomlForRun's reader.
// `align` carries the beam-energy window and per-strip alignment factors;
// pass nullptr when neither has been computed yet (branches are written as
// zero).
void WriteCalibrationTree(TFile *dst, const std::vector<ChannelCal> &chans,
                          const StripAlignmentResult *align) {
  dst->cd();
  if (TObject *old = dst->Get("calibration"))
    old->Delete();
  TTree *cal = new TTree("calibration", "Per-channel normMUSIC calibration");
  Float_t gain[kMaxChannels] = {0};
  Float_t fit_adc[kMaxChannels] = {0}, fit_sigma[kMaxChannels] = {0};
  Long64_t fit_n[kMaxChannels] = {0};
  Bool_t ok[kMaxChannels] = {0};
  // Per-strip gains laid out to match the events tree exactly: GainLeft[s]
  // multiplies Left_0_17_dE[s] (s=0/17 are the single-ended guards, s=1..16 the
  // left ends), GainRight[s] multiplies RightdE[s] (0 at the guards). This is
  // what EnergyView reads to calibrate on the fly -- no per-event a.u. is
  // stored.
  Float_t gain_left[18] = {0}, gain_right[18] = {0};
  Float_t gain_cathode = 0.0f;
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    ok[k] = IsCalibrated(c) || c.gain > 0;
    gain[k] = ok[k] ? Float_t(Gain(c)) : 0.0f;
    fit_adc[k] = Float_t(c.fit_adc);
    fit_sigma[k] = Float_t(c.fit_sigma_adc);
    fit_n[k] = c.n_samples;
    if (c.side == 'S' && c.strip >= 0 && c.strip <= 17)
      gain_left[c.strip] = gain[k];
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
      gain_left[c.strip] = gain[k];
    else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
      gain_right[c.strip] = gain[k];
    else if (c.side == 'C')
      gain_cathode = gain[k];
  }
  cal->Branch("Gain", gain, Form("Gain[%d]/F", kMaxChannels));
  cal->Branch("Ok", ok, Form("Ok[%d]/O", kMaxChannels));
  cal->Branch("FitADC", fit_adc, Form("FitADC[%d]/F", kMaxChannels));
  cal->Branch("FitSigmaADC", fit_sigma,
              Form("FitSigmaADC[%d]/F", kMaxChannels));
  cal->Branch("FitN", fit_n, Form("FitN[%d]/L", kMaxChannels));
  cal->Branch("GainLeft", gain_left, "GainLeft[18]/F");
  cal->Branch("GainRight", gain_right, "GainRight[18]/F");
  cal->Branch("GainCathode", &gain_cathode, "GainCathode/F");

  // Beam-energy window and per-strip multiplicative alignment factors
  // matching the notebook approach (pol3 reference / centroid).
  // EnergyView applies total_corrected = factor * total after the
  // per-channel gain. Default factor = 1.0 (identity).
  Float_t beam_e_min = 0.0f, beam_e_max = 0.0f;
  Float_t strip_factor[18];
  for (Int_t s = 0; s < 18; s++)
    strip_factor[s] = 1.0f;
  if (align) {
    beam_e_min = Float_t(align->beam_e_min);
    beam_e_max = Float_t(align->beam_e_max);
    if (align->ok) {
      for (Int_t s = 0; s < 18; s++)
        strip_factor[s] = Float_t(align->factors[s]);
    }
  }
  cal->Branch("BeamEMin", &beam_e_min, "BeamEMin/F");
  cal->Branch("BeamEMax", &beam_e_max, "BeamEMax/F");
  cal->Branch("StripFactor", strip_factor, "StripFactor[18]/F");
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
}

// Writes the per-channel gain table (tree "calibration") into the subfile's own
// events file. No per-event calibrated tree is produced: downstream readers
// recover a.u. on the fly via gain x raw ADC (EnergyView), so the raw events
// tree plus this one-row gain table fully determine every calibrated value.
void WriteCalibrationToEvents(const FileSpec &spec,
                              const std::vector<ChannelCal> &chans,
                              const StripAlignmentResult *align) {
  TString events_subpath = FileSet::EventsName(spec) + ".root";
  TFile *f = IO::OpenForWriting(events_subpath, "UPDATE");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << " to write calibration"
              << std::endl;
    if (f)
      delete f;
    return;
  }
  WriteCalibrationTree(f, chans, align);
  std::cout << "  wrote calibration into " << events_subpath << std::endl;
  f->Close();
  delete f;
}

// Derives the beam-energy window from the Strip0 (entrance guard) beam-peak
// fit. The notebook (37Cl_an.ipynb cell 10) fits a Gaussian to raw stp0 ADC
// and takes mu ± 3*sigma. Here, Strip0's beam peak is already fitted in
// ReduceToAnchors (fit_adc / fit_sigma_adc), so we convert to a.u. via the
// channel's own gain. In a.u. the peak sits at 1.0 by construction, so the
// window is 1.0 ± 3*sigma/fit_adc.
void DeriveBeamEnergyWindow(const std::vector<ChannelCal> &chans,
                            StripAlignmentResult &align) {
  const Double_t kBeamNSigma = 3.0;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (c.side == 'S' && c.strip == 0 && IsCalibrated(c)) {
      Double_t g = Gain(c);
      align.beam_e_min = g * (c.fit_adc - kBeamNSigma * c.fit_sigma_adc);
      align.beam_e_max = g * (c.fit_adc + kBeamNSigma * c.fit_sigma_adc);
      std::cout << "  beam energy window (from Strip0): [" << align.beam_e_min
                << ", " << align.beam_e_max << "] a.u." << std::endl;
      return;
    }
  }
  std::cerr << "  beam energy window: Strip0 not calibrated, using [0, 0]"
            << std::endl;
}

// Per-strip multiplicative alignment, matching the notebook
// (37Cl_an.ipynb cell 5). Decodes events with the per-channel gains already on
// disk, finds each strip's beam-peak centroid from the eSum 2D histogram,
// fits a robust degree-3 polynomial reference trend through the centroids
// (strips 1-16 only — the single-ended anodes sit at a different scale), and
// derives a multiplicative factor = reference[s] / centroid[s] for 1-16.
// Strips 0/17 get factor = 1.0 / centroid (push beam peak to 1.0 a.u.).
//
// Uses ALL events (not beam-gated): the beam dominates every strip's
// histogram by a wide margin.
StripAlignmentResult FindStripCentroidAlignment(const FileSpec &spec,
                                                const TString &plot_subdir,
                                                const TString &file_label) {
  const Int_t kNStrips = 18;
  const Int_t kNHistBins = 200;
  const Double_t kHistMin = 0.0;
  const Double_t kHistMax = 10.0;
  const Int_t kSmoothTimes = 5;
  const Double_t kMisalignPct = 1.5;
  const Int_t kMaxIter = 4;
  const Int_t kPolyDeg = 3;
  const Double_t kGausFitHalfWidth = 0.15;

  StripAlignmentResult result;
  for (Int_t s = 0; s < kNStrips; s++) {
    result.factors[s] = 1.0;
    result.centroids[s] = 0.0;
  }

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    delete sf;
    return result;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return result;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    std::cerr << "  " << file_label
              << ": calibration tree not found -- skipping alignment"
              << std::endl;
    sf->Close();
    delete sf;
    return result;
  }

  TH2D *h2 = new TH2D(Form("h2_strip_align_%s", file_label.Data()),
                      ";Strip number;#DeltaE [a.u.]", kNStrips, -0.5,
                      kNStrips - 0.5, kNHistBins, kHistMin, kHistMax);
  h2->SetDirectory(nullptr);

  Long64_t n = tree->GetEntries();
  Long64_t n_used = 0;
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    Bool_t any = kFALSE;
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.total[s];
      if (v <= 0)
        continue;
      h2->Fill(Double_t(s), v);
      any = kTRUE;
    }
    if (any)
      n_used++;
  }
  sf->Close();
  delete sf;

  std::cout << "  strip alignment: " << n_used << " events decoded"
            << std::endl;

  Double_t beam_centroids[kNStrips] = {0};
  Bool_t beam_ok[kNStrips] = {kFALSE};

  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t bin_ix = s + 1;
    TH1D *proj = h2->ProjectionY(
        Form("hproj_align_%s_s%d", file_label.Data(), s), bin_ix, bin_ix);
    proj->SetDirectory(nullptr);
    Long64_t n_entries = Long64_t(proj->GetEntries());
    if (n_entries < kMinSamples) {
      std::cerr << "  strip " << s << ": too few entries for alignment ("
                << n_entries << ")" << std::endl;
      delete proj;
      continue;
    }
    proj->Smooth(kSmoothTimes);

    // Strip 0/17: single-ended anode, beam at ~1.0, pileup at ~2.0.
    // Strips 1-16: total is bimodal when the short side doesn't fire
    // (long-only ≈ 0.5) — pick the peak nearest 1.0, which is the
    // full-strip beam peak (both sides contributing, eSum ≈ 1.0).
    Double_t peak_target = 1.0;
    Double_t search_lo, search_hi;
    if (s == 0 || s == 17) {
      search_lo = 0.3;
      search_hi = 1.6;
    } else {
      search_lo = 0.6;
      search_hi = 1.5;
    }
    Int_t b_lo = proj->FindBin(search_lo);
    Int_t b_hi = proj->FindBin(search_hi);
    Int_t b_max = -1;
    Double_t best_dist = 1e9;
    for (Int_t b = b_lo; b <= b_hi; b++) {
      Double_t v = proj->GetBinContent(b);
      if (v <= 0)
        continue;
      // Local maximum check: higher than neighbours
      if (b > b_lo && proj->GetBinContent(b - 1) >= v)
        continue;
      if (b < b_hi && proj->GetBinContent(b + 1) > v)
        continue;
      Double_t bc = proj->GetBinCenter(b);
      Double_t d = TMath::Abs(bc - peak_target);
      if (d < best_dist) {
        best_dist = d;
        b_max = b;
      }
    }
    // Fallback: global maximum if no local peaks found near 1.0
    if (b_max < 0) {
      Double_t val_max = 0;
      for (Int_t b = b_lo; b <= b_hi; b++) {
        Double_t v = proj->GetBinContent(b);
        if (v > val_max) {
          val_max = v;
          b_max = b;
        }
      }
    }
    if (b_max < 0) {
      delete proj;
      continue;
    }
    Double_t seed_peak = proj->GetBinCenter(b_max);
    Double_t seed_h = proj->GetBinContent(b_max);
    if (seed_h <= 0 || seed_peak <= 0) {
      delete proj;
      continue;
    }

    // Sub-bin precision: Gaussian fit around the smoothed max-bin seed
    Double_t fit_lo = seed_peak - kGausFitHalfWidth;
    Double_t fit_hi = seed_peak + kGausFitHalfWidth;
    if (fit_lo < kHistMin)
      fit_lo = kHistMin;
    if (fit_hi > kHistMax)
      fit_hi = kHistMax;
    TF1 *fg = new TF1("f_align_peak_refine", "gaus", fit_lo, fit_hi);
    fg->SetParameters(seed_h, seed_peak, 0.05);
    fg->SetParLimits(1, fit_lo, fit_hi);
    TFitResultPtr r = proj->Fit(fg, "QSRN");
    Double_t beam_peak = fg->GetParameter(1);
    if (!(r.Get() && r->IsValid() && beam_peak > 0))
      beam_peak = seed_peak;
    delete fg;

    beam_centroids[s] = beam_peak;
    beam_ok[s] = kTRUE;
    std::cout << "  strip " << s << " beam=" << Form("%.4f", beam_peak)
              << " a.u.  (n=" << n_entries << ")" << std::endl;
    delete proj;
  }

  // Robust pol3 reference trend through strips 1-16 centroids.
  // Iteratively drop the strip with the worst residual > kMisalignPct.
  TGraph *g_cent = new TGraph(kNStrips);
  Int_t np = 0;
  for (Int_t s = 1; s <= 16; s++) {
    if (beam_ok[s]) {
      g_cent->SetPoint(np, Double_t(s), beam_centroids[s]);
      np++;
    }
  }
  g_cent->Set(np);

  std::set<Int_t> outliers;
  for (Int_t iter = 0; iter < kMaxIter; iter++) {
    if (g_cent->GetN() <= kPolyDeg + 1)
      break;
    TF1 *fpol = new TF1(Form("fpol_align_%s_iter%d", file_label.Data(), iter),
                        "pol3", -0.5, kNStrips - 0.5);
    TFitResultPtr r = g_cent->Fit(fpol, "QSRN");
    if (!r.Get() || !r->IsValid()) {
      delete fpol;
      break;
    }
    Double_t worst_pct = 0;
    Int_t worst_idx = -1;
    for (Int_t i = 0; i < g_cent->GetN(); i++) {
      Double_t x, y;
      g_cent->GetPoint(i, x, y);
      Double_t pred = fpol->Eval(x);
      if (pred <= 0)
        continue;
      Double_t resid_pct = TMath::Abs(y - pred) / pred * 100.0;
      if (resid_pct > worst_pct) {
        worst_pct = resid_pct;
        worst_idx = i;
      }
    }
    delete fpol;
    if (worst_pct < kMisalignPct || worst_idx < 0)
      break;
    Double_t rx, ry;
    g_cent->GetPoint(worst_idx, rx, ry);
    std::cout << "  strip " << Int_t(rx) << " outlier "
              << Form("%.1f", worst_pct) << "% — dropped" << std::endl;
    outliers.insert(Int_t(rx));
    g_cent->RemovePoint(worst_idx);
  }

  TF1 *fbeam = nullptr;
  if (g_cent->GetN() > kPolyDeg + 1) {
    fbeam = new TF1(Form("fbeam_align_%s", file_label.Data()), "pol3", -0.5,
                    kNStrips - 0.5);
    g_cent->Fit(fbeam, "QSRN");
    std::cout << "  strip alignment reference: pol3 fitted through "
              << g_cent->GetN() << " strips" << std::endl;
  }

  // Notebook (37Cl_an.ipynb cell 5): factors = reference / centroid,
  // where reference is the pol3 trend through strips 1-16 centroids.
  // Strips 0/17 are single-ended anodes not in the pol3 fit; push to 1.0.
  Int_t valid_strips = 0;
  for (Int_t s = 0; s < kNStrips; s++) {
    if (!beam_ok[s])
      continue;
    Double_t centro = beam_centroids[s];
    if (centro <= 0)
      continue;
    result.centroids[s] = centro;
    if (s >= 1 && s <= 16 && fbeam) {
      Double_t ref = fbeam->Eval(Double_t(s));
      if (ref > 0)
        result.factors[s] = ref / centro;
      else
        result.factors[s] = 1.0 / centro;
    } else {
      result.factors[s] = 1.0 / centro;
    }
    valid_strips++;
    std::cout << "  strip " << s << " centroid=" << Form("%.4f", centro)
              << " factor=" << Form("%.4f", result.factors[s]) << std::endl;
  }
  result.ok = (valid_strips >= 4) ? kTRUE : kFALSE;

  // Diagnostic plot
  {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2, cv);
    TGraph *g_beam_plot = new TGraph(kNStrips);
    Int_t nb = 0;
    for (Int_t s = 0; s < kNStrips; s++) {
      if (beam_ok[s]) {
        g_beam_plot->SetPoint(nb, Double_t(s), beam_centroids[s]);
        nb++;
      }
    }
    g_beam_plot->Set(nb);
    if (nb > 0) {
      g_beam_plot->SetMarkerStyle(20);
      g_beam_plot->SetMarkerColor(kOrange);
      g_beam_plot->Draw("P SAME");
    }
    if (fbeam) {
      fbeam->SetLineColor(kViolet + 2);
      fbeam->SetLineWidth(2);
      fbeam->Draw("SAME");
    }
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "strip_alignment_check", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete g_beam_plot;
  }

  delete fbeam;
  delete g_cent;
  delete h2;

  return result;
}

void CalibrateBeam::CalibrateBeamOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &chans_template) {
  TString file_label = FileSet::FileLabel(spec);
  TString plot_subdir = "beam_calibration/" + file_label;
  std::cout << "Beam calibration: " << file_label << std::endl;

  BeamFit2D beam;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beam = FindBeamGateStp2VsStp1(spec, file_label, plot_subdir);
  }
  if (!beam.ok) {
    std::cerr << "  " << file_label << ": Strip2-vs-Strip1 beam gate failed"
              << std::endl;
    return;
  }

  BeamFit2D beam0vGrid;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beam0vGrid = FindBeamGateStp0VsGrid(spec, file_label, plot_subdir, beam);
  }
  if (!beam0vGrid.ok)
    std::cerr << "  " << file_label
              << ": Strip0-vs-Grid beam gate failed — continuing without it"
              << std::endl;

  std::vector<ChannelCal> chans = chans_template;
  std::vector<std::vector<Float_t>> samples;
  StripPairSamples pairs[18];
  CollectAnchorSamplesOneSubfile(spec, chans, beam, beam0vGrid, samples, pairs);
  std::vector<TF1 *> peak_fits;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, peak_fits, file_label, pairs);
    SaveBeamPeakChannelHistograms(chans, samples, peak_fits, plot_subdir);
  }
  for (Int_t i = 0; i < Int_t(peak_fits.size()); i++)
    delete peak_fits[i];
  peak_fits.clear();

  // Every channel that failed to calibrate is silently forced to gain 0 (reads
  // 0 a.u. and drops out of the strip total), so spell out which ones and why.
  // The long/short tag makes a long-side failure -- the dominant signal, which
  // should never starve -- easy to spot: grep "[uncalibrated long]".
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (IsCalibrated(c))
      continue;
    TString kind;
    if (c.side == 'C')
      kind = "cathode";
    else if (c.side == 'S')
      kind = "guard";
    else
      kind = (c.side == LongSide(c.strip)) ? "long" : "short";
    TString why;
    if (c.n_samples < kMinSamples)
      why =
          Form("too few beam samples (%lld < %lld)", c.n_samples, kMinSamples);
    else
      why = Form("bad exp anchor (fit_adc=%.1f)", c.fit_adc);
    std::cerr << "  [uncalibrated " << kind << "] " << c.name << " -> gain 0; "
              << why << " (beam n=" << c.n_samples << ")" << std::endl;
  }

  for (Int_t i = 0; i < Int_t(chans.size()); i++)
    if (IsCalibrated(chans[i]))
      std::cout << "  " << chans[i].name << " gain=" << Gain(chans[i])
                << " a.u./ADC  resolution="
                << Form("%.2f", ResolutionFWHMPercent(chans[i])) << "% FWHM"
                << std::endl;

  // Derive beam energy window from Strip0 (mu ± 3*sigma in a.u.).
  StripAlignmentResult align;
  DeriveBeamEnergyWindow(chans, align);

  // Write initial calibration tree: per-channel gains + beam window. The
  // alignment step needs this on disk so EnergyView can decode events in a.u.
  WriteCalibrationToEvents(spec, chans, &align);

  // Per-strip multiplicative alignment: decode events with the per-channel
  // gains just written, find each strip's beam-peak centroid, fit a pol3
  // reference trend, and derive factors = reference / centroid. Stored in
  // the calibration tree and applied by EnergyView as total *= factor.
  // gains just written, find each strip's beam-peak AND pileup-peak centroids,
  // and derive a linear correction (slope + intercept) that maps beam→1.0 and
  // pileup→2.0, flattening the Bragg curve. The polynomial trend is used only
  // Preserve beam energy window from DeriveBeamEnergyWindow across the
  // alignment call (which returns a fresh StripAlignmentResult).
  Double_t beam_e_min = align.beam_e_min;
  Double_t beam_e_max = align.beam_e_max;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    align = FindStripCentroidAlignment(spec, plot_subdir, file_label);
  }
  align.beam_e_min = beam_e_min;
  align.beam_e_max = beam_e_max;
  if (align.ok) {
    WriteCalibrationToEvents(spec, chans, &align);
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveDynamicRangeOverlay(spec, chans, plot_subdir, file_label);
  }
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveCalibSampleOverlay(chans, samples, plot_subdir, file_label);
  }

  std::cout << "  " << file_label << " calibration complete." << std::endl;
}

void WriteEresTomlRaw(const TString &out_subpath,
                      const Double_t eres_vals[35]) {
  toml::table eres_tbl;
  eres_tbl.insert("Cathode", eres_vals[0]);
  eres_tbl.insert("S0", eres_vals[1]);
  eres_tbl.insert("S17", eres_vals[2]);
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "L" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[3 + (s - 1)]);
  }
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "R" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[19 + (s - 1)]);
  }
  toml::table detector_tbl;
  detector_tbl.insert("eres", eres_tbl);
  toml::table root_tbl;
  root_tbl.insert("detector", detector_tbl);

  // The eres calibration TOML is a small, version-controlled input (a control
  // file), not bulk output: write it into the repo's control/ dir alongside the
  // other Calibration_Run*_eres.toml, regardless of where root_files point.
  TString out_dir = Paths::DatasetDir() + "/sim_control";
  gSystem->mkdir(out_dir, kTRUE);
  TString out_full = out_dir + "/" + out_subpath;
  std::ofstream f(out_full.Data());
  if (!f) {
    std::cerr << "Cannot write eres TOML: " << out_full << std::endl;
    return;
  }
  f << root_tbl << std::endl;
  std::cout << "  wrote eres TOML: " << out_full << std::endl;
}

void CalibrateBeam::AggregateEresTomlForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  const Int_t n_eres = 35;
  std::vector<std::vector<Double_t>> fwhm_per_chan(n_eres);

  std::vector<ChannelCal> tmpl = CalibrateBeam::BuildChannels();
  Int_t n_chans = Int_t(tmpl.size());

  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    // The calibration tree now lives inside each subfile's events file.
    TString cal_sub = FileSet::EventsName(specs[s]) + ".root";
    TFile *cf = IO::OpenForReading(cal_sub);
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t fit_adc[kMaxChannels] = {0};
    Float_t fit_sigma[kMaxChannels] = {0};
    Bool_t ok[kMaxChannels] = {0};
    t->SetBranchAddress("FitADC", fit_adc);
    t->SetBranchAddress("Ok", ok);
    t->SetBranchAddress("FitSigmaADC", fit_sigma);
    if (t->GetEntries() < 1) {
      cf->Close();
      delete cf;
      continue;
    }
    t->GetEntry(0);

    for (Int_t i = 0; i < n_chans && i < kMaxChannels; i++) {
      if (!ok[i])
        continue;
      Double_t sig_adc = fit_sigma[i];
      if (sig_adc <= 0 || fit_adc[i] <= 0)
        continue;
      // Relative resolution in % FWHM, straight from the raw-ADC peak fit —
      // independent of the normMUSIC gain/normalization by construction.
      const Double_t kFwhmPerSigma = 2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0));
      Double_t fwhm_pct =
          100.0 * kFwhmPerSigma * Double_t(sig_adc) / Double_t(fit_adc[i]);
      Int_t idx = ChannelToEresIndex(tmpl[i]);
      if (idx >= 0 && idx < n_eres)
        fwhm_per_chan[idx].push_back(fwhm_pct);
    }
    cf->Close();
    delete cf;
  }

  Double_t eres_vals[35];
  for (Int_t i = 0; i < n_eres; i++)
    eres_vals[i] = -1.0;
  for (Int_t i = 0; i < n_eres; i++) {
    std::vector<Double_t> &v = fwhm_per_chan[i];
    if (v.empty())
      continue;
    std::sort(v.begin(), v.end());
    Int_t m = Int_t(v.size());
    eres_vals[i] = (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
  }
  std::cout << "Run " << run << ": writing per-channel %FWHM medians"
            << std::endl;
  WriteEresTomlRaw(Form("Calibration_Run%d_eres.toml", run), eres_vals);
}

void CalibrateBeam::Run(const TString &file_label) {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = FileSet::BuildProcessedFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from FileSet::BuildProcessedFileSpecs()"
                << std::endl;
      return;
    }
  } else {
    FileSpec s = FileSet::ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  std::vector<ChannelCal> chans = CalibrateBeam::BuildChannels();

  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> runs;
  for (Int_t k = 0; k < n_specs; k++)
    runs.insert(specs[k].run);

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::cfg.MAX_FUSED_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "calibrate-beam: " << n_specs << " subfiles on " << n_workers
            << " workers" << std::endl;

  std::queue<Int_t> work;
  for (Int_t k = 0; k < n_specs; k++)
    work.push(k);
  std::mutex work_mutex;

  std::vector<std::thread> workers;
  for (Int_t w = 0; w < n_workers; w++) {
    workers.emplace_back([&]() {
      while (true) {
        Int_t k;
        {
          std::lock_guard<std::mutex> lk(work_mutex);
          if (work.empty())
            return;
          k = work.front();
          work.pop();
        }
        CalibrateBeamOneSubfile(specs[k], chans);
      }
    });
  }
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();

  for (std::set<Int_t>::const_iterator it = runs.begin(); it != runs.end();
       ++it) {
    std::vector<FileSpec> run_specs;
    for (Int_t k = 0; k < n_specs; k++)
      if (specs[k].run == *it)
        run_specs.push_back(specs[k]);
    AggregateEresTomlForRun(*it, run_specs);
  }
}
