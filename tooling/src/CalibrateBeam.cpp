#include "CalibrateBeam.hpp"

const Int_t kMaxChannels = 35;
// Axis-aligned gate half-widths for the Strip1-vs-Strip0 beam ellipse: tight
// horizontally (Strip0), loose vertically (Strip1).
const Double_t kEllipseNSigmaX = 1.0;
const Double_t kEllipseNSigmaY = 3.0;
const Long64_t kMinSamples = 200;
const Long64_t kSampleCap = 20000;

TString CalibrateBeam::DefaultSimBeamPath(const TString &project_root) {
  TString p = project_root + "/sim_root_files/" + Constants::SIM_BEAM_FILE;
  return TString(gSystem->ExpandPathName(p.Data()));
}

std::vector<ChannelCal> BuildChannels() {
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

Char_t LongSide(Int_t strip) { return (strip % 2 == 0) ? 'R' : 'L'; }

Bool_t IsBeamdEChannel(const ChannelCal &c) {
  if (c.side == 'S')
    return kTRUE;
  if (c.side != 'L' && c.side != 'R')
    return kFALSE;
  if (c.strip < 1 || c.strip > 16)
    return kFALSE;
  return c.side == LongSide(c.strip);
}

// Standard ("normMUSIC") calibration: each beam-dE channel's beam-peak ADC is
// scaled to the common reference NORM_MUSIC_MEV by a single gain through the
// origin. Calibrated dE is normalized to a common beam value, not absolute --
// which is what the PID / summed-dE analysis needs. sim_mu_mev>0 selects the
// channels that have a sim beam anchor (the long sides); the sim energy itself
// is no longer the calibration target.
Bool_t IsCalibrated(const ChannelCal &c) {
  return c.sim_mu_mev > 0 && c.fit_adc > 0;
}

Double_t Gain(const ChannelCal &c) {
  return Constants::NORM_MUSIC_MEV / c.fit_adc;
}

// Energy resolution at the beam peak as percent FWHM. Scale-invariant (the gain
// cancels), so it is computed straight from the ADC-domain fit: FWHM = 2.3548 *
// sigma, %FWHM = 100 * FWHM / centroid.
Double_t ResolutionFWHMPercent(const ChannelCal &c) {
  if (c.fit_adc <= 0)
    return 0.0;
  const Double_t kFwhmPerSigma = 2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0));
  return 100.0 * kFwhmPerSigma * c.fit_sigma_adc / c.fit_adc;
}

inline Double_t ApplyCal(const ChannelCal &c, Double_t adc) {
  return Gain(c) * adc;
}

// Beam selection: the beam is the dominant population in Strip1 dE vs Strip0
// dE, so we take the largest peak there and gate on a 2D ellipse built from
// local moments (no L-strip-sum gate, no TF2).
BeamFit2D FindBeamGateStp1VsStp0(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir,
                                 Bool_t save_plot = kTRUE) {
  BeamFit2D out;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      delete sf;
    return out;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return out;
  }
  // Raw ADC, read before any calibration exists (this fit produces it), so read
  // the branches directly. Strip0 = Left_0_17_dE[0]; strip1 total = L1 + R1.
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);

  TH2F *h = new TH2F(Form("h2_stp1_vs_stp0_%s", run_label.Data()),
                     ";Strip0 #DeltaE [ADC];Strip1 #DeltaE [ADC]", 256, 0.0,
                     16384.0, 256, 0.0, 16384.0);
  h->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    Int_t stp0 = Int_t(left_0_17_adc[0]);
    Int_t stp1 = Int_t(left_0_17_adc[1]) + Int_t(rightdE_adc[1]);
    if (stp0 > 0 && stp1 > 0)
      h->Fill(Double_t(stp0), Double_t(stp1));
  }
  sf->Close();
  delete sf;

  if (h->GetEntries() < 100) {
    std::cerr << "  " << run_label
              << ": too few events for Strip1-vs-Strip0 beam gate" << std::endl;
    delete h;
    return out;
  }

  const Double_t kSeedFrac = 0.30;
  const Int_t kSeedHalfBins = 40;
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  " << run_label << ": no bins above beam seed threshold"
              << std::endl;
    delete h;
    return out;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;
  std::cout << "  beam gate (Strip1 vs Strip0): mu=(" << out.mu_x << ","
            << out.mu_y << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;

  if (save_plot) {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, cv);
    TEllipse *e =
        new TEllipse(out.mu_x, out.mu_y, kEllipseNSigmaX * out.sigma_x,
                     kEllipseNSigmaY * out.sigma_y);
    e->SetFillStyle(0);
    e->SetLineColor(kRed + 1);
    e->SetLineWidth(2);
    e->Draw();
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "beam_gate_stp1_vs_stp0", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
  }
  delete h;
  return out;
}

// Upstream normEsegment peak extraction: locate the beam peak with
// TSpectrum::Search(hist, sigma=2, "", threshold=0.9), then fit a plain "gaus"
// in a window [peak*0.9, peak*1.1] around the found position. Returns the
// fitted centroid (peak_adc, = upstream gCalib->GetParameter(1)) and sigma
// (sigma_adc, = GetParameter(2)). Returns kFALSE on failure; the caller falls
// back to median/IQR. The samples are already beam-gated by the 2D ellipse, so
// the threshold=0.9 search returns essentially the single beam peak.
Bool_t FitTSpectrumGaussianPeak(const std::vector<Float_t> &v,
                                const TString &fname, Double_t &peak_adc,
                                Double_t &sigma_adc, TF1 *&fit_out) {
  fit_out = nullptr;
  if (v.size() < kMinSamples)
    return kFALSE;
  Float_t lo = v[0], hi = v[0];
  for (Int_t j = 1; j < Int_t(v.size()); j++) {
    if (v[j] < lo)
      lo = v[j];
    if (v[j] > hi)
      hi = v[j];
  }
  Double_t pad = 0.05 * (Double_t(hi) - Double_t(lo));
  if (pad < 1.0)
    pad = 1.0;
  const Int_t nbins = 100;
  Double_t xlo = Double_t(lo) - pad;
  Double_t xhi = Double_t(hi) + pad;
  TH1F h(fname + "_h", "", nbins, xlo, xhi);
  h.SetDirectory(nullptr);
  for (Int_t j = 0; j < Int_t(v.size()); j++)
    h.Fill(Double_t(v[j]));

  // Upstream peak finder: sigma=2 bins, threshold=0.9 (only peaks within 90% of
  // the tallest survive). "nodraw" suppresses the marker polymarker so nothing
  // leaks onto the batch canvases.
  TSpectrum spec;
  Int_t npeaks = spec.Search(&h, 2, "nodraw", 0.9);
  if (npeaks < 1)
    return kFALSE;
  Double_t *xpeaks = spec.GetPositionX();
  // Take the tallest of the returned peaks as the beam peak. (Upstream takes
  // index 0; with threshold=0.9 there is usually only one, but selecting by
  // bin content is robust to TSpectrum's position ordering.)
  Int_t best = 0;
  Double_t best_val = h.GetBinContent(h.FindBin(xpeaks[0]));
  for (Int_t p = 1; p < npeaks; p++) {
    Double_t val = h.GetBinContent(h.FindBin(xpeaks[p]));
    if (val > best_val) {
      best_val = val;
      best = p;
    }
  }
  Double_t peak_pos = xpeaks[best];
  if (!(peak_pos > 0))
    return kFALSE;

  // Gaussian fit in the upstream [peak*0.9, peak*1.1] window around the peak.
  Double_t fit_lo = peak_pos * 0.9;
  Double_t fit_hi = peak_pos * 1.1;
  if (fit_lo < xlo)
    fit_lo = xlo;
  if (fit_hi > xhi)
    fit_hi = xhi;
  Double_t bw = h.GetBinWidth(1);
  Double_t sigma_seed = (fit_hi - fit_lo) / 4.0;
  if (sigma_seed < bw)
    sigma_seed = bw;

  TF1 *f = new TF1(fname, "gaus", fit_lo, fit_hi);
  f->SetNpx(1000);
  f->SetParameters(best_val, peak_pos, sigma_seed);
  f->SetParLimits(1, fit_lo, fit_hi);
  f->SetParLimits(2, bw, fit_hi - fit_lo);
  TFitResultPtr r = h.Fit(f, "QSRNL");
  if (!r.Get() || !r->IsValid()) {
    delete f;
    return kFALSE;
  }
  peak_adc = f->GetParameter(1);
  sigma_adc = std::fabs(f->GetParameter(2));
  if (!(peak_adc > 0) || !(sigma_adc > 0)) {
    delete f;
    return kFALSE;
  }
  fit_out = f;
  return kTRUE;
}

inline Double_t Median(std::vector<Float_t> &v) {
  if (v.empty())
    return 0.0;
  size_t n = v.size();
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  Double_t med = Double_t(v[n / 2]);
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.end());
    med = 0.5 * (med + Double_t(v[n / 2 - 1]));
  }
  return med;
}

// IQR = Q3 - Q1; divide by 1.349 outside this helper for the Gaussian-sigma
// approximation when needed.
inline Double_t InterquartileRange(std::vector<Float_t> &v) {
  if (v.size() < 4)
    return 0.0;
  size_t n = v.size();
  size_t i1 = n / 4;
  size_t i3 = (3 * n) / 4;
  std::nth_element(v.begin(), v.begin() + i1, v.end());
  Double_t q1 = Double_t(v[i1]);
  std::nth_element(v.begin(), v.begin() + i3, v.end());
  Double_t q3 = Double_t(v[i3]);
  return q3 - q1;
}

// Fit a Gaussian to the bucket and return (mu, sigma). Falls back to
// (median, IQR/1.349) on fit failure. Sim per-channel deposits are
// well-approximated by a Gaussian, so a direct fit gives a cleaner
// (mu, sigma) than median/IQR estimators.
Bool_t FitGaussianMuSigma(const std::vector<Float_t> &v, const TString &fname,
                          Double_t &mu, Double_t &sigma) {
  if (v.size() < 50)
    return kFALSE;
  Float_t lo = v[0], hi = v[0];
  for (Int_t j = 1; j < Int_t(v.size()); j++) {
    if (v[j] < lo)
      lo = v[j];
    if (v[j] > hi)
      hi = v[j];
  }
  Double_t pad = 0.05 * (Double_t(hi) - Double_t(lo));
  if (pad < 1e-6)
    pad = 1e-6;
  const Int_t nbins = 75;
  TH1F h(fname + "_h", "", nbins, Double_t(lo) - pad, Double_t(hi) + pad);
  h.SetDirectory(nullptr);
  for (Int_t j = 0; j < Int_t(v.size()); j++)
    h.Fill(Double_t(v[j]));
  TF1 fg(fname, "gaus", Double_t(lo) - pad, Double_t(hi) + pad);
  Int_t pb = h.GetMaximumBin();
  fg.SetParameters(h.GetBinContent(pb), h.GetBinCenter(pb), h.GetRMS());
  TFitResultPtr r = h.Fit(&fg, "QSRN");
  if (!r.Get() || !r->IsValid())
    return kFALSE;
  mu = fg.GetParameter(1);
  sigma = std::fabs(fg.GetParameter(2));
  return mu > 0 && sigma > 0;
}

// Load sim per-channel anchors from the events_MeV tree by Gaussian-fitting
// the nonzero per-strip MeV deposit. Falls back to median/IQR on fit failure.
Bool_t LoadSimMeans(const TString &sim_path, std::vector<ChannelCal> &chans) {
  TFile *f = TFile::Open(sim_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open sim beam file: " << sim_path << std::endl;
    if (f)
      delete f;
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "Sim file has no 'events_MeV' tree" << std::endl;
    f->Close();
    delete f;
    return kFALSE;
  }
  Float_t leftdE[18], rightdE[18], totaldE[18];
  Float_t cathode = 0;
  t->SetBranchAddress("LeftdE", leftdE);
  t->SetBranchAddress("RightdE", rightdE);
  t->SetBranchAddress("TotaldE", totaldE);
  t->SetBranchAddress("Cathode", &cathode);

  std::vector<std::vector<Float_t>> buckets(chans.size());
  Long64_t n = t->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    t->GetEntry(j);
    for (Int_t i = 0; i < Int_t(chans.size()); i++) {
      const ChannelCal &c = chans[i];
      Float_t v = 0.0f;
      if (c.side == 'S')
        v = totaldE[c.strip];
      else if (c.side == 'L')
        v = leftdE[c.strip];
      else if (c.side == 'R')
        v = rightdE[c.strip];
      else if (c.side == 'C')
        v = cathode;
      if (v > 0.0f)
        buckets[i].push_back(v);
    }
  }
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    ChannelCal &c = chans[i];
    if (buckets[i].size() < 50) {
      std::cerr << "  sim " << c.name << ": only " << buckets[i].size()
                << " nonzero entries; no anchor" << std::endl;
      continue;
    }
    Double_t mu = 0, sigma = 0;
    TString fname = Form("f_sim_gaus_%s", c.name.Data());
    if (FitGaussianMuSigma(buckets[i], fname, mu, sigma)) {
      c.sim_mu_mev = mu;
      c.sim_sigma_mev = sigma;
    } else {
      std::cerr << "  sim " << c.name
                << ": gaussian fit failed; falling back to median/IQR"
                << std::endl;
      c.sim_mu_mev = Median(buckets[i]);
      c.sim_sigma_mev = InterquartileRange(buckets[i]) / 1.349;
    }
    std::cout << "  sim " << c.name << ": mu=" << c.sim_mu_mev
              << " MeV, sigma=" << c.sim_sigma_mev
              << " (n=" << buckets[i].size() << ")" << std::endl;
  }
  f->Close();
  delete f;
  return kTRUE;
}

std::vector<ChannelCal> CalibrateBeam::LoadSimChans(const TString &sim_path) {
  std::vector<ChannelCal> chans = BuildChannels();
  if (!LoadSimMeans(sim_path, chans))
    return std::vector<ChannelCal>();
  return chans;
}

void CollectAnchorSamplesOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &chans,
    const BeamFit2D &beam, std::vector<std::vector<Float_t>> &samples) {
  Int_t n_chans = Int_t(chans.size());
  samples.assign(n_chans, std::vector<Float_t>());
  if (!beam.ok)
    return;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return;
  }
  // Raw ADC, pre-calibration. Guard strips (S) and left ends (L) live in
  // Left_0_17_dE; right ends in RightdE. Strip totals are L+R; the gate uses
  // Strip0 (=Left_0_17_dE[0]) and the strip1 total (L1+R1).
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t cathode_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    Double_t x = Double_t(left_0_17_adc[0]);
    Double_t y = Double_t(left_0_17_adc[1]) + Double_t(rightdE_adc[1]);
    if (x <= 0 || y <= 0)
      continue;
    if (!BeamFitUtils::InEllipseXY(beam, x, y, kEllipseNSigmaX,
                                   kEllipseNSigmaY))
      continue;
    for (Int_t i = 0; i < n_chans; i++) {
      if (Long64_t(samples[i].size()) >= kSampleCap)
        continue;
      const ChannelCal &c = chans[i];
      Int_t v = 0;
      if (c.side == 'S' || c.side == 'L')
        v = Int_t(left_0_17_adc[c.strip]);
      else if (c.side == 'R')
        v = Int_t(rightdE_adc[c.strip]);
      else if (c.side == 'C')
        v = Int_t(cathode_adc);
      if (v > 0)
        samples[i].push_back(Float_t(v));
    }
  }
  sf->Close();
  delete sf;
}

// Cathode uses median + IQR/1.349 (asymmetric tail not as clean and the user
// prefers to keep cathode on the existing approach). All other channels
// (S guard strips + L/R long anodes) use the upstream normEsegment recipe:
// TSpectrum locates the beam peak, a Gaussian is fit in a [peak*0.9, peak*1.1]
// window, and the fitted centroid is the gain anchor while the sigma is the
// detector resolution. Fall back to mean + IQR/1.349 on fit failure.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     std::vector<TF1 *> &fits_out, const TString &run_label) {
  Int_t n_chans = Int_t(chans.size());
  fits_out.assign(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    ChannelCal &c = chans[i];
    std::vector<Float_t> &v = samples[i];
    c.n_samples = Long64_t(v.size());
    if (Long64_t(v.size()) < kMinSamples) {
      c.fit_adc = 0;
      c.fit_sigma_adc = 0;
      continue;
    }
    // Sample mean kept only as the fallback anchor if the TSpectrum+Gaussian
    // fit fails.
    Double_t mean_adc = 0.0;
    for (Int_t j = 0; j < Int_t(v.size()); j++)
      mean_adc += Double_t(v[j]);
    mean_adc /= Double_t(v.size());

    if (c.side == 'C') {
      // Cathode: median + IQR (asymmetric tail, no clean peak).
      c.fit_adc = Median(v);
      c.fit_sigma_adc = InterquartileRange(v) / 1.349;
    } else {
      // Upstream recipe: anchor = Gaussian centroid around the TSpectrum peak;
      // sigma = Gaussian width. Mean + IQR fallback on fit failure.
      TString fname =
          Form("f_tspec_gaus_%s_%s", c.name.Data(), run_label.Data());
      Double_t peak = 0, sig = 0;
      TF1 *fit = nullptr;
      if (FitTSpectrumGaussianPeak(v, fname, peak, sig, fit)) {
        c.fit_adc = peak;
        c.fit_sigma_adc = sig;
        fits_out[i] = fit;
      } else {
        c.fit_adc = mean_adc;
        c.fit_sigma_adc = InterquartileRange(v) / 1.349;
      }
    }
    std::cout << "  " << c.name << " anchor[ADC]=" << c.fit_adc
              << " sig=" << c.fit_sigma_adc << " (n=" << c.n_samples << ")"
              << std::endl;
  }
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
  TString out_dir = Paths::DatasetDir() + "/control";
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

// Writes a one-row `calibration` tree into the open file `dst` (typically a
// per-subfile .cal.root). Layout matches AggregateEresTomlForRun's reader.
void WriteCalibrationTree(TFile *dst, const std::vector<ChannelCal> &chans) {
  dst->cd();
  if (TObject *old = dst->Get("calibration"))
    old->Delete();
  TTree *cal = new TTree("calibration", "Per-channel normMUSIC calibration");
  Float_t gain[kMaxChannels] = {0};
  Float_t sim_mu[kMaxChannels] = {0}, sim_sigma[kMaxChannels] = {0};
  Float_t fit_adc[kMaxChannels] = {0}, fit_sigma[kMaxChannels] = {0};
  Long64_t fit_n[kMaxChannels] = {0};
  Bool_t ok[kMaxChannels] = {0};
  // Per-strip gains laid out to match the events tree exactly: GainLeft[s]
  // multiplies Left_0_17_dE[s] (s=0/17 are the single-ended guards, s=1..16 the
  // left ends), GainRight[s] multiplies RightdE[s] (0 at the guards). This is
  // what EnergyView reads to calibrate on the fly -- no per-event MeV is
  // stored.
  Float_t gain_left[18] = {0}, gain_right[18] = {0};
  Float_t gain_cathode = 0.0f;
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    ok[k] = IsCalibrated(c);
    gain[k] = ok[k] ? Float_t(Gain(c)) : 0.0f;
    sim_mu[k] = Float_t(c.sim_mu_mev);
    sim_sigma[k] = Float_t(c.sim_sigma_mev);
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
  cal->Branch("SimMu_MeV", sim_mu, Form("SimMu_MeV[%d]/F", kMaxChannels));
  cal->Branch("SimSigma_MeV", sim_sigma,
              Form("SimSigma_MeV[%d]/F", kMaxChannels));
  cal->Branch("FitADC", fit_adc, Form("FitADC[%d]/F", kMaxChannels));
  cal->Branch("FitSigmaADC", fit_sigma,
              Form("FitSigmaADC[%d]/F", kMaxChannels));
  cal->Branch("FitN", fit_n, Form("FitN[%d]/L", kMaxChannels));
  cal->Branch("GainLeft", gain_left, "GainLeft[18]/F");
  cal->Branch("GainRight", gain_right, "GainRight[18]/F");
  cal->Branch("GainCathode", &gain_cathode, "GainCathode/F");
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
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
    Float_t lo = v[0], hi = v[0];
    for (Int_t j = 1; j < Int_t(v.size()); j++) {
      if (v[j] < lo)
        lo = v[j];
      if (v[j] > hi)
        hi = v[j];
    }
    Double_t pad = 0.05 * (Double_t(hi) - Double_t(lo));
    if (pad < 1.0)
      pad = 1.0;
    const Int_t nbins = 75;
    TH1F *h = new TH1F(Form("h_beam_peak_%s", c.name.Data()),
                       Form(";%s #DeltaE [ADC];Counts", c.name.Data()), nbins,
                       Double_t(lo) - pad, Double_t(hi) + pad);
    h->SetDirectory(nullptr);
    for (Int_t j = 0; j < Int_t(v.size()); j++)
      h->Fill(Double_t(v[j]));
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h, kBlack);
    TF1 *fit = fits[i];
    if (fit) {
      fit->SetLineColor(kRed + 1);
      fit->SetLineWidth(2);
      fit->Draw("L SAME");
    }
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, Form("beam_peak_%s", c.name.Data()), subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete h;
  }
}

// Writes the per-channel gain table (tree "calibration") into the subfile's own
// events file. No per-event calibrated tree is produced: downstream readers
// recover MeV on the fly via gain x raw ADC (EnergyView), so the raw events
// tree plus this one-row gain table fully determine every calibrated value.
void WriteCalibrationToEvents(const FileSpec &spec,
                              const std::vector<ChannelCal> &chans) {
  TString events_subpath = FileSet::EventsName(spec) + ".root";
  TFile *f = IO::OpenForWriting(events_subpath, "UPDATE");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << " to write calibration"
              << std::endl;
    if (f)
      delete f;
    return;
  }
  WriteCalibrationTree(f, chans);
  std::cout << "  wrote calibration into " << events_subpath << std::endl;
  f->Close();
  delete f;
}

// Per-channel calibrated overlay for one subfile, via AttachCalSidecar (the
// same path downstream macros use). The sidecar must already be on disk.
void SaveDynamicRangeOverlay(const FileSpec &spec,
                             const std::vector<ChannelCal> &chans,
                             const TString &plot_subdir,
                             const TString &file_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = Constants::STRIP_E_MIN_MEV;
  const Double_t emax = Constants::STRIP_E_MAX_MEV;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    if (!IsBeamdEChannel(chans[i]))
      continue;
    TString hname =
        Form("h_dynrange_%s_%s", file_label.Data(), chans[i].name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [MeV];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
  }
  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    for (Int_t i = 0; i < n_chans; i++)
      delete h[i];
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    for (Int_t i = 0; i < n_chans; i++)
      delete h[i];
    return;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_mev) {
    sf->Close();
    delete sf;
    for (Int_t i = 0; i < n_chans; i++)
      delete h[i];
    return;
  }
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    for (Int_t i = 0; i < n_chans; i++) {
      if (!h[i])
        continue;
      const ChannelCal &c = chans[i];
      Double_t v = 0.0;
      if (c.side == 'S')
        v = ev.total[c.strip];
      else if (c.side == 'L')
        v = ev.left[c.strip];
      else
        v = ev.right[c.strip];
      if (v > 0)
        h[i]->Fill(v);
    }
  }
  sf->Close();
  delete sf;
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Double_t m = h[i]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Int_t color = colors[i % Int_t(colors.size())];
    h[i]->SetLineColor(color);
    h[i]->SetLineWidth(2);
    h[i]->SetMaximum(1.15 * y_top);
    h[i]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    leg->AddEntry(h[i], chans[i].name.Data(), "l");
  }
  leg->Draw();
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}

// Overlay (one color per channel, log-y) of ONLY the events used for
// calibration: the beam anchor samples, converted to MeV via each channel's
// gain. Same axes/style as SaveDynamicRangeOverlay but restricted to
// calibration events rather than the full spectrum.
void CalibrateBeam::SaveCalibSampleOverlay(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<Float_t>> &samples,
    const TString &plot_subdir, const TString &file_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = Constants::STRIP_E_MIN_MEV;
  const Double_t emax = Constants::STRIP_E_MAX_MEV;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    const ChannelCal &c = chans[i];
    if (!IsBeamdEChannel(c) || !IsCalibrated(c))
      continue;
    TString hname =
        Form("h_calibrange_%s_%s", file_label.Data(), c.name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [MeV];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
    const std::vector<Float_t> &v = samples[i];
    for (Int_t j = 0; j < Int_t(v.size()); j++) {
      Double_t mev = ApplyCal(c, Double_t(v[j]));
      if (mev > 0)
        h[i]->Fill(mev);
    }
  }

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Double_t m = h[i]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Int_t color = colors[i % Int_t(colors.size())];
    h[i]->SetLineColor(color);
    h[i]->SetLineWidth(2);
    h[i]->SetMaximum(1.15 * y_top);
    h[i]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    leg->AddEntry(h[i], chans[i].name.Data(), "l");
  }
  leg->Draw();
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_calib_events", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}

void CalibrateBeam::CalibrateBeamOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &sim_chans) {
  TString file_label = FileSet::FileLabel(spec);
  TString plot_subdir = "beam_calibration/" + file_label;
  std::cout << "Beam calibration: " << file_label << std::endl;

  if (sim_chans.empty()) {
    std::cerr << "  " << file_label
              << ": sim anchors empty; skipping calibration" << std::endl;
    return;
  }

  BeamFit2D beam;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beam = FindBeamGateStp1VsStp0(spec, file_label, plot_subdir);
  }
  if (!beam.ok) {
    std::cerr << "  " << file_label << ": Strip1-vs-Strip0 beam gate failed"
              << std::endl;
    return;
  }

  std::vector<ChannelCal> chans = sim_chans;
  std::vector<std::vector<Float_t>> samples;
  CollectAnchorSamplesOneSubfile(spec, chans, beam, samples);
  std::vector<TF1 *> peak_fits;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, peak_fits, file_label);
    SaveBeamPeakChannelHistograms(chans, samples, peak_fits, plot_subdir);
  }
  for (Int_t i = 0; i < Int_t(peak_fits.size()); i++)
    delete peak_fits[i];
  peak_fits.clear();

  for (Int_t i = 0; i < Int_t(chans.size()); i++)
    if (IsCalibrated(chans[i]))
      std::cout << "  " << chans[i].name << " gain=" << Gain(chans[i])
                << " MeV/ADC  resolution="
                << Form("%.2f", ResolutionFWHMPercent(chans[i])) << "% FWHM"
                << std::endl;

  WriteCalibrationToEvents(spec, chans);

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

void CalibrateBeam::AggregateEresTomlForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  const Int_t n_eres = 35;
  std::vector<std::vector<Double_t>> sigma_per_chan(n_eres);

  std::vector<ChannelCal> tmpl = BuildChannels();
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
    Float_t gain[kMaxChannels] = {0};
    Float_t fit_sigma[kMaxChannels] = {0};
    Bool_t ok[kMaxChannels] = {0};
    t->SetBranchAddress("Gain", gain);
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
      if (sig_adc <= 0 || gain[i] <= 0)
        continue;
      Double_t sigma_mev = Double_t(sig_adc) * Double_t(gain[i]);
      Int_t idx = ChannelToEresIndex(tmpl[i]);
      if (idx >= 0 && idx < n_eres)
        sigma_per_chan[idx].push_back(sigma_mev);
    }
    cf->Close();
    delete cf;
  }

  Double_t eres_vals[35];
  for (Int_t i = 0; i < n_eres; i++)
    eres_vals[i] = -1.0;
  for (Int_t i = 0; i < n_eres; i++) {
    std::vector<Double_t> &v = sigma_per_chan[i];
    if (v.empty())
      continue;
    std::sort(v.begin(), v.end());
    size_t m = v.size();
    eres_vals[i] = (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
  }
  std::cout << "Run " << run << ": writing per-channel sigma (MeV) medians"
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
    specs = FileSet::BuildFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from FileSet::BuildFileSpecs()" << std::endl;
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

  TString sim_beam_path = DefaultSimBeamPath(project_root);
  std::vector<ChannelCal> sim_chans = LoadSimChans(sim_beam_path);
  if (sim_chans.empty()) {
    std::cerr << "Sim anchors unavailable; aborting" << std::endl;
    return;
  }

  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> runs;
  for (Int_t k = 0; k < n_specs; k++)
    runs.insert(specs[k].run);

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::MAX_FUSED_WORKERS);
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
        CalibrateBeamOneSubfile(specs[k], sim_chans);
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
