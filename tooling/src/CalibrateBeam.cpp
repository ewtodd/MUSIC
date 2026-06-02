#include "CalibrateBeam.hpp"

const Int_t kMaxChannels = 35;
const Double_t kEllipseNSigma = 2.5;
// Axis-aligned gate half-widths for the Strip1-vs-Strip0 beam/pileup ellipses:
// tight horizontally (Strip0), loose vertically (Strip1).
const Double_t kEllipseNSigmaX = 1.0;
const Double_t kEllipseNSigmaY = 3.0;
const Long64_t kMinSamplesPerAnchor = 200;
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

TH2F *ExtractH2(TFile *f, const TString &name) {
  TObject *obj = f->Get(name);
  if (!obj)
    return nullptr;
  if (obj->InheritsFrom(TH2F::Class()))
    return static_cast<TH2F *>(obj);
  if (!obj->InheritsFrom(TCanvas::Class()))
    return nullptr;
  TCanvas *canv = static_cast<TCanvas *>(obj);
  TList *prims = canv->GetListOfPrimitives();
  for (Int_t i = 0; i < prims->GetSize(); i++) {
    TObject *p = prims->At(i);
    if (p && p->InheritsFrom(TH2F::Class()))
      return static_cast<TH2F *>(p);
  }
  return nullptr;
}

// Beam selection for the standard single-scale calibration: the beam is the
// dominant population in Strip1 dE vs Strip0 dE, so we take the largest peak
// there and gate on a 2D ellipse built from local moments (no L-strip-sum gate,
// no TF2). Only the beam peak is needed -- no pileup (k=2) anchor.
std::vector<BeamFit2D> FindBeamGateStp1VsStp0(const FileSpec &spec,
                                              const TString &run_label,
                                              const TString &plot_subdir,
                                              Bool_t save_plot = kTRUE) {
  std::vector<BeamFit2D> out(kNPeaks);

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
  Int_t totaldE_adc[18];
  tree->SetBranchAddress("TotaldE", totaldE_adc);

  TH2F *h = new TH2F(Form("h2_stp1_vs_stp0_%s", run_label.Data()),
                     ";Strip0 #DeltaE [ADC];Strip1 #DeltaE [ADC]", 256, 0.0,
                     16384.0, 256, 0.0, 16384.0);
  h->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    if (totaldE_adc[0] > 0 && totaldE_adc[1] > 0)
      h->Fill(Double_t(totaldE_adc[0]), Double_t(totaldE_adc[1]));
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
  out[0].amp = peak_val;
  out[0].mu_x = m.mu_x;
  out[0].mu_y = m.mu_y;
  out[0].sigma_x = m.sigma_x;
  out[0].sigma_y = m.sigma_y;
  out[0].rho = m.rho;
  out[0].ok = kTRUE;
  std::cout << "  beam gate (Strip1 vs Strip0): mu=(" << out[0].mu_x << ","
            << out[0].mu_y << ") sigma=(" << out[0].sigma_x << ","
            << out[0].sigma_y << ") rho=" << out[0].rho << std::endl;

  // Standard calibration uses only the beam (k=1) peak -> single gain to
  // NORM_MUSIC_MEV, so no pileup (k=2) anchor is gated or collected.

  if (save_plot) {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, cv);
    for (Int_t k = 0; k < kNPeaks; k++) {
      if (!out[k].ok)
        continue;
      TEllipse *e = new TEllipse(out[k].mu_x, out[k].mu_y,
                                 kEllipseNSigmaX * out[k].sigma_x,
                                 kEllipseNSigmaY * out[k].sigma_y);
      e->SetFillStyle(0);
      e->SetLineColor(kRed + 1);
      e->SetLineWidth(2);
      e->Draw();
    }
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "beam_gate_stp1_vs_stp0", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
  }
  delete h;
  return out;
}

// Landau (intrinsic energy-loss fluctuation) convolved with Gaussian (detector
// resolution). par[0]=Landau width, par[1]=Landau MPV (after the standard
// -0.22278298*width shift, so par[1] is the true MPV of the unsmeared Landau),
// par[2]=area, par[3]=Gaussian sigma. Standard CERN langaus midpoint-rule
// implementation.
Double_t LangausFunction(Double_t *x, Double_t *par) {
  const Double_t kInvSq2Pi = 0.3989422804014;
  const Double_t kMpShift = -0.22278298;
  const Int_t kNSteps = 100;
  const Double_t kRangeSigmas = 5.0;
  Double_t mpc = par[1] - kMpShift * par[0];
  Double_t xlow = x[0] - kRangeSigmas * par[3];
  Double_t xupp = x[0] + kRangeSigmas * par[3];
  Double_t step = (xupp - xlow) / Double_t(kNSteps);
  Double_t sum = 0.0;
  for (Int_t i = 1; i <= kNSteps / 2; i++) {
    Double_t xx = xlow + (Double_t(i) - 0.5) * step;
    Double_t fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (Double_t(i) - 0.5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * kInvSq2Pi / par[3];
}

// Fit a langaus to the sample distribution and return the convolved peak
// position (peak_adc) and the Gaussian sigma (sigma_gauss_adc). Returns
// kFALSE on failure; the caller should fall back to median/IQR.
Bool_t FitLangausPeak(const std::vector<Float_t> &v, const TString &fname,
                      Double_t &peak_adc, Double_t &sigma_gauss_adc,
                      TF1 *&fit_out, Double_t mpv_seed_override = -1.0) {
  fit_out = nullptr;
  if (v.size() < kMinSamplesPerAnchor)
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

  // k=2 fits are seeded at 2x the k=1 MPV (passed in) so the pileup fit can't
  // drift onto the much taller beam peak that contaminates the k=2 gate; k=1
  // self-seeds from the sample histogram maximum.
  Int_t peak_bin;
  if (mpv_seed_override > 0.0) {
    peak_bin = h.FindBin(mpv_seed_override);
    if (peak_bin < 1)
      peak_bin = 1;
    if (peak_bin > nbins)
      peak_bin = nbins;
  } else {
    peak_bin = h.GetMaximumBin();
  }
  Double_t peak_val = h.GetBinContent(peak_bin);
  Double_t mpv_seed = h.GetBinCenter(peak_bin);
  Int_t lo_half = peak_bin;
  while (lo_half > 1 && h.GetBinContent(lo_half) > 0.5 * peak_val)
    lo_half--;
  Int_t hi_half = peak_bin;
  while (hi_half < nbins && h.GetBinContent(hi_half) > 0.5 * peak_val)
    hi_half++;
  Double_t fwhm_lhalf = mpv_seed - h.GetBinCenter(lo_half);
  Double_t fwhm_rhalf = h.GetBinCenter(hi_half) - mpv_seed;
  Double_t fwhm = fwhm_lhalf + fwhm_rhalf;
  if (fwhm <= 0)
    fwhm = (xhi - xlo) / 10.0;
  if (fwhm_lhalf <= 0)
    fwhm_lhalf = 0.4 * fwhm;
  if (fwhm_rhalf <= 0)
    fwhm_rhalf = 0.6 * fwhm;
  // Asymmetric seeding: left side is Gaussian-dominated (HWHM_l ~ 1.177*sigma),
  // right-side excess over left is Landau-dominated. This breaks the seed
  // degeneracy where both widths get pulled from the same FWHM number.
  Double_t gauss_sigma_seed = fwhm_lhalf / 1.177;
  Double_t landau_excess = fwhm_rhalf - fwhm_lhalf;
  Double_t landau_width_seed =
      (landau_excess > 0) ? landau_excess / 2.39 : fwhm / 8.0;
  Double_t area_seed = Double_t(v.size()) * h.GetBinWidth(1);
  Double_t fit_lo = mpv_seed - 1.5 * fwhm;
  Double_t fit_hi = mpv_seed + 3.0 * fwhm;
  if (fit_lo < xlo)
    fit_lo = xlo;
  if (fit_hi > xhi)
    fit_hi = xhi;

  TF1 *f = new TF1(fname, LangausFunction, fit_lo, fit_hi, 4);
  f->SetNpx(1000);
  f->SetParNames("LandauWidth", "MPV", "Area", "GaussSigma");
  f->SetParameters(landau_width_seed, mpv_seed, area_seed, gauss_sigma_seed);
  // Floor both width params at the bin width: a Landau-or-Gauss narrower than
  // one bin lets the convolution spike across bin centers and chase noise.
  Double_t bw = h.GetBinWidth(1);
  f->SetParLimits(0, bw, 5.0 * fwhm);
  f->SetParLimits(1, mpv_seed - fwhm, mpv_seed + fwhm);
  f->SetParLimits(2, 0.0, 1e6 * area_seed);
  f->SetParLimits(3, bw, 5.0 * fwhm);
  TFitResultPtr r = h.Fit(f, "QSRNL");
  if (!r.Get() || !r->IsValid()) {
    delete f;
    return kFALSE;
  }
  peak_adc = f->GetMaximumX(fit_lo, fit_hi);
  sigma_gauss_adc = f->GetParameter(3);
  if (!(peak_adc > 0) || !(sigma_gauss_adc > 0)) {
    delete f;
    return kFALSE;
  }
  fit_out = f;
  return kTRUE;
}

// Fit a plain Gaussian to the sample distribution around its peak and return
// the mean (peak_adc) and sigma. Windowed around the peak so the asymmetric
// energy-loss tail doesn't bias the centroid. Returns kFALSE on failure; the
// caller falls back to median/IQR. mpv_seed_override (k=2) seeds and
// constrains the mean near 2x the k=1 anchor; k=1 self-seeds from the maximum.
Bool_t FitGaussianPeak(const std::vector<Float_t> &v, const TString &fname,
                       Double_t &peak_adc, Double_t &sigma_adc, TF1 *&fit_out,
                       Double_t mpv_seed_override = -1.0) {
  fit_out = nullptr;
  if (v.size() < kMinSamplesPerAnchor)
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

  Int_t peak_bin;
  if (mpv_seed_override > 0.0) {
    peak_bin = h.FindBin(mpv_seed_override);
    if (peak_bin < 1)
      peak_bin = 1;
    if (peak_bin > nbins)
      peak_bin = nbins;
  } else {
    peak_bin = h.GetMaximumBin();
  }
  Double_t peak_val = h.GetBinContent(peak_bin);
  Double_t mean_seed = h.GetBinCenter(peak_bin);
  Int_t lo_half = peak_bin;
  while (lo_half > 1 && h.GetBinContent(lo_half) > 0.5 * peak_val)
    lo_half--;
  Int_t hi_half = peak_bin;
  while (hi_half < nbins && h.GetBinContent(hi_half) > 0.5 * peak_val)
    hi_half++;
  Double_t fwhm = h.GetBinCenter(hi_half) - h.GetBinCenter(lo_half);
  if (fwhm <= 0)
    fwhm = (xhi - xlo) / 10.0;
  Double_t bw = h.GetBinWidth(1);
  Double_t sigma_seed = fwhm / 2.355;
  if (sigma_seed < bw)
    sigma_seed = bw;
  // Narrow core fit (this supplies only the resolution sigma + overlay curve;
  // the gain anchor is the sample mean, set by the caller).
  Double_t fit_lo = mean_seed - 1.5 * fwhm;
  Double_t fit_hi = mean_seed + 1.5 * fwhm;
  if (fit_lo < xlo)
    fit_lo = xlo;
  if (fit_hi > xhi)
    fit_hi = xhi;

  TF1 *f = new TF1(fname, "gaus", fit_lo, fit_hi);
  f->SetNpx(1000);
  f->SetParameters(peak_val, mean_seed, sigma_seed);
  f->SetParLimits(1, mean_seed - fwhm, mean_seed + fwhm);
  f->SetParLimits(2, bw, 5.0 * fwhm);
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
    const std::vector<BeamFit2D> &beams,
    std::vector<std::vector<std::vector<Float_t>>> &samples) {
  Int_t n_chans = Int_t(chans.size());
  samples.assign(n_chans, std::vector<std::vector<Float_t>>(kNPeaks));

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
  Int_t leftdE_adc[18], rightdE_adc[18], totaldE_adc[18];
  Int_t cathode_adc = 0;
  tree->SetBranchAddress("LeftdE", leftdE_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("TotaldE", totaldE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    Double_t x = Double_t(totaldE_adc[0]);
    Double_t y = Double_t(totaldE_adc[1]);
    if (x <= 0 || y <= 0)
      continue;

    for (Int_t k = 1; k <= kNPeaks; k++) {
      if (!beams[k - 1].ok)
        continue;
      if (!BeamFitUtils::InEllipseXY(beams[k - 1], x, y, kEllipseNSigmaX,
                                     kEllipseNSigmaY))
        continue;
      for (Int_t i = 0; i < n_chans; i++) {
        if (Long64_t(samples[i][k - 1].size()) >= kSampleCap)
          continue;
        const ChannelCal &c = chans[i];
        Int_t v = 0;
        if (c.side == 'S')
          v = totaldE_adc[c.strip];
        else if (c.side == 'L')
          v = leftdE_adc[c.strip];
        else if (c.side == 'R')
          v = rightdE_adc[c.strip];
        else if (c.side == 'C')
          v = cathode_adc;
        if (v > 0)
          samples[i][k - 1].push_back(Float_t(v));
      }
      break;
    }
  }
  sf->Close();
  delete sf;
}

// Cathode uses median + IQR/1.349 (asymmetric tail not as clean and the user
// prefers to keep cathode on the existing approach). All other channels
// (S guard strips + L/R long anodes) fit a plain Gaussian per k=1, k=2 in a
// window around the peak; the fitted mean is the gain anchor and the sigma is
// the detector resolution. Fall back to median/IQR on fit failure.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<std::vector<Float_t>>> &samples,
                     std::vector<std::vector<TF1 *>> &fits_out,
                     const TString &run_label) {
  Int_t n_chans = Int_t(chans.size());
  fits_out.assign(n_chans, std::vector<TF1 *>(kNPeaks, nullptr));
  for (Int_t i = 0; i < n_chans; i++) {
    ChannelCal &c = chans[i];
    for (Int_t k = 0; k < kNPeaks; k++) {
      std::vector<Float_t> &v = samples[i][k];
      c.n_samples[k] = Long64_t(v.size());
      if (Long64_t(v.size()) < kMinSamplesPerAnchor) {
        c.anchor_adc[k] = 0;
        c.anchor_sigma_adc[k] = 0;
        continue;
      }
      // Anchor = the averaged value of the beam-gated samples (matches the
      // paper / upstream normMUSIC: the beam is normalized to its AVERAGE dE).
      // A Gaussian fit returns the mode regardless of window, which leaves the
      // summed-dE ridge ~2% high; the mean is what "averaged value" means.
      Double_t mean_adc = 0.0;
      for (Int_t j = 0; j < Int_t(v.size()); j++)
        mean_adc += Double_t(v[j]);
      mean_adc /= Double_t(v.size());

      if (c.side == 'C') {
        // Cathode: median + IQR (asymmetric tail, no clean peak).
        c.anchor_adc[k] = Median(v);
        c.anchor_sigma_adc[k] = InterquartileRange(v) / 1.349;
      } else {
        c.anchor_adc[k] = mean_adc;
        // Width + overlay curve from a Gaussian fit; IQR fallback on failure.
        TString fname =
            Form("f_gaus_%s_k%d_%s", c.name.Data(), k + 1, run_label.Data());
        Double_t peak = 0, sig = 0;
        TF1 *fit = nullptr;
        if (FitGaussianPeak(v, fname, peak, sig, fit)) {
          c.anchor_sigma_adc[k] = sig;
          fits_out[i][k] = fit;
        } else {
          c.anchor_sigma_adc[k] = InterquartileRange(v) / 1.349;
        }
      }
    }
    std::cout << "  " << c.name << " anchors[ADC]:";
    for (Int_t k = 0; k < kNPeaks; k++)
      std::cout << " k" << (k + 1) << "=" << c.anchor_adc[k]
                << " sig=" << c.anchor_sigma_adc[k] << " (n=" << c.n_samples[k]
                << ")";
    std::cout << std::endl;
  }
}

// Linear y = a + b*x through (adc_k1, sim_mu) and (adc_k2, 2*sim_mu). Each
// channel is calibrated independently from its own ADC anchors, so no
// cross-channel gain matching is needed. If only the beam (k=1) anchor is
// available (pileup peak too sparse), fall back to a single scale through the
// origin so the channel is still calibrated.
// Standard ("normMUSIC") calibration: scale each beam-dE channel's beam peak to
// the common reference NORM_MUSIC_MEV via a single gain through the origin
// (lin_a = 0, lin_b = NORM_MUSIC_MEV / beam_peak_adc). Calibrated dE is thus
// normalized to a common beam value, not absolute -- which is what the PID /
// summed-dE analysis needs; the reaction energy axis comes from the beam energy
// + stopping-power degradation, not per-strip dE. sim_mu_mev>0 is used only to
// select the channels that have a sim beam anchor (the long sides); the sim
// energy itself is no longer the calibration target.
void FitLinear(ChannelCal &c) {
  if (c.sim_mu_mev <= 0) {
    c.lin_ok = kFALSE;
    return;
  }
  Double_t adc = c.anchor_adc[0];
  if (adc <= 0) {
    c.lin_ok = kFALSE;
    std::cerr << "  " << c.name << ": no beam anchor; no calibration"
              << std::endl;
    return;
  }
  c.lin_a = 0.0;
  c.lin_b = Constants::NORM_MUSIC_MEV / adc;
  c.lin_chi2_ndf = -1;
  c.lin_ok = kTRUE;
}

inline Double_t ApplyLinear(const ChannelCal &c, Double_t adc) {
  return c.lin_a + c.lin_b * adc;
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

  TString out_full = IO::GetRootFilesBaseDir() + "/" + out_subpath;
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
  Float_t lin_b[kMaxChannels] = {0};
  Float_t sim_mu[kMaxChannels] = {0}, sim_sigma[kMaxChannels] = {0};
  Float_t anchor_adc[kNPeaks][kMaxChannels] = {{0}};
  Float_t anchor_sigma[kNPeaks][kMaxChannels] = {{0}};
  Long64_t anchor_n[kNPeaks][kMaxChannels] = {{0}};
  Bool_t lin_ok[kMaxChannels] = {0};
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    lin_b[k] = Float_t(c.lin_b);
    sim_mu[k] = Float_t(c.sim_mu_mev);
    sim_sigma[k] = Float_t(c.sim_sigma_mev);
    for (Int_t p = 0; p < kNPeaks; p++) {
      anchor_adc[p][k] = Float_t(c.anchor_adc[p]);
      anchor_sigma[p][k] = Float_t(c.anchor_sigma_adc[p]);
      anchor_n[p][k] = c.n_samples[p];
    }
    lin_ok[k] = c.lin_ok;
  }
  cal->Branch("LinB", lin_b, Form("LinB[%d]/F", kMaxChannels));
  cal->Branch("LinOK", lin_ok, Form("LinOK[%d]/O", kMaxChannels));
  cal->Branch("SimMu_MeV", sim_mu, Form("SimMu_MeV[%d]/F", kMaxChannels));
  cal->Branch("SimSigma_MeV", sim_sigma,
              Form("SimSigma_MeV[%d]/F", kMaxChannels));
  for (Int_t p = 0; p < kNPeaks; p++) {
    cal->Branch(Form("AnchorMPV_k%d_ADC", p + 1), anchor_adc[p],
                Form("AnchorMPV_k%d_ADC[%d]/F", p + 1, kMaxChannels));
    cal->Branch(Form("AnchorSigma_k%d_ADC", p + 1), anchor_sigma[p],
                Form("AnchorSigma_k%d_ADC[%d]/F", p + 1, kMaxChannels));
    cal->Branch(Form("AnchorN_k%d", p + 1), anchor_n[p],
                Form("AnchorN_k%d[%d]/L", p + 1, kMaxChannels));
  }
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
}

// Per-channel ADC histogram of the samples that fed each multiplicity anchor.
// One file per (channel, k) under <plot_subdir>/beam_peak, named
// beam_peak_k<n>_<channel>.
void SaveBeamPeakChannelHistograms(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<std::vector<Float_t>>> &samples,
    const std::vector<std::vector<TF1 *>> &fits, const TString &plot_subdir) {
  TString subdir = plot_subdir + "/beam_peak";
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    for (Int_t k = 0; k < kNPeaks; k++) {
      const std::vector<Float_t> &v = samples[i][k];
      if (Long64_t(v.size()) < kMinSamplesPerAnchor)
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
      TH1F *h = new TH1F(Form("h_beam_peak_k%d_%s", k + 1, c.name.Data()),
                         Form(";%s #DeltaE [ADC];Counts", c.name.Data()), nbins,
                         Double_t(lo) - pad, Double_t(hi) + pad);
      h->SetDirectory(nullptr);
      for (Int_t j = 0; j < Int_t(v.size()); j++)
        h->Fill(Double_t(v[j]));
      TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h, kBlack);
      TF1 *fit = fits[i][k];
      if (fit) {
        fit->SetLineColor(kRed + 1);
        fit->SetLineWidth(2);
        fit->Draw("L SAME");
      }
      if (Constants::SAVE_PLOTS)
        PlottingUtils::SaveFigure(
            cv, Form("beam_peak_k%d_%s", k + 1, c.name.Data()), subdir,
            PlotSaveOptions::kLINEAR);
      delete cv;
      delete h;
    }
  }
}

void SaveLinearFitPlots(const std::vector<ChannelCal> &chans,
                        const TString &plot_subdir) {
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (!c.lin_ok)
      continue;
    TGraphErrors *g = new TGraphErrors();
    Int_t np = 0;
    Double_t xhi = 0;
    for (Int_t k = 1; k <= kNPeaks; k++) {
      Double_t adc = c.anchor_adc[k - 1];
      if (adc <= 0)
        continue;
      Double_t sigma_adc = c.anchor_sigma_adc[k - 1];
      g->SetPoint(np, adc, Double_t(k) * c.sim_mu_mev);
      g->SetPointError(np, std::max(1.0, sigma_adc), 0.0);
      np++;
      if (adc > xhi)
        xhi = adc;
    }
    if (np == 0) {
      delete g;
      continue;
    }
    xhi *= 1.1;
    if (xhi <= 0)
      xhi = 16384.0;
    PlottingUtils::ConfigureGraph(
        g, kBlack, Form(";%s #DeltaE [ADC];#DeltaE [MeV]", c.name.Data()));
    g->SetMarkerStyle(20);
    TF1 *f = new TF1(Form("f_lin_%s", c.name.Data()), "pol1", 0.0, xhi);
    f->SetParameters(c.lin_a, c.lin_b);
    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    g->Draw("AP");
    f->Draw("L SAME");
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "lin_" + c.name, plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete f;
    delete g;
  }
}

void WriteCalSidecar(const FileSpec &spec,
                     const std::vector<ChannelCal> &chans) {
  TString events_subpath = FileSet::EventsName(spec) + ".root";
  TFile *src = IO::OpenForReading(events_subpath);
  if (!src || src->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << std::endl;
    if (src)
      delete src;
    return;
  }
  TTree *events = static_cast<TTree *>(src->Get("events"));
  if (!events) {
    std::cerr << "No events tree in " << events_subpath << std::endl;
    src->Close();
    delete src;
    return;
  }
  Int_t leftdE_adc[18], rightdE_adc[18], totaldE_adc[18];
  Int_t cathode_adc = 0;
  events->SetBranchAddress("LeftdE", leftdE_adc);
  events->SetBranchAddress("RightdE", rightdE_adc);
  events->SetBranchAddress("TotaldE", totaldE_adc);
  events->SetBranchAddress("Cathode", &cathode_adc);

  const ChannelCal *cal_S[18] = {nullptr};
  const ChannelCal *cal_L[18] = {nullptr};
  const ChannelCal *cal_R[18] = {nullptr};
  const ChannelCal *cal_C = nullptr;
  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    const ChannelCal &c = chans[k];
    if (c.side == 'S')
      cal_S[c.strip] = &c;
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
      cal_L[c.strip] = &c;
    else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
      cal_R[c.strip] = &c;
    else if (c.side == 'C')
      cal_C = &c;
  }

  TString sidecar_subpath = FileSet::CalSidecarName(spec);
  TFile *dst = IO::OpenForWriting(sidecar_subpath, "RECREATE");
  if (!dst || dst->IsZombie()) {
    std::cerr << "Cannot create cal sidecar " << sidecar_subpath << std::endl;
    if (dst)
      delete dst;
    src->Close();
    delete src;
    return;
  }
  dst->cd();
  TTree *cal_tree = new TTree("events_cal", "Per-event calibrated energies");
  Float_t leftdE_MeV[18], rightdE_MeV[18], totaldE_MeV[18];
  Float_t cathode_MeV;
  cal_tree->Branch("LeftdEMeV", leftdE_MeV, "LeftdEMeV[18]/F");
  cal_tree->Branch("RightdEMeV", rightdE_MeV, "RightdEMeV[18]/F");
  cal_tree->Branch("TotaldEMeV", totaldE_MeV, "TotaldEMeV[18]/F");
  cal_tree->Branch("CathodeMeV", &cathode_MeV, "CathodeMeV/F");

  Long64_t n = events->GetEntries();
  std::cout << "  writing cal sidecar " << sidecar_subpath << " (" << n
            << " events)" << std::endl;
  for (Long64_t j = 0; j < n; j++) {
    events->GetEntry(j);
    for (Int_t s = 0; s < 18; s++) {
      leftdE_MeV[s] = 0.0f;
      rightdE_MeV[s] = 0.0f;
      totaldE_MeV[s] = 0.0f;
    }
    for (Int_t s = 0; s < 18; s++) {
      if (s == 0 || s == 17) {
        const ChannelCal *cS = cal_S[s];
        if (cS && cS->lin_ok && totaldE_adc[s] > 0)
          totaldE_MeV[s] = Float_t(ApplyLinear(*cS, Double_t(totaldE_adc[s])));
      } else {
        const ChannelCal *cL = cal_L[s];
        const ChannelCal *cR = cal_R[s];
        if (cL && cL->lin_ok && leftdE_adc[s] > 0)
          leftdE_MeV[s] = Float_t(ApplyLinear(*cL, Double_t(leftdE_adc[s])));
        if (cR && cR->lin_ok && rightdE_adc[s] > 0)
          rightdE_MeV[s] = Float_t(ApplyLinear(*cR, Double_t(rightdE_adc[s])));
        totaldE_MeV[s] = leftdE_MeV[s] + rightdE_MeV[s];
      }
    }
    cathode_MeV = (cal_C && cal_C->lin_ok && cathode_adc > 0)
                      ? Float_t(ApplyLinear(*cal_C, Double_t(cathode_adc)))
                      : 0.0f;
    cal_tree->Fill();
  }
  cal_tree->Write("events_cal", TObject::kOverwrite);
  WriteCalibrationTree(dst, chans);
  dst->Close();
  delete dst;
  src->Close();
  delete src;
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
  TFile *cal_file = FileSet::AttachCalSidecar(tree, spec);
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_mev) {
    sf->Close();
    delete sf;
    if (cal_file) {
      cal_file->Close();
      delete cal_file;
    }
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
  if (cal_file) {
    cal_file->Close();
    delete cal_file;
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
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}

// Overlay (one color per channel, log-y) of ONLY the events used for
// calibration: the k=1 + k=2 anchor samples, converted to MeV via each
// channel's linear fit. Each beam-dE channel shows its beam (k=1) and
// single-pileup (k=2) peaks. Same axes/style as SaveDynamicRangeOverlay but
// restricted to calibration events rather than the full spectrum.
void CalibrateBeam::SaveCalibSampleOverlay(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<std::vector<Float_t>>> &samples,
    const TString &plot_subdir, const TString &file_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = Constants::STRIP_E_MIN_MEV;
  const Double_t emax = Constants::STRIP_E_MAX_MEV;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    const ChannelCal &c = chans[i];
    if (!IsBeamdEChannel(c) || !c.lin_ok)
      continue;
    TString hname =
        Form("h_calibrange_%s_%s", file_label.Data(), c.name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [MeV];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
    for (Int_t k = 0; k < kNPeaks; k++) {
      const std::vector<Float_t> &v = samples[i][k];
      for (Int_t j = 0; j < Int_t(v.size()); j++) {
        Double_t mev = ApplyLinear(c, Double_t(v[j]));
        if (mev > 0)
          h[i]->Fill(mev);
      }
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

  std::vector<BeamFit2D> beams;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beams = FindBeamGateStp1VsStp0(spec, file_label, plot_subdir);
  }
  if (beams.empty() || !beams[0].ok) {
    std::cerr << "  " << file_label << ": Strip1-vs-Strip0 beam gate failed"
              << std::endl;
    return;
  }

  std::vector<ChannelCal> chans = sim_chans;
  std::vector<std::vector<std::vector<Float_t>>> samples;
  CollectAnchorSamplesOneSubfile(spec, chans, beams, samples);
  std::vector<std::vector<TF1 *>> peak_fits;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, peak_fits, file_label);
    SaveBeamPeakChannelHistograms(chans, samples, peak_fits, plot_subdir);
  }
  for (Int_t i = 0; i < Int_t(peak_fits.size()); i++)
    for (Int_t k = 0; k < Int_t(peak_fits[i].size()); k++)
      delete peak_fits[i][k];
  peak_fits.clear();

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    FitLinear(chans[i]);
    if (chans[i].lin_ok)
      std::cout << "  " << chans[i].name << " lin: a=" << chans[i].lin_a
                << ", b=" << chans[i].lin_b
                << ", chi2/ndf=" << chans[i].lin_chi2_ndf << std::endl;
  }

  WriteCalSidecar(spec, chans);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveLinearFitPlots(chans, plot_subdir);
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

void CalibrateBeam::AggregateEresTomlForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  const Int_t n_eres = 35;
  std::vector<std::vector<Double_t>> sigma_per_chan(n_eres);

  std::vector<ChannelCal> tmpl = BuildChannels();
  Int_t n_chans = Int_t(tmpl.size());

  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    TString cal_sub = FileSet::CalSidecarName(specs[s]);
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
    Float_t lin_b[kMaxChannels] = {0};
    Float_t anchor_sigma[kNPeaks][kMaxChannels] = {{0}};
    Bool_t lin_ok[kMaxChannels] = {0};
    t->SetBranchAddress("LinB", lin_b);
    t->SetBranchAddress("LinOK", lin_ok);
    for (Int_t p = 0; p < kNPeaks; p++)
      t->SetBranchAddress(Form("AnchorSigma_k%d_ADC", p + 1), anchor_sigma[p]);
    if (t->GetEntries() < 1) {
      cf->Close();
      delete cf;
      continue;
    }
    t->GetEntry(0);

    for (Int_t i = 0; i < n_chans && i < kMaxChannels; i++) {
      if (!lin_ok[i])
        continue;
      Double_t sig_adc = anchor_sigma[0][i];
      if (sig_adc <= 0 || lin_b[i] <= 0)
        continue;
      Double_t sigma_mev = Double_t(sig_adc) * Double_t(lin_b[i]);
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
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
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
