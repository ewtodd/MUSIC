#include "CalibrateBeamInternal.hpp"

// Per-channel calibrated overlay for one subfile, via AttachCalSidecar (the
// same path downstream macros use). The sidecar must already be on disk.
void SaveDynamicRangeOverlay(const FileSpec &spec,
                             const std::vector<ChannelCal> &chans,
                             const TString &plot_subdir,
                             const TString &file_label) {
  const Int_t kNStrips = 18;
  const Int_t nbins = 300;
  const Double_t emin = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t emax = Constants::cfg.STRIP_DE_MAX_NORMED;
  TH1D *h[kNStrips];
  for (Int_t s = 0; s < kNStrips; s++) {
    h[s] = new TH1D(Form("h_dynrange_%s_S%d", file_label.Data(), s),
                    ";#DeltaE [a.u.];Counts", nbins, emin, emax);
    h[s]->SetDirectory(nullptr);
  }
  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.total[s];
      if (v > 0)
        h[s]->Fill(v);
    }
  }
  sf->Close();
  delete sf;
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t s = 0; s < kNStrips; s++) {
    Double_t m = h[s]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t color = colors[s % Int_t(colors.size())];
    h[s]->SetLineColor(color);
    h[s]->SetLineWidth(2);
    h[s]->SetMaximum(1.15 * y_top);
    h[s]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t s = 0; s < kNStrips; s++)
    leg->AddEntry(h[s], Form("S%d", s), "l");
  leg->Draw();
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t s = 0; s < kNStrips; s++)
    delete h[s];
}

// Log-y overlay (one color per channel) of only the calibration samples
// (beam anchors in a.u. via channel gain); same axes as the dynamic-range one.
void CalibrateBeam::SaveCalibSampleOverlay(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<Float_t>> &samples,
    const TString &plot_subdir, const TString &file_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t emax = Constants::cfg.STRIP_DE_MAX_NORMED;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    const ChannelCal &c = chans[i];
    if (!IsBeamdEChannel(c) || !IsCalibrated(c))
      continue;
    TString hname =
        Form("h_calibrange_%s_%s", file_label.Data(), c.name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [a.u.];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
    const std::vector<Float_t> &v = samples[i];
    for (Int_t j = 0; j < Int_t(v.size()); j++) {
      Double_t val = ApplyCal(c, Double_t(v[j]));
      if (val > 0)
        h[i]->Fill(val);
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
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_calib_events", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}
