#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <vector>

namespace GainDrift {

const Int_t N_CHUNKS = 20;
const Int_t HIST_BINS = 400;
const Double_t HIST_MIN = 0.0;
const Double_t HIST_MAX = 10000;

// Fit windows for the lowest peak, by channel group. The skew-normal fit
// looks for the max bin inside [lo, hi] and seeds mu there.
const Double_t FIT_LO_S0 = 1200.0;
const Double_t FIT_HI_S0 = 6000.0;
const Double_t FIT_LO_LR = 500.0;
const Double_t FIT_HI_LR = 1200.0;
const Double_t FIT_LO_S17 = 700.0;
const Double_t FIT_HI_S17 = 1800.0;

// Skew-normal + linear background.
//   f(x) = A * exp(-z^2/2) * (1 + erf(alpha*z/sqrt(2))) + bkg_const +
//   bkg_slope*x
// where z = (x - xi) / omega.
// Reduces to a Gaussian when alpha = 0; alpha controls one-sided skew.
//   par[0] = A          (amplitude scale; peak height ~ 2*A near alpha=0)
//   par[1] = xi         (location)
//   par[2] = omega      (scale, > 0)
//   par[3] = alpha      (shape; 0 = symmetric)
//   par[4] = bkg const
//   par[5] = bkg slope
inline double SkewNormalPlusBkg(double *x, double *p) {
  double bkg = p[4] + p[5] * x[0];
  if (p[2] <= 0)
    return bkg;
  double z = (x[0] - p[1]) / p[2];
  return p[0] * TMath::Exp(-0.5 * z * z) *
             (1.0 + TMath::Erf(p[3] * z / TMath::Sqrt(2.0))) +
         bkg;
}

struct PeakFit {
  Double_t mu;
  Double_t sigma;
  Bool_t valid;
};

// Linear fit (mu vs run-time-seconds) per channel, used to apply a
// multiplicative gain-drift correction.
struct ChannelDrift {
  Double_t intercept; // mu at t = 0 [ADC]
  Double_t slope;     // mu drift rate [ADC/s]
  Bool_t valid;
};

inline PeakFit FitLowestPeak(TH1D *h, Double_t fit_lo, Double_t fit_hi,
                             const TString &save_name = "",
                             const TString &save_subdir = "") {
  PeakFit pf{0, 0, kFALSE};

  Int_t lo_bin = h->FindBin(fit_lo);
  Int_t hi_bin = h->FindBin(fit_hi);
  if (lo_bin < 1)
    lo_bin = 1;
  if (hi_bin > h->GetNbinsX())
    hi_bin = h->GetNbinsX();

  Int_t max_bin = lo_bin;
  Double_t max_count = h->GetBinContent(lo_bin);
  for (Int_t b = lo_bin + 1; b <= hi_bin; b++) {
    if (h->GetBinContent(b) > max_count) {
      max_count = h->GetBinContent(b);
      max_bin = b;
    }
  }

  Double_t mu_seed = h->GetBinCenter(max_bin);
  Double_t omega_seed = (fit_hi - fit_lo) * 0.05;
  Double_t bkg_seed = h->GetBinContent(lo_bin);
  Double_t A_seed = TMath::Max(1.0, max_count - bkg_seed);

  TF1 f(PlottingUtils::GetRandomName(), SkewNormalPlusBkg, fit_lo, fit_hi, 6);
  f.SetParameters(A_seed, mu_seed, omega_seed, 0.0, bkg_seed, 0.0);
  f.SetParLimits(1, fit_lo, fit_hi);
  f.SetParLimits(2, omega_seed * 0.05, (fit_hi - fit_lo) * 0.5);
  f.SetParLimits(3, -10.0, 10.0);

  TFitResultPtr r = h->Fit(&f, "QSRN");
  pf.valid = (r.Get() && r->IsValid());
  pf.mu = f.GetParameter(1);
  pf.sigma = TMath::Abs(f.GetParameter(2));

  if (!save_name.IsNull()) {
    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h, kBlack);
    h->GetXaxis()->SetRangeUser(fit_lo, fit_hi);
    f.SetLineColor(kRed + 1);
    f.SetLineWidth(2);
    f.SetNpx(1000);
    f.DrawCopy("SAME");

    TLatex *txt = PlottingUtils::AddText(
        Form("#mu = %.1f, #sigma = %.1f, #alpha = %.2f%s", pf.mu, pf.sigma,
             f.GetParameter(3), pf.valid ? "" : " (FAILED)"),
        0.88, 0.85);

    PlottingUtils::SaveFigure(canvas, save_name, save_subdir,
                              PlotSaveOptions::kLOG);
    delete txt;
    delete canvas;
  }

  return pf;
}

// Long-tap channels: side=='S' uses TotaldE (Strip0/17, single readout),
// 'L' uses LeftdE (odd strips), 'R' uses RightdE (even strips). Matches the
// rule in BaselineNormalization.
struct Channel {
  TString name;
  Int_t strip;
  Char_t side;
};

inline std::vector<Channel> LongChannels() {
  std::vector<Channel> chans;
  chans.push_back({"Strip0", 0, 'S'});
  for (Int_t s = 1; s <= 15; s += 2)
    chans.push_back({Form("L%d", s), s, 'L'});
  for (Int_t s = 2; s <= 16; s += 2)
    chans.push_back({Form("R%d", s), s, 'R'});
  chans.push_back({"Strip17", 17, 'S'});
  return chans;
}

inline ULong64_t MinNonzeroTimestamp(const ULong64_t *all_ts, Int_t n) {
  ULong64_t best = 0;
  for (Int_t k = 0; k < n; k++) {
    if (all_ts[k] == 0)
      continue;
    if (best == 0 || all_ts[k] < best)
      best = all_ts[k];
  }
  return best;
}

void PlotOverlay(const std::vector<Channel> &chans,
                 const std::vector<std::vector<TH1D *>> &hists,
                 const TString &plot_subdir) {
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
    canvas->SetRightMargin(0.20);

    Double_t ymax = 0;
    for (Int_t k = 0; k < N_CHUNKS; k++) {
      Double_t m = hists[i][k]->GetMaximum();
      if (m > ymax)
        ymax = m;
    }

    for (Int_t k = 0; k < N_CHUNKS; k++) {
      Int_t color = colors[k % Int_t(colors.size())];
      PlottingUtils::ConfigureHistogram(hists[i][k], color);
      hists[i][k]->SetMaximum(1.15 * ymax);
      hists[i][k]->Draw(k == 0 ? "HIST" : "HIST SAME");
    }

    TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.15, 0.92);
    leg->SetTextSize(TMath::Min(20, 320 / TMath::Max(N_CHUNKS, 1)));
    for (Int_t k = 0; k < N_CHUNKS; k++)
      leg->AddEntry(hists[i][k], Form("chunk %d/%d", k + 1, N_CHUNKS), "l");
    leg->Draw();

    PlottingUtils::SaveFigure(canvas, "drift_" + chans[i].name, plot_subdir,
                              PlotSaveOptions::kLOG);
    delete leg;
    delete canvas;
  }
}

std::vector<ChannelDrift>
PlotDriftVsTime(const std::vector<Channel> &chans,
                const std::vector<std::vector<TH1D *>> &hists,
                const std::vector<Double_t> &chunk_t_mid,
                const std::vector<Double_t> &chunk_t_half,
                const TString &file_label, const TString &plot_subdir) {
  std::vector<ChannelDrift> drifts(chans.size(), {0, 0, kFALSE});

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    TGraphErrors *graph = new TGraphErrors();
    graph->SetTitle(Form("%s lowest-peak drift (%s);Run time [s];"
                         "Peak position [ADC]",
                         chans[i].name.Data(), file_label.Data()));

    for (Int_t k = 0; k < N_CHUNKS; k++) {
      TH1D *h = hists[i][k];
      if (h->GetEntries() < 50) {
        std::cerr << "Skipping " << chans[i].name << " chunk " << k << " (only "
                  << h->GetEntries() << " entries)" << std::endl;
        continue;
      }

      Double_t fit_lo, fit_hi;
      if (chans[i].side == 'S' && chans[i].strip == 0) {
        fit_lo = FIT_LO_S0;
        fit_hi = FIT_HI_S0;
      } else if (chans[i].side == 'S' && chans[i].strip == 17) {
        fit_lo = FIT_LO_S17;
        fit_hi = FIT_HI_S17;
      } else {
        fit_lo = FIT_LO_LR;
        fit_hi = FIT_HI_LR;
      }

      TString fit_name = Form("fit_%s_chunk%02d", chans[i].name.Data(), k);
      TString fit_subdir = plot_subdir + "/fits/" + chans[i].name;
      PeakFit pf = FitLowestPeak(h, fit_lo, fit_hi, fit_name, fit_subdir);
      if (!pf.valid) {
        std::cerr << "Fit failed for " << chans[i].name << " chunk " << k
                  << std::endl;
        continue;
      }

      Int_t n = graph->GetN();
      graph->SetPoint(n, chunk_t_mid[k], pf.mu);
      graph->SetPointError(n, chunk_t_half[k], pf.sigma);
    }

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureGraph(graph, kBlue + 1);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.9);
    graph->Draw("AP");

    TLatex *txt = nullptr;
    if (graph->GetN() >= 2) {
      TF1 *line = new TF1(PlottingUtils::GetRandomName(), "pol1");
      line->SetLineColor(kRed + 1);
      line->SetLineWidth(2);
      TFitResultPtr lr = graph->Fit(line, "QS");
      if (lr.Get() && lr->IsValid()) {
        Double_t slope = line->GetParameter(1);
        Double_t slope_err = line->GetParError(1);
        drifts[i].intercept = line->GetParameter(0);
        drifts[i].slope = slope;
        drifts[i].valid = kTRUE;
        txt = PlottingUtils::AddText(
            Form("slope = %.3g #pm %.3g ADC/s", slope, slope_err), 0.88, 0.85);
      }
    }

    PlottingUtils::SaveFigure(canvas, "drift_vs_time_" + chans[i].name,
                              plot_subdir, PlotSaveOptions::kLINEAR);
    if (txt)
      delete txt;
    delete canvas;
    delete graph;
  }

  return drifts;
}

void DriftOneFile(const TString &events_name, const TString &file_label) {
  TString filepath = events_name + ".root";
  TFile *file = IO::OpenForReading(filepath);
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << filepath << std::endl;
    return;
  }
  TTree *tree = static_cast<TTree *>(file->Get("event"));
  if (!tree) {
    std::cerr << "No 'event' tree in " << filepath << std::endl;
    file->Close();
    return;
  }

  Int_t leftdE[18], rightdE[18], totaldE[18];
  ULong64_t allTs[36];
  tree->SetBranchAddress("LeftdE", leftdE);
  tree->SetBranchAddress("RightdE", rightdE);
  tree->SetBranchAddress("TotaldE", totaldE);
  tree->SetBranchAddress("AllTimestamps", allTs);

  Long64_t n_entries = tree->GetEntries();
  if (n_entries < N_CHUNKS) {
    std::cerr << "Only " << n_entries << " entries in " << filepath
              << ", need at least " << N_CHUNKS << std::endl;
    file->Close();
    return;
  }

  std::vector<Channel> chans = LongChannels();

  std::vector<std::vector<TH1D *>> hists(
      chans.size(), std::vector<TH1D *>(N_CHUNKS, nullptr));
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    for (Int_t c = 0; c < N_CHUNKS; c++) {
      TString hname = Form("h_drift_%s_%s_chunk%d", file_label.Data(),
                           chans[i].name.Data(), c);
      TString title = Form("%s gain drift (%s);#DeltaE [ADC];Counts",
                           chans[i].name.Data(), file_label.Data());
      hists[i][c] = new TH1D(hname, title, HIST_BINS, HIST_MIN, HIST_MAX);
    }
  }

  std::vector<ULong64_t> chunk_t_first(N_CHUNKS, 0);
  std::vector<ULong64_t> chunk_t_last(N_CHUNKS, 0);
  ULong64_t run_t0 = 0;

  Long64_t chunk_size = n_entries / N_CHUNKS;
  std::cout << "Filling " << n_entries << " events into " << N_CHUNKS
            << " time chunks (~" << chunk_size << " events each)" << std::endl;

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    Int_t chunk = TMath::Min(Int_t(j / chunk_size), N_CHUNKS - 1);

    ULong64_t et = MinNonzeroTimestamp(allTs, 36);
    if (et > 0) {
      if (run_t0 == 0)
        run_t0 = et;
      if (chunk_t_first[chunk] == 0)
        chunk_t_first[chunk] = et;
      chunk_t_last[chunk] = et;
    }

    for (Int_t i = 0; i < Int_t(chans.size()); i++) {
      const Channel &c = chans[i];
      Int_t v = (c.side == 'S')
                    ? totaldE[c.strip]
                    : (c.side == 'L' ? leftdE[c.strip] : rightdE[c.strip]);
      if (v > 0)
        hists[i][chunk]->Fill(Double_t(v));
    }
  }

  std::vector<Double_t> chunk_t_mid(N_CHUNKS, 0);
  std::vector<Double_t> chunk_t_half(N_CHUNKS, 0);
  const Double_t PS_PER_S = 1e12;
  for (Int_t k = 0; k < N_CHUNKS; k++) {
    if (chunk_t_first[k] == 0 || chunk_t_last[k] == 0 || run_t0 == 0) {
      chunk_t_mid[k] = k + 0.5;
      chunk_t_half[k] = 0.5;
      continue;
    }
    Double_t t_first = (chunk_t_first[k] - run_t0) / PS_PER_S;
    Double_t t_last = (chunk_t_last[k] - run_t0) / PS_PER_S;
    chunk_t_mid[k] = 0.5 * (t_first + t_last);
    chunk_t_half[k] = 0.5 * (t_last - t_first);
  }

  TString plot_subdir = "gain_drift/" + file_label;

  PlotOverlay(chans, hists, plot_subdir);
  std::vector<ChannelDrift> drifts = PlotDriftVsTime(
      chans, hists, chunk_t_mid, chunk_t_half, file_label, plot_subdir);

  // Apply per-channel linear drift correction:
  //   corrected = raw * mu_ref / (intercept + slope * t)
  // with mu_ref taken at the run midpoint so corrected values stay in roughly
  // the same range as the raw histograms.
  std::vector<std::vector<TH1D *>> corrected_hists(
      chans.size(), std::vector<TH1D *>(N_CHUNKS, nullptr));
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    for (Int_t c = 0; c < N_CHUNKS; c++) {
      TString hname = Form("h_corr_%s_%s_chunk%d", file_label.Data(),
                           chans[i].name.Data(), c);
      TString title = Form("%s gain drift corrected (%s);#DeltaE [ADC];Counts",
                           chans[i].name.Data(), file_label.Data());
      corrected_hists[i][c] = new TH1D(hname, title, HIST_BINS, HIST_MIN, HIST_MAX);
    }
  }

  Double_t t_run_mid = 0.5 * (chunk_t_mid.front() + chunk_t_mid.back());
  std::vector<Double_t> mu_ref(chans.size(), 0);
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    if (drifts[i].valid)
      mu_ref[i] = drifts[i].intercept + drifts[i].slope * t_run_mid;
  }

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    Int_t chunk = TMath::Min(Int_t(j / chunk_size), N_CHUNKS - 1);

    ULong64_t et = MinNonzeroTimestamp(allTs, 36);
    Double_t t_s = (et > 0 && run_t0 > 0) ? (et - run_t0) / PS_PER_S : -1;

    for (Int_t i = 0; i < Int_t(chans.size()); i++) {
      const Channel &c = chans[i];
      Int_t v = (c.side == 'S')
                    ? totaldE[c.strip]
                    : (c.side == 'L' ? leftdE[c.strip] : rightdE[c.strip]);
      if (v <= 0)
        continue;

      Double_t corrected = Double_t(v);
      if (drifts[i].valid && t_s >= 0 && mu_ref[i] > 0) {
        Double_t mu_t = drifts[i].intercept + drifts[i].slope * t_s;
        if (mu_t > 0)
          corrected = v * mu_ref[i] / mu_t;
      }
      corrected_hists[i][chunk]->Fill(corrected);
    }
  }

  PlotOverlay(chans, corrected_hists, plot_subdir + "/corrected");

  for (Int_t i = 0; i < Int_t(chans.size()); i++)
    for (Int_t k = 0; k < N_CHUNKS; k++) {
      delete hists[i][k];
      delete corrected_hists[i][k];
    }

  file->Close();
  std::cout << "Saved gain-drift plots under plots/" << plot_subdir
            << std::endl;
}

} // namespace GainDrift

void CheckGainDrift() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString events_name = EventsName(specs[k]);
    TString file_label = FileLabel(specs[k]);
    std::cout << "Checking gain drift in " << events_name << std::endl;
    GainDrift::DriftOneFile(events_name, file_label);
  }
}
