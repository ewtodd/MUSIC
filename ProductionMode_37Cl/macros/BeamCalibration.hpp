#ifndef BEAM_CALIBRATION_HPP
#define BEAM_CALIBRATION_HPP

#include "Constants.hpp"
#include "DrawCut.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mutex>
#include <vector>

namespace BeamCal {

inline TString SimBeamPath(const TString &project_root) {
  TString p = project_root + "/../Remix_37Cl/root_files/traces_37Cl_beam.root";
  return TString(gSystem->ExpandPathName(p.Data()));
}

const Int_t kSimBins = 300;
const Int_t kBeamTagBinsX = 400;
const Int_t kBeamTagBinsY = 400;
const Double_t kBeamTagMaxStripADC = 4000.0;
const Double_t kBeamTagMaxTotalADC = 40000.0;
const TString kBeamCutName = "beam";
const TString kBeamTagHistName = "h2_TotalE_vs_StripE_s1";

// Channel mapping (matches data event-builder convention):
//   Strip0  -> b.totaldE[0]                       (single-ended guard)
//   Strip17 -> b.totaldE[17]                      (single-ended guard)
//   L_odd   -> b.leftdE[s]   for s in {1,3,...,15}  (long-side anode)
//   R_even  -> b.rightdE[s]  for s in {2,4,...,16}  (long-side anode)
//   Cathode -> b.cathode
inline Char_t LongSide(Int_t strip) { return (strip % 2 == 1) ? 'L' : 'R'; }

inline TString LongName(Int_t strip) {
  return Form("%c%d", LongSide(strip), strip);
}

// Mu/sigma extracted from a Gaussian peak fit (data ADC or sim MeV).
struct PeakFit {
  Double_t mean = 0;
  Double_t mean_err = 0;
  Double_t sigma = 0;
  Double_t sigma_err = 0;
  Double_t amp = 0;
  Double_t reduced_chi2 = -1;
  Bool_t ok = kFALSE;
};

struct ChannelCal {
  TString name; // "Strip0", "Strip17", "L1".."R16", "Cathode"
  Char_t side;  // 'S' / 'L' / 'R' / 'C'
  Int_t strip;  // 0..17 for strip channels; -1 for cathode

  TH1D *data_hist_adc; // per-channel ADC histogram, beam-cut events only
  TH1D *sim_hist_mev;  // per-channel sim MeV reference distribution

  PeakFit data_fit; // Gaussian fit to data_hist_adc (ADC)
  PeakFit sim_fit;  // Gaussian fit to sim_hist_mev (MeV)
  Double_t slope;   // sim_fit.mean / data_fit.mean (MeV/ADC)
  Bool_t fit_ok;    // both fits valid AND non-zero data mean
  Long64_t data_n_events;
};

struct DataBranches {
  Int_t leftdE[18], rightdE[18], totaldE[18];
  Int_t cathode;
  Int_t hits[36];
};

inline void AttachReadBranches(TTree *t, DataBranches &b) {
  t->SetBranchAddress("LeftdE", b.leftdE);
  t->SetBranchAddress("RightdE", b.rightdE);
  t->SetBranchAddress("TotaldE", b.totaldE);
  t->SetBranchAddress("Cathode", &b.cathode);
  t->SetBranchAddress("Hits", b.hits);
}

inline Int_t DataValueFor(const ChannelCal &c, const DataBranches &b) {
  if (c.side == 'S')
    return b.totaldE[c.strip];
  if (c.side == 'L')
    return b.leftdE[c.strip];
  if (c.side == 'R')
    return b.rightdE[c.strip];
  if (c.side == 'C')
    return b.cathode;
  return 0;
}

inline Double_t StripTotal(const DataBranches &b, Int_t s) {
  if (!Constants::IGNORE_SHORT_STRIPS || s < 1 || s > 16)
    return Double_t(b.totaldE[s]);
  return (Constants::LongAnodeSide(s) == 'L') ? Double_t(b.leftdE[s])
                                              : Double_t(b.rightdE[s]);
}

inline Double_t EventTotal(const DataBranches &b) {
  Double_t total = 0;
  for (Int_t s = 0; s < 18; s++)
    total += StripTotal(b, s);
  return total;
}

inline Int_t IdxStrip0() { return 0; }
inline Int_t IdxStrip17() { return 1; }
inline Int_t IdxLong(Int_t s) { return 1 + s; }
inline Int_t IdxCathode() { return 18; }
inline Int_t NChannels() { return 19; }

inline std::vector<ChannelCal> BuildChannels(const TString &file_label) {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  c.name = "Strip0";
  c.side = 'S';
  c.strip = 0;
  c.data_hist_adc = nullptr;
  c.sim_hist_mev = nullptr;
  chans.push_back(c);
  c.name = "Strip17";
  c.side = 'S';
  c.strip = 17;
  chans.push_back(c);
  for (Int_t s = 1; s <= 16; s++) {
    c.name = LongName(s);
    c.side = LongSide(s);
    c.strip = s;
    chans.push_back(c);
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  chans.push_back(c);
  return chans;
}

inline Bool_t LoadSimHistograms(const TString &sim_path,
                                std::vector<ChannelCal> &chans,
                                const TString &file_label) {
  TFile *f = TFile::Open(sim_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open sim beam file: " << sim_path << std::endl;
    if (f)
      f->Close();
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "Sim file has no 'events_MeV' tree: " << sim_path << std::endl;
    f->Close();
    return kFALSE;
  }

  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    ChannelCal &c = chans[k];
    TString expr;
    TString cut;
    if (c.side == 'S') {
      expr = Form("TotaldE[%d]", c.strip);
      cut = Form("TotaldE[%d]>0", c.strip);
    } else if (c.side == 'L') {
      expr = Form("LeftdE[%d]", c.strip);
      cut = Form("LeftdE[%d]>0", c.strip);
    } else if (c.side == 'R') {
      expr = Form("RightdE[%d]", c.strip);
      cut = Form("RightdE[%d]>0", c.strip);
    } else if (c.side == 'C') {
      expr = "Cathode";
      cut = "Cathode>0";
    } else {
      ::Fatal("BeamCal::LoadSimHistograms",
              "Unknown channel side '%c' for '%s'", c.side, c.name.Data());
    }

    TString probe_name = "h_simprobe_" + PlottingUtils::GetRandomName();
    t->Draw(expr + ">>" + probe_name + "(200,0,0)", cut, "goff");
    TH1F *probe = static_cast<TH1F *>(gDirectory->Get(probe_name));
    if (!probe || probe->GetEntries() == 0) {
      std::cerr << "  sim Draw produced no histogram for " << c.name
                << std::endl;
      delete probe;
      continue;
    }
    Double_t hi = probe->GetXaxis()->GetXmax() * 1.5;
    delete probe;

    TString hname = "sim_" + file_label + "_" + c.name;
    TH1D *h = new TH1D(hname,
                       Form("Sim %s (%s);#DeltaE [MeV];Counts", c.name.Data(),
                            file_label.Data()),
                       kSimBins, 0.0, hi);

    t->Project(hname, expr, cut);
    h->SetDirectory(nullptr);
    c.sim_hist_mev = h;
  }
  f->Close();
  return kTRUE;
}

const Double_t kFitHalfWindowSigmas = 5.0;
const Int_t kFitMinHalfWindowBins = 15;
const Long64_t kMinEntriesForFit = 500;

inline Bool_t EstimateFitRange(TH1 *h, Double_t &lo, Double_t &hi) {
  if (h->GetEntries() < kMinEntriesForFit)
    return kFALSE;
  Int_t bmax = h->GetMaximumBin();
  Double_t xmax = h->GetBinCenter(bmax);
  Double_t ymax = h->GetBinContent(bmax);
  Double_t bin_w = h->GetBinWidth(bmax);
  Double_t half = ymax * 0.5;
  Int_t bL = bmax;
  while (bL > 1 && h->GetBinContent(bL) > half)
    bL--;
  Int_t bR = bmax;
  while (bR < h->GetNbinsX() && h->GetBinContent(bR) > half)
    bR++;
  Double_t fwhm = h->GetBinCenter(bR) - h->GetBinCenter(bL);
  Double_t sigma = std::max(fwhm / 2.355, bin_w * 2.0);
  Double_t halfwin =
      std::max(kFitHalfWindowSigmas * sigma, kFitMinHalfWindowBins * bin_w);
  lo = std::max(h->GetXaxis()->GetXmin(), xmax - halfwin);
  hi = std::min(h->GetXaxis()->GetXmax(), xmax + halfwin);
  return kTRUE;
}

inline PeakFit PackFitResult(const FitResult &r) {
  PeakFit p;
  if (r.valid && !r.peaks.empty()) {
    p.mean = r.peaks[0].mu;
    p.mean_err = r.peaks[0].mu_error;
    p.sigma = r.peaks[0].sigma;
    p.sigma_err = r.peaks[0].sigma_error;
    p.amp = r.peaks[0].gaus_amplitude;
    p.reduced_chi2 = r.reduced_chi2;
    p.ok = kTRUE;
  }
  return p;
}

inline PeakFit FitADCInteractive(TH1 *h, const TString &peak_name,
                                 const TString &input_name) {
  PeakFit p;
  Double_t lo = 0, hi = 0;
  if (!EstimateFitRange(h, lo, hi))
    return p;
  FittingUtils fitter(h, lo, hi,
                      /*use_flat_background=*/kFALSE,
                      /*use_step=*/kFALSE,
                      /*use_low_exp_tail=*/kFALSE,
                      /*use_low_lin_tail=*/kFALSE,
                      /*use_high_exp_tail=*/kFALSE);
  fitter.SetInteractive(kTRUE);
  return PackFitResult(fitter.FitSinglePeak(input_name, peak_name));
}

inline PeakFit FitSimGaussian(TH1 *h, const TString &peak_name,
                              const TString &input_name) {
  PeakFit p;
  Double_t lo = 0, hi = 0;
  if (!EstimateFitRange(h, lo, hi))
    return p;
  FittingUtils fitter(h, lo, hi,
                      /*use_flat_background=*/kTRUE,
                      /*use_step=*/kFALSE,
                      /*use_low_exp_tail=*/kFALSE,
                      /*use_low_lin_tail=*/kFALSE,
                      /*use_high_exp_tail=*/kFALSE);
  return PackFitResult(fitter.FitSinglePeak(input_name, peak_name));
}

inline void WriteBeamTagHist(TFile *file, TTree *tree,
                             const TString &file_label) {
  TH2D *h2 = new TH2D(kBeamTagHistName,
                      Form(";Strip 1 #DeltaE [ADC];Total #DeltaE [ADC] (%s)",
                           file_label.Data()),
                      kBeamTagBinsX, 0.0, kBeamTagMaxStripADC, kBeamTagBinsY,
                      0.0, kBeamTagMaxTotalADC);

  DataBranches b;
  AttachReadBranches(tree, b);
  Long64_t n = tree->GetEntries();
  std::cout << "[" << file_label << "] Writing " << kBeamTagHistName
            << " (Total vs Strip 1 in ADC) over " << n << " events..."
            << std::endl;
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    if (b.cathode == -1 || b.hits[34] != 1)
      continue;
    h2->Fill(StripTotal(b, 1), EventTotal(b));
  }

  file->cd();
  h2->Write(kBeamTagHistName, TObject::kOverwrite);
  delete h2;
}

inline TCutG *LoadBeamCut(TFile *file) {
  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir)
    return nullptr;
  return static_cast<TCutG *>(cuts_dir->Get(kBeamCutName));
}

inline void SaveCalibrationPlots(const std::vector<ChannelCal> &chans,
                                 const TString &subdir,
                                 const TString &file_label) {
  Int_t n = Int_t(chans.size());
  TGraph *g_slope = new TGraph(n);
  TGraph *g_data_res = new TGraph(n);
  TGraph *g_sim_res = new TGraph(n);
  TGraph *g_eres = new TGraph(n);
  for (Int_t k = 0; k < n; k++) {
    const ChannelCal &c = chans[k];
    g_slope->SetPoint(k, k, c.slope);
    Double_t data_res = (c.data_fit.mean > 0)
                            ? 2.355 * c.data_fit.sigma / c.data_fit.mean
                            : 0.0;
    Double_t sim_res =
        (c.sim_fit.mean > 0) ? 2.355 * c.sim_fit.sigma / c.sim_fit.mean : 0.0;
    Double_t e_res =
        std::sqrt(std::max(0.0, data_res * data_res - sim_res * sim_res));
    g_data_res->SetPoint(k, k, 100.0 * data_res);
    g_sim_res->SetPoint(k, k, 100.0 * sim_res);
    g_eres->SetPoint(k, k, 100.0 * e_res);
  }

  auto save_summary = [&](TGraph *g, Color_t col, const TString &ytitle,
                          const TString &fname) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureGraph(g, col,
                                  Form("%s (%s);Channel index;%s", fname.Data(),
                                       file_label.Data(), ytitle.Data()));
    g->Draw("ALP");
    PlottingUtils::SaveFigure(c, fname, subdir, PlotSaveOptions::kLINEAR);
    delete c;
  };
  save_summary(g_slope, kBlack, "slope [MeV/ADC]", "slope");
  save_summary(g_data_res, kRed + 1, "data FWHM/#mu [%]", "res_data");
  save_summary(g_sim_res, kBlue + 1, "sim FWHM/#mu [%]", "res_sim");
  save_summary(g_eres, kMagenta + 1, "implied e_res FWHM/#mu [%]", "e_res");
  delete g_slope;
  delete g_data_res;
  delete g_sim_res;
  delete g_eres;

  for (Int_t k = 0; k < n; k++) {
    const ChannelCal &c = chans[k];
    if (!c.data_hist_adc)
      continue;

    TCanvas *cv_adc = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(c.data_hist_adc, kBlack);
    PlottingUtils::SaveFigure(cv_adc, "data_adc_" + c.name, subdir,
                              PlotSaveOptions::kLINEAR);
    delete cv_adc;

    if (!c.sim_hist_mev || c.slope <= 0)
      continue;
    TH1D *h_data = new TH1D(
        Form("h_dataMeV_%s_%s", file_label.Data(), c.name.Data()),
        Form("%s sim vs data (slope=%.3e MeV/ADC);#DeltaE [MeV];Counts",
             c.name.Data(), c.slope),
        c.sim_hist_mev->GetNbinsX(), c.sim_hist_mev->GetXaxis()->GetXmin(),
        c.sim_hist_mev->GetXaxis()->GetXmax());
    h_data->SetDirectory(nullptr);
    Int_t nbins_adc = c.data_hist_adc->GetNbinsX();
    for (Int_t b = 1; b <= nbins_adc; b++) {
      Double_t adc = c.data_hist_adc->GetBinCenter(b);
      Double_t cnt = c.data_hist_adc->GetBinContent(b);
      if (cnt > 0)
        h_data->Fill(c.slope * adc, cnt);
    }
    Double_t scale =
        c.sim_hist_mev->Integral() / std::max(1.0, h_data->Integral());
    h_data->Scale(scale);

    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c.sim_hist_mev->SetLineColor(kBlue + 1);
    c.sim_hist_mev->SetLineWidth(2);
    h_data->SetLineColor(kRed + 1);
    h_data->SetLineWidth(2);
    c.sim_hist_mev->Draw("HIST");
    h_data->Draw("HIST SAME");
    PlottingUtils::SaveFigure(cv, "shape_" + c.name, subdir,
                              PlotSaveOptions::kLINEAR);
    delete cv;
    delete h_data;
  }
}

} // namespace BeamCal

inline void BuildBeamCalibration(std::vector<TString> input_filenames,
                                 std::vector<TString> file_labels,
                                 const TString &sim_beam_path,
                                 Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = Int_t(input_filenames.size());
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filenames[i] + ".root";
    TString file_label = file_labels[i];
    TString plot_subdir = "beam_calibration/" + file_label;

    {
      TFile *file = IO::OpenForWriting(input_filepath, "UPDATE");
      if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << input_filepath << std::endl;
        continue;
      }
      TTree *tree = static_cast<TTree *>(file->Get("events"));
      if (!tree) {
        std::cerr << "No events tree in " << input_filepath << std::endl;
        file->Close();
        continue;
      }
      BeamCal::WriteBeamTagHist(file, tree, file_label);
      file->Close();
    }

    Bool_t had_cut = kFALSE;
    {
      TFile *file = IO::OpenForReading(input_filepath);
      had_cut = (BeamCal::LoadBeamCut(file) != nullptr);
      file->Close();
    }
    if (!had_cut) {
      std::cout << "[" << file_label
                << "] No 'beam' cut yet; opening DrawCut interactively..."
                << std::endl;
      {
        std::lock_guard<std::mutex> lock(g_plot_mutex);
        Bool_t was_batch = gROOT->IsBatch();
        DrawCut(BeamCal::kBeamTagHistName, BeamCal::kBeamCutName, file_label);
        gROOT->SetBatch(was_batch);
      }
    }

    TFile *file = IO::OpenForWriting(input_filepath, "UPDATE");
    if (!file || file->IsZombie()) {
      std::cerr << "Cannot re-open " << input_filepath << std::endl;
      continue;
    }
    TTree *tree = static_cast<TTree *>(file->Get("events"));
    if (!tree) {
      std::cerr << "No events tree in " << input_filepath << std::endl;
      file->Close();
      continue;
    }
    TCutG *beam_cut = BeamCal::LoadBeamCut(file);
    if (!beam_cut) {
      std::cerr << "[" << file_label
                << "] 'beam' cut still missing after DrawCut; aborting."
                << " Re-run BeamCalibration once a cut is drawn." << std::endl;
      file->Close();
      continue;
    }
    std::cout << "[" << file_label << "] Using 'beam' TCutG ("
              << beam_cut->GetN() << " vertices)" << std::endl;

    std::vector<BeamCal::ChannelCal> chans = BeamCal::BuildChannels(file_label);
    if (!BeamCal::LoadSimHistograms(sim_beam_path, chans, file_label)) {
      std::cerr << "Skipping " << file_label << " (sim histograms unavailable)"
                << std::endl;
      for (Int_t k = 0; k < Int_t(chans.size()); k++)
        delete chans[k].sim_hist_mev;
      file->Close();
      continue;
    }

    BeamCal::DataBranches b;
    BeamCal::AttachReadBranches(tree, b);
    std::cout << "Loading baskets into memory..." << std::endl;
    tree->LoadBaskets();
    Long64_t n = tree->GetEntries();

    // Pass 1: collect per-channel ADC values for events passing 'beam' cut.
    std::vector<std::vector<Double_t>> data_adc(chans.size());
    Long64_t n_pass = 0;
    std::cout << "Pass 1/2: collecting beam-cut data over " << n
              << " entries..." << std::endl;
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      if (b.cathode == -1 || b.hits[34] != 1)
        continue;
      if (!beam_cut->IsInside(BeamCal::StripTotal(b, 1),
                              BeamCal::EventTotal(b)))
        continue;
      n_pass++;
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        Int_t v = BeamCal::DataValueFor(chans[k], b);
        if (v > 0)
          data_adc[k].push_back(Double_t(v));
      }
      if ((j + 1) % 1000000 == 0)
        std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }
    std::cout << "Beam-cut events: " << n_pass << " / " << n << " ("
              << (n > 0 ? 100.0 * n_pass / n : 0.0) << "%)" << std::endl;

    if (n_pass == 0) {
      std::cerr << "[" << file_label
                << "] No events passed the 'beam' cut; aborting." << std::endl;
      for (Int_t k = 0; k < Int_t(chans.size()); k++)
        delete chans[k].sim_hist_mev;
      file->Close();
      continue;
    }

    for (Int_t k = 0; k < Int_t(chans.size()); k++) {
      BeamCal::ChannelCal &c = chans[k];
      Int_t n_v = Int_t(data_adc[k].size());
      if (n_v == 0) {
        c.data_hist_adc = nullptr;
        continue;
      }
      Double_t hi_adc = (c.side == 'L' || c.side == 'R')
                            ? BeamCal::kBeamTagMaxStripADC
                            : 2 * BeamCal::kBeamTagMaxStripADC;
      c.data_hist_adc =
          new TH1D(Form("h_dataADC_%s_%s", file_label.Data(), c.name.Data()),
                   Form("%s beam-cut data (%lld evts);#DeltaE [ADC];Counts",
                        c.name.Data(), Long64_t(n_v)),
                   200, 0.0, hi_adc);
      c.data_hist_adc->SetDirectory(nullptr);
      for (Int_t j = 0; j < n_v; j++)
        c.data_hist_adc->Fill(data_adc[k][j]);
    }

    Bool_t all_fits_ok = kTRUE;
    std::cout << "Per-channel Gaussian fits:" << std::endl;
    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        BeamCal::ChannelCal &c = chans[k];
        c.data_n_events = Long64_t(data_adc[k].size());
        if (!c.data_hist_adc || !c.sim_hist_mev) {
          c.fit_ok = kFALSE;
          c.slope = 0;
          all_fits_ok = kFALSE;
          std::cout << "  " << c.name << ": SKIPPED (no data or no sim hist)"
                    << std::endl;
          continue;
        }
        c.sim_fit = BeamCal::FitSimGaussian(c.sim_hist_mev, "sim_" + c.name,
                                            file_label);
        c.data_fit =
            BeamCal::FitADCInteractive(c.data_hist_adc, c.name, file_label);

        if (c.sim_fit.ok && c.data_fit.ok && c.data_fit.mean > 0) {
          c.slope = c.sim_fit.mean / c.data_fit.mean;
          c.fit_ok = kTRUE;
        } else {
          c.slope = 0;
          c.fit_ok = kFALSE;
          all_fits_ok = kFALSE;
        }

        Double_t data_res_frac =
            (c.data_fit.mean > 0) ? 2.355 * c.data_fit.sigma / c.data_fit.mean
                                  : 0.0;
        Double_t sim_res_frac = (c.sim_fit.mean > 0)
                                    ? 2.355 * c.sim_fit.sigma / c.sim_fit.mean
                                    : 0.0;
        Double_t e_res_frac = std::sqrt(std::max(
            0.0, data_res_frac * data_res_frac - sim_res_frac * sim_res_frac));
        std::cout << "  " << c.name << ": data_mu=" << c.data_fit.mean
                  << " ADC (FWHM/mu=" << 100.0 * data_res_frac
                  << "%, chi2/ndf=" << c.data_fit.reduced_chi2
                  << "), sim_mu=" << c.sim_fit.mean
                  << " MeV (FWHM/mu=" << 100.0 * sim_res_frac
                  << "%), slope=" << c.slope
                  << ", implied e_res=" << 100.0 * e_res_frac << "%"
                  << (c.data_fit.ok ? "" : " [DATA FIT FAILED]")
                  << (c.sim_fit.ok ? "" : " [SIM FIT FAILED]") << std::endl;
      }
    }

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      BeamCal::SaveCalibrationPlots(chans, plot_subdir, file_label);
    }

    if (!all_fits_ok) {
      std::cerr << "[" << file_label
                << "] One or more channel fits FAILED; skipping Pass 2 "
                << "(MeV branches not written)." << std::endl;
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        delete chans[k].sim_hist_mev;
        delete chans[k].data_hist_adc;
      }
      file->Close();
      continue;
    }

    Float_t leftdE_MeV[18], rightdE_MeV[18], totaldE_MeV[18];
    Float_t cathode_MeV;
    file->cd();
    tree->Branch("LeftdEMeV", leftdE_MeV, "LeftdEMeV[18]/F");
    tree->Branch("RightdEMeV", rightdE_MeV, "RightdEMeV[18]/F");
    tree->Branch("TotaldEMeV", totaldE_MeV, "TotaldEMeV[18]/F");
    tree->Branch("CathodeMeV", &cathode_MeV, "CathodeMeV/F");

    Double_t s_long[18] = {0};
    s_long[0] = chans[BeamCal::IdxStrip0()].slope;
    s_long[17] = chans[BeamCal::IdxStrip17()].slope;
    for (Int_t s = 1; s <= 16; s++) {
      s_long[s] = chans[BeamCal::IdxLong(s)].slope;
    }
    Double_t s_C = chans[BeamCal::IdxCathode()].slope;

    std::cout << "Pass 2/2: writing calibrated MeV branches..." << std::endl;
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      for (Int_t s = 0; s < 18; s++) {
        leftdE_MeV[s] = 0.0f;
        rightdE_MeV[s] = 0.0f;
        totaldE_MeV[s] = 0.0f;
      }
      totaldE_MeV[0] = Float_t(s_long[0] * Double_t(b.totaldE[0]));
      totaldE_MeV[17] = Float_t(s_long[17] * Double_t(b.totaldE[17]));
      for (Int_t s = 1; s <= 16; s++) {
        if (BeamCal::LongSide(s) == 'L') {
          leftdE_MeV[s] = Float_t(s_long[s] * Double_t(b.leftdE[s]));
          totaldE_MeV[s] = leftdE_MeV[s];
        } else {
          rightdE_MeV[s] = Float_t(s_long[s] * Double_t(b.rightdE[s]));
          totaldE_MeV[s] = rightdE_MeV[s];
        }
      }
      cathode_MeV = Float_t(s_C * Double_t(b.cathode));
      tree->GetBranch("LeftdEMeV")->Fill();
      tree->GetBranch("RightdEMeV")->Fill();
      tree->GetBranch("TotaldEMeV")->Fill();
      tree->GetBranch("CathodeMeV")->Fill();
      if ((j + 1) % 1000000 == 0)
        std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }
    tree->Write("events", TObject::kOverwrite);

    if (TObject *old = file->Get("calibration"))
      old->Delete();
    TTree *cal = new TTree("calibration", "Per-channel beam calibration");

    const Int_t N = BeamCal::NChannels();
    Float_t slope_arr[19];
    Float_t data_mu_arr[19], data_mu_err_arr[19];
    Float_t data_sigma_arr[19], data_sigma_err_arr[19];
    Float_t data_chi2_arr[19];
    Float_t sim_mu_arr[19], sim_mu_err_arr[19];
    Float_t sim_sigma_arr[19], sim_sigma_err_arr[19];
    Float_t sim_chi2_arr[19];
    Float_t sigma_data_mev_arr[19];
    Float_t e_res_mev_arr[19];
    Float_t data_res_arr[19], sim_res_arr[19], e_res_frac_arr[19];
    Long64_t data_n_arr[19];
    Bool_t data_fit_ok_arr[19], sim_fit_ok_arr[19], fit_ok_arr[19];
    for (Int_t k = 0; k < N; k++) {
      const BeamCal::ChannelCal &c = chans[k];
      slope_arr[k] = Float_t(c.slope);
      data_mu_arr[k] = Float_t(c.data_fit.mean);
      data_mu_err_arr[k] = Float_t(c.data_fit.mean_err);
      data_sigma_arr[k] = Float_t(c.data_fit.sigma);
      data_sigma_err_arr[k] = Float_t(c.data_fit.sigma_err);
      data_chi2_arr[k] = Float_t(c.data_fit.reduced_chi2);
      sim_mu_arr[k] = Float_t(c.sim_fit.mean);
      sim_mu_err_arr[k] = Float_t(c.sim_fit.mean_err);
      sim_sigma_arr[k] = Float_t(c.sim_fit.sigma);
      sim_sigma_err_arr[k] = Float_t(c.sim_fit.sigma_err);
      sim_chi2_arr[k] = Float_t(c.sim_fit.reduced_chi2);
      Double_t sigma_data_mev = c.data_fit.sigma * c.slope;
      sigma_data_mev_arr[k] = Float_t(sigma_data_mev);
      e_res_mev_arr[k] = Float_t(
          std::sqrt(std::max(0.0, sigma_data_mev * sigma_data_mev -
                                      c.sim_fit.sigma * c.sim_fit.sigma)));
      // FWHM/mu fractional resolutions (the cross-channel-comparable form).
      Double_t data_res_frac = (c.data_fit.mean > 0)
                                   ? 2.355 * c.data_fit.sigma / c.data_fit.mean
                                   : 0.0;
      Double_t sim_res_frac =
          (c.sim_fit.mean > 0) ? 2.355 * c.sim_fit.sigma / c.sim_fit.mean : 0.0;
      data_res_arr[k] = Float_t(data_res_frac);
      sim_res_arr[k] = Float_t(sim_res_frac);
      e_res_frac_arr[k] = Float_t(std::sqrt(std::max(
          0.0, data_res_frac * data_res_frac - sim_res_frac * sim_res_frac)));
      data_n_arr[k] = c.data_n_events;
      data_fit_ok_arr[k] = c.data_fit.ok;
      sim_fit_ok_arr[k] = c.sim_fit.ok;
      fit_ok_arr[k] = c.fit_ok;
    }
    cal->Branch("Slope", slope_arr, "Slope[19]/F");
    cal->Branch("DataMu_ADC", data_mu_arr, "DataMu_ADC[19]/F");
    cal->Branch("DataMuErr_ADC", data_mu_err_arr, "DataMuErr_ADC[19]/F");
    cal->Branch("DataSigma_ADC", data_sigma_arr, "DataSigma_ADC[19]/F");
    cal->Branch("DataSigmaErr_ADC", data_sigma_err_arr,
                "DataSigmaErr_ADC[19]/F");
    cal->Branch("DataChi2NDF", data_chi2_arr, "DataChi2NDF[19]/F");
    cal->Branch("SimMu_MeV", sim_mu_arr, "SimMu_MeV[19]/F");
    cal->Branch("SimMuErr_MeV", sim_mu_err_arr, "SimMuErr_MeV[19]/F");
    cal->Branch("SimSigma_MeV", sim_sigma_arr, "SimSigma_MeV[19]/F");
    cal->Branch("SimSigmaErr_MeV", sim_sigma_err_arr, "SimSigmaErr_MeV[19]/F");
    cal->Branch("SimChi2NDF", sim_chi2_arr, "SimChi2NDF[19]/F");
    cal->Branch("DataSigma_MeV", sigma_data_mev_arr, "DataSigma_MeV[19]/F");
    cal->Branch("ImpliedERes_MeV", e_res_mev_arr, "ImpliedERes_MeV[19]/F");
    cal->Branch("DataRes_FWHM_over_mu", data_res_arr,
                "DataRes_FWHM_over_mu[19]/F");
    cal->Branch("SimRes_FWHM_over_mu", sim_res_arr,
                "SimRes_FWHM_over_mu[19]/F");
    cal->Branch("ImpliedERes_FWHM_over_mu", e_res_frac_arr,
                "ImpliedERes_FWHM_over_mu[19]/F");
    cal->Branch("DataNEvents", data_n_arr, "DataNEvents[19]/L");
    cal->Branch("DataFitOK", data_fit_ok_arr, "DataFitOK[19]/O");
    cal->Branch("SimFitOK", sim_fit_ok_arr, "SimFitOK[19]/O");
    cal->Branch("FitOK", fit_ok_arr, "FitOK[19]/O");
    cal->Fill();
    cal->Write("calibration", TObject::kOverwrite);

    TDirectory *hist_dir =
        file->mkdir("beam_calibration_hists", "beam_calibration_hists", kTRUE);
    hist_dir->cd();
    for (Int_t k = 0; k < N; k++) {
      const BeamCal::ChannelCal &c = chans[k];
      if (c.sim_hist_mev)
        c.sim_hist_mev->Write(Form("h_sim_%s", c.name.Data()),
                              TObject::kOverwrite);
    }

    for (Int_t k = 0; k < Int_t(chans.size()); k++) {
      delete chans[k].sim_hist_mev;
      delete chans[k].data_hist_adc;
    }
    file->Close();
    std::cout << "Wrote Gaussian-fit calibration for " << file_label
              << std::endl;
  }
}

#endif
