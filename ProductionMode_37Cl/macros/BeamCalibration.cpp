#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <mutex>
#include <vector>

namespace BeamCal {

inline TString SimBeamPath(const TString &project_root) {
  TString p = project_root + "/../Remix_37Cl/root_files/traces_37Cl_beam.root";
  return TString(gSystem->ExpandPathName(p.Data()));
}

const Double_t kCathodeWindowSigmas = 3.0;
const Double_t kFitHalfWindowSigmas = 5.0;
const Int_t kFitMinHalfWindowBins = 15;
const Long64_t kMinEntriesForFit = 500;

struct PeakFit {
  Double_t mean = 0;
  Double_t mean_err = 0;
  Double_t sigma = 0;
  Double_t sigma_err = 0;
  Double_t amp = 0;
  Double_t reduced_chi2 = -1;
  Bool_t ok = kFALSE;
};

Bool_t EstimateFitRange(TH1 *h, Double_t &lo, Double_t &hi) {
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

PeakFit PackFitResult(const FitResult &r) {
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

PeakFit FitADCInteractive(TH1 *h, const TString &peak_name,
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
                      /*use_high_exp_tail=*/kTRUE);
  fitter.SetInteractive(kTRUE);
  return PackFitResult(fitter.FitSinglePeak(input_name, peak_name));
}

PeakFit FitSimGaussian(TH1 *h, const TString &peak_name,
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

struct ChannelCal {
  TString name;     // "Strip0", "L1", "R2", ..., "R16", "Strip17", "Cathode"
  Char_t side;      // 'S' / 'L' / 'R' / 'C'
  Int_t strip;      // 0..17 for strip channels; -1 for cathode
  TH1D *hist;       // beam-window-filtered ADC histogram (data)
  PeakFit data_fit; // ADC peak fit
  PeakFit sim_fit;  // Sim Gaussian fit in MeV
  Double_t gain;    // sim_fit.mean / data_fit.mean; 0 if either side unfit
};

struct DataBranches {
  Int_t leftdE[18], rightdE[18], totaldE[18];
  ULong64_t allTimestamps[36];
  UInt_t allFlags[36];
  Int_t hits[36];
  Int_t cathode, grid;
  Bool_t isComplete;
};

void AttachReadBranches(TTree *t, DataBranches &b) {
  t->SetBranchAddress("LeftdE", b.leftdE);
  t->SetBranchAddress("RightdE", b.rightdE);
  t->SetBranchAddress("TotaldE", b.totaldE);
  t->SetBranchAddress("Cathode", &b.cathode);
}

inline Char_t LongSide(Int_t strip) {
  // L_odd / R_even is the dominant ("long") side; the asymmetric split puts
  // the bulk of the per-strip charge on this channel.
  return (strip % 2 == 1) ? 'L' : 'R';
}

inline TString LongName(Int_t strip) {
  return Form("%c%d", LongSide(strip), strip);
}

std::vector<ChannelCal> BuildChannels(const TString &file_label) {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  c.name = "Strip0";
  c.side = 'S';
  c.strip = 0;
  c.hist = new TH1D(Form("h_calraw_%s_Strip0", file_label.Data()),
                    Form("Strip0 (%s);#DeltaE [ADC];Counts", file_label.Data()),
                    Constants::BASELINE_HIST_BINS, 0,
                    Constants::BASELINE_HIST_MAX_ADC);
  chans.push_back(c);
  c.name = "Strip17";
  c.side = 'S';
  c.strip = 17;
  c.hist = new TH1D(
      Form("h_calraw_%s_Strip17", file_label.Data()),
      Form("Strip17 (%s);#DeltaE [ADC];Counts", file_label.Data()),
      Constants::BASELINE_HIST_BINS, 0, Constants::BASELINE_HIST_MAX_ADC);
  chans.push_back(c);
  for (Int_t s = 1; s <= 16; s++) {
    c.name = LongName(s);
    c.side = LongSide(s);
    c.strip = s;
    c.hist = new TH1D(
        Form("h_calraw_%s_%s", file_label.Data(), c.name.Data()),
        Form("%s (%s);#DeltaE [ADC];Counts", c.name.Data(), file_label.Data()),
        Constants::BASELINE_HIST_BINS, 0, Constants::BASELINE_HIST_MAX_ADC);
    chans.push_back(c);
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  c.hist = new TH1D(
      Form("h_calraw_%s_Cathode", file_label.Data()),
      Form("Cathode (%s);#DeltaE [ADC];Counts", file_label.Data()),
      Constants::BASELINE_HIST_BINS, 0, Constants::BASELINE_HIST_MAX_ADC);
  chans.push_back(c);
  return chans;
}

inline Int_t IdxStrip0() { return 0; }
inline Int_t IdxStrip17() { return 1; }
inline Int_t IdxLong(Int_t s) { return 1 + s; }
inline Int_t IdxCathode() { return 18; }

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

Bool_t LoadSimPeaks(const TString &sim_path, std::vector<ChannelCal> &chans,
                    const TString &file_label) {
  TFile *f = TFile::Open(sim_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open sim beam file: " << sim_path << std::endl;
    if (f)
      f->Close();
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("event_MeV"));
  if (!t) {
    std::cerr << "Sim file has no 'event_MeV' tree: " << sim_path
              << " (sim rename pending?)" << std::endl;
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
      ::Fatal("BeamCal::LoadSimPeaks", "Unknown channel side '%c' for '%s'",
              c.side, c.name.Data());
    }
    TString hname = "h_simfit_" + PlottingUtils::GetRandomName();
    TString cmd = expr + ">>" + hname + "(500,0,0)";
    t->Draw(cmd, cut, "goff");
    TH1F *h = static_cast<TH1F *>(gDirectory->Get(hname));
    if (!h) {
      std::cerr << "  sim Draw produced no histogram for " << c.name
                << std::endl;
      continue;
    }
    c.sim_fit = FitSimGaussian(h, "sim_" + c.name, file_label);
    if (!c.sim_fit.ok)
      std::cerr << "  sim fit FAILED for " << c.name << std::endl;
    delete h;
  }
  f->Close();
  return kTRUE;
}

void SaveCalibrationSummary(const std::vector<ChannelCal> &chans,
                            const TString &subdir, const TString &file_label) {
  Int_t n = Int_t(chans.size());
  TGraph *g_data = new TGraph(n);
  TGraph *g_sim = new TGraph(n);
  TGraph *g_gain = new TGraph(n);
  for (Int_t k = 0; k < n; k++) {
    g_data->SetPoint(k, k, chans[k].data_fit.mean);
    g_sim->SetPoint(k, k, chans[k].sim_fit.mean);
    g_gain->SetPoint(k, k, chans[k].gain);
  }
  TCanvas *c1 = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      g_data, kBlack,
      Form("Per-channel beam peak (%s);Channel index;Peak [ADC]",
           file_label.Data()));
  g_data->Draw("ALP");
  PlottingUtils::SaveFigure(c1, "peak_data_adc", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c1;

  TCanvas *c2 = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      g_sim, kBlue + 1,
      Form("Sim beam peak (%s);Channel index;Peak [MeV]", file_label.Data()));
  g_sim->Draw("ALP");
  PlottingUtils::SaveFigure(c2, "peak_sim_mev", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c2;

  TCanvas *c3 = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      g_gain, kRed + 1,
      Form("Per-channel gain (%s);Channel index;gain [MeV/ADC]",
           file_label.Data()));
  g_gain->Draw("ALP");
  PlottingUtils::SaveFigure(c3, "gain", subdir, PlotSaveOptions::kLINEAR);
  delete c3;

  delete g_data;
  delete g_sim;
  delete g_gain;
}

} // namespace BeamCal

void BuildBeamCalibration(std::vector<TString> input_filenames,
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

    std::vector<BeamCal::ChannelCal> chans = BeamCal::BuildChannels(file_label);
    if (!BeamCal::LoadSimPeaks(sim_beam_path, chans, file_label)) {
      std::cerr << "Skipping " << file_label << " (sim peaks unavailable)"
                << std::endl;
      for (Int_t k = 0; k < Int_t(chans.size()); k++)
        delete chans[k].hist;
      file->Close();
      continue;
    }

    BeamCal::DataBranches b;
    BeamCal::AttachReadBranches(tree, b);
    std::cout << "Loading baskets into memory..." << std::endl;
    tree->LoadBaskets();
    Long64_t n = tree->GetEntries();

    TH1D *h_cath_raw = new TH1D(
        Form("h_cathraw_%s", file_label.Data()),
        Form("Cathode raw (%s);Cathode [ADC];Counts", file_label.Data()),
        Constants::BASELINE_HIST_BINS, 0, Constants::BASELINE_HIST_MAX_ADC);
    std::cout << "Pass 1/3: cathode histogram (" << n << " entries)..."
              << std::endl;
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      if (b.cathode > 0)
        h_cath_raw->Fill(Double_t(b.cathode));
      if ((j + 1) % 1000000 == 0)
        std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }

    BeamCal::PeakFit cath =
        BeamCal::FitADCInteractive(h_cath_raw, "Cathode_raw", file_label);
    if (!cath.ok)
      std::cerr << "Cathode fit failed; falling back to no beam mask."
                << std::endl;
    else
      std::cout << "Cathode peak: mu=" << cath.mean << " sigma=" << cath.sigma
                << " chi2/ndf=" << cath.reduced_chi2 << std::endl;

    Double_t cath_lo =
        cath.ok ? cath.mean - BeamCal::kCathodeWindowSigmas * cath.sigma : 0.0;
    Double_t cath_hi =
        cath.ok ? cath.mean + BeamCal::kCathodeWindowSigmas * cath.sigma : 1e12;

    Long64_t beam_event_count = 0;
    std::cout << "Pass 2/3: per-channel histograms..." << std::endl;
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      Double_t cv = Double_t(b.cathode);
      if (!(cv >= cath_lo && cv <= cath_hi))
        continue;
      beam_event_count++;
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        Int_t v = BeamCal::DataValueFor(chans[k], b);
        if (v > 0)
          chans[k].hist->Fill(Double_t(v));
      }
      if ((j + 1) % 1000000 == 0)
        std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }
    std::cout << "Beam-window events: " << beam_event_count << " / " << n
              << " (" << (n > 0 ? 100.0 * beam_event_count / n : 0.0) << "%)"
              << std::endl;

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        BeamCal::ChannelCal &c = chans[k];
        c.data_fit = BeamCal::FitADCInteractive(c.hist, c.name, file_label);
        if (c.data_fit.ok && c.sim_fit.ok && c.data_fit.mean > 0 &&
            c.sim_fit.mean > 0)
          c.gain = c.sim_fit.mean / c.data_fit.mean;
        else
          c.gain = 0.0;
        std::cout << "  " << c.name << ": data_mu=" << c.data_fit.mean
                  << " ADC (sigma=" << c.data_fit.sigma
                  << ", chi2/ndf=" << c.data_fit.reduced_chi2
                  << "), sim_mean=" << c.sim_fit.mean << " MeV"
                  << ", gain=" << c.gain
                  << (c.data_fit.ok ? "" : " [DATA FIT FAILED]")
                  << (c.sim_fit.ok ? "" : " [SIM FIT FAILED]") << std::endl;
      }
      BeamCal::SaveCalibrationSummary(chans, plot_subdir, file_label);
    }

    Float_t leftdE_MeV[18], rightdE_MeV[18], totaldE_MeV[18];
    Float_t cathode_MeV;

    file->cd();

    tree->Branch("LeftdEMeV", leftdE_MeV, "LeftdEMeV[18]/F");
    tree->Branch("RightdEMeV", rightdE_MeV, "RightdEMeV[18]/F");
    tree->Branch("TotaldEMeV", totaldE_MeV, "TotaldEMeV[18]/F");
    tree->Branch("CathodeMeV", &cathode_MeV, "CathodeMeV/F");

    Double_t g_long[18] = {0};
    g_long[0] = chans[BeamCal::IdxStrip0()].gain;
    g_long[17] = chans[BeamCal::IdxStrip17()].gain;
    for (Int_t s = 1; s <= 16; s++)
      g_long[s] = chans[BeamCal::IdxLong(s)].gain;
    Double_t g_C = chans[BeamCal::IdxCathode()].gain;

    std::cout << "Pass 3/3: writing updated event tree..." << std::endl;
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      for (Int_t s = 0; s < 18; s++) {
        leftdE_MeV[s] = 0.0f;
        rightdE_MeV[s] = 0.0f;
        totaldE_MeV[s] = 0.0f;
      }
      totaldE_MeV[0] = Float_t(g_long[0] * Double_t(b.totaldE[0]));
      totaldE_MeV[17] = Float_t(g_long[17] * Double_t(b.totaldE[17]));
      for (Int_t s = 1; s <= 16; s++) {
        if (BeamCal::LongSide(s) == 'L') {
          leftdE_MeV[s] = Float_t(g_long[s] * Double_t(b.leftdE[s]));
          totaldE_MeV[s] = leftdE_MeV[s];
        } else {
          rightdE_MeV[s] = Float_t(g_long[s] * Double_t(b.rightdE[s]));
          totaldE_MeV[s] = rightdE_MeV[s];
        }
      }
      cathode_MeV = Float_t(g_C * Double_t(b.cathode));
      tree->GetBranch("LeftdEMeV")->Fill();
      tree->GetBranch("RightdEMeV")->Fill();
      tree->GetBranch("TotaldEMeV")->Fill();
      tree->GetBranch("CathodeMeV")->Fill();
      if ((j + 1) % 1000000 == 0)
        std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }
    tree->Write("event", TObject::kOverwrite);

    if (TObject *old = file->Get("calibration"))
      old->Delete();
    TTree *cal = new TTree("calibration", "Per-channel beam calibration");

    Float_t mu_adc_arr[19], mu_adc_err_arr[19];
    Float_t sigma_adc_arr[19], sigma_adc_err_arr[19];
    Float_t chi2_adc_arr[19];
    Float_t mean_mev_arr[19], mean_mev_err_arr[19];
    Float_t sigma_mev_arr[19], sigma_mev_err_arr[19];
    Float_t chi2_mev_arr[19];
    Float_t gain_arr[19];
    Bool_t fit_data_ok_arr[19], fit_sim_ok_arr[19];
    for (Int_t k = 0; k < Int_t(chans.size()); k++) {
      const BeamCal::ChannelCal &c = chans[k];
      mu_adc_arr[k] = Float_t(c.data_fit.mean);
      mu_adc_err_arr[k] = Float_t(c.data_fit.mean_err);
      sigma_adc_arr[k] = Float_t(c.data_fit.sigma);
      sigma_adc_err_arr[k] = Float_t(c.data_fit.sigma_err);
      chi2_adc_arr[k] = Float_t(c.data_fit.reduced_chi2);
      mean_mev_arr[k] = Float_t(c.sim_fit.mean);
      mean_mev_err_arr[k] = Float_t(c.sim_fit.mean_err);
      sigma_mev_arr[k] = Float_t(c.sim_fit.sigma);
      sigma_mev_err_arr[k] = Float_t(c.sim_fit.sigma_err);
      chi2_mev_arr[k] = Float_t(c.sim_fit.reduced_chi2);
      gain_arr[k] = Float_t(c.gain);
      fit_data_ok_arr[k] = c.data_fit.ok;
      fit_sim_ok_arr[k] = c.sim_fit.ok;
    }
    cal->Branch("Mu_ADC", mu_adc_arr, "Mu_ADC[19]/F");
    cal->Branch("MuError_ADC", mu_adc_err_arr, "MuError_ADC[19]/F");
    cal->Branch("Sigma_ADC", sigma_adc_arr, "Sigma_ADC[19]/F");
    cal->Branch("SigmaError_ADC", sigma_adc_err_arr, "SigmaError_ADC[19]/F");
    cal->Branch("ReducedChi2_ADC", chi2_adc_arr, "ReducedChi2_ADC[19]/F");
    cal->Branch("Mean_MeV", mean_mev_arr, "Mean_MeV[19]/F");
    cal->Branch("MeanError_MeV", mean_mev_err_arr, "MeanError_MeV[19]/F");
    cal->Branch("Sigma_MeV", sigma_mev_arr, "Sigma_MeV[19]/F");
    cal->Branch("SigmaError_MeV", sigma_mev_err_arr, "SigmaError_MeV[19]/F");
    cal->Branch("ReducedChi2_MeV", chi2_mev_arr, "ReducedChi2_MeV[19]/F");
    cal->Branch("Gain", gain_arr, "Gain[19]/F");
    cal->Branch("FitDataOK", fit_data_ok_arr, "FitDataOK[19]/O");
    cal->Branch("FitSimOK", fit_sim_ok_arr, "FitSimOK[19]/O");
    Float_t cath_window[2] = {Float_t(cath_lo), Float_t(cath_hi)};
    cal->Branch("CathodeWindow", cath_window, "CathodeWindow[2]/F");
    cal->Fill();
    cal->Write("calibration", TObject::kOverwrite);

    TDirectory *hist_dir =
        file->mkdir("beam_calibration_hists", "beam_calibration_hists", kTRUE);
    hist_dir->cd();
    h_cath_raw->Write("h_cathode_raw", TObject::kOverwrite);
    for (Int_t k = 0; k < Int_t(chans.size()); k++)
      chans[k].hist->Write(Form("h_%s", chans[k].name.Data()),
                           TObject::kOverwrite);

    delete h_cath_raw;
    for (Int_t k = 0; k < Int_t(chans.size()); k++)
      delete chans[k].hist;
    file->Close();
    std::cout << "Wrote updated events + calibration for " << file_label
              << std::endl;
  }
}

void BeamCalibration() {
  std::vector<TString> input_filenames, file_labels;
  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    input_filenames.push_back(EventsName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  TString sim_beam_path = BeamCal::SimBeamPath(project_root);
  std::cout << "Sim beam reference: " << sim_beam_path << std::endl;
  BuildBeamCalibration(input_filenames, file_labels, sim_beam_path, kTRUE);
}
