#include "Constants.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

void SiCalibration() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *out_file = IO::OpenForWriting(Constants::FIT_RESULTS_FILE);

  TTree *results = new TTree("FitResults", "Per-run Si peak fit results");
  Int_t run_number;
  Double_t mu, mu_err, sigma, sigma_err, reduced_chi2;
  results->Branch("RunNumber", &run_number, "RunNumber/I");
  results->Branch("Mu", &mu, "Mu/D");
  results->Branch("MuError", &mu_err, "MuError/D");
  results->Branch("Sigma", &sigma, "Sigma/D");
  results->Branch("SigmaError", &sigma_err, "SigmaError/D");
  results->Branch("ReducedChi2", &reduced_chi2, "ReducedChi2/D");

  for (Int_t run : Constants::RUN_NUMBERS) {
    TString filepath = Constants::GetRunFilepath(run);
    TFile *file = IO::OpenForReading(filepath);
    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file for run " << run << std::endl;
      continue;
    }

    TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
    if (!tree) {
      std::cerr << "Error getting tree for run " << run << std::endl;
      file->Close();
      delete file;
      continue;
    }

    UShort_t channel, energy;
    tree->SetBranchAddress("Channel", &channel);
    tree->SetBranchAddress("Energy", &energy);

    TString hist_name = Form("Run_%d", run);
    TH1F *hist =
        new TH1F(hist_name,
                 Form("Run %d ; ADC Channel; Counts / %g ADC", run,
                      Constants::BIN_WIDTH_ADC),
                 Constants::N_ADC_BINS, Constants::ADC_MIN, Constants::ADC_MAX);

    Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; i++) {
      tree->GetEntry(i);
      if (channel == Constants::SI_DETECTOR_CHANNEL) {
        hist->Fill(energy);
      }
    }

    Double_t mu_guess = Constants::MU_GUESSES.at(run);
    Double_t fit_min = TMath::Max(Constants::ADC_MIN,
                                  mu_guess - Constants::FIT_HALF_WIDTH_ADC);
    Double_t fit_max = TMath::Min(Constants::ADC_MAX,
                                  mu_guess + Constants::FIT_HALF_WIDTH_ADC);

    FittingUtils fitter(hist, fit_min, fit_max,
                        /*use_flat_background=*/kTRUE,
                        /*use_step=*/kFALSE,
                        /*use_low_exp_tail=*/kFALSE,
                        /*use_low_lin_tail=*/kFALSE,
                        /*use_high_exp_tail=*/kTRUE);
    fitter.SetInteractive(kTRUE);

    FitResult fit_result = fitter.FitSinglePeak(hist_name, "alpha");

    run_number = run;
    if (fit_result.valid && !fit_result.peaks.empty()) {
      mu = fit_result.peaks[0].mu;
      mu_err = fit_result.peaks[0].mu_error;
      sigma = fit_result.peaks[0].sigma;
      sigma_err = fit_result.peaks[0].sigma_error;
      reduced_chi2 = fit_result.reduced_chi2;
    } else {
      std::cerr << "WARNING: FitSinglePeak failed for run " << run << std::endl;
      mu = mu_guess;
      mu_err = 0;
      sigma = 200;
      sigma_err = 0;
      reduced_chi2 = -1;
    }
    results->Fill();

    out_file->cd();
    hist->Write(hist_name, TObject::kOverwrite);
    if (TF1 *fit_func = fitter.GetFitFunction()) {
      fit_func->Write(Form("FitFunc_Run_%d", run), TObject::kOverwrite);
    }

    delete hist;
    file->Close();
    delete file;
  }

  out_file->cd();
  results->Write("FitResults", TObject::kOverwrite);
  out_file->Close();
  delete out_file;
}
