#include "SiFits.hpp"

void SiFits::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *out_file = IO::OpenForWriting(SiCalib::FIT_RESULTS_FILE);

  TTree *results = new TTree("FitResults", "Per-run Si peak fit results");
  Int_t run_number;
  Double_t mu, mu_err, sigma, sigma_err, reduced_chi2;
  results->Branch("RunNumber", &run_number, "RunNumber/I");
  results->Branch("Mu", &mu, "Mu/D");
  results->Branch("MuError", &mu_err, "MuError/D");
  results->Branch("Sigma", &sigma, "Sigma/D");
  results->Branch("SigmaError", &sigma_err, "SigmaError/D");
  results->Branch("ReducedChi2", &reduced_chi2, "ReducedChi2/D");

  for (size_t i = 0; i < SiCalib::RUN_NUMBERS.size(); i++) {
    Int_t run = SiCalib::RUN_NUMBERS[i];
    TString filepath = SiCalib::GetRunFilepath(run);
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
                      SiCalib::BIN_WIDTH_ADC),
                 SiCalib::N_ADC_BINS, SiCalib::ADC_MIN, SiCalib::ADC_MAX);

    Long64_t n_entries = tree->GetEntries();
    for (Long64_t j = 0; j < n_entries; j++) {
      tree->GetEntry(j);
      if (channel == SiCalib::SI_DETECTOR_CHANNEL) {
        hist->Fill(energy);
      }
    }

    Double_t mu_guess = SiCalib::MU_GUESSES.at(run);
    Double_t fit_min =
        TMath::Max(SiCalib::ADC_MIN, mu_guess - SiCalib::FIT_HALF_WIDTH_ADC);
    Double_t fit_max =
        TMath::Min(SiCalib::ADC_MAX, mu_guess + SiCalib::FIT_HALF_WIDTH_ADC);

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
    TF1 *fit_func = fitter.GetFitFunction();
    if (fit_func) {
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
