#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstddef>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <ostream>
#include <thread>
#include <vector>

static std::mutex bintoroot_print_mutex;

Bool_t ConvertOneFile(TString filepath, TString output_name,
                      UShort_t global_header) {
  UShort_t returned =
      InitUtils::ConvertCoMPASSBinToROOT(filepath, output_name, global_header);
  std::lock_guard<std::mutex> lock(bintoroot_print_mutex);
  std::cout << "Converted: " << output_name << " (header 0x" << std::hex
            << returned << std::dec << ")" << std::endl;
  return kTRUE;
}

void BinaryToRoot() {
  ROOT::EnableThreadSafety();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");

  TString base_path = "/home/e-work/LabData/MUSIC/37Cl/";

  struct Job {
    TString filepath;
    TString output_name;
    Int_t run;
  };

  std::vector<Job> file0_jobs;
  std::vector<Job> remaining_jobs;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (std::size_t k = 0; k < specs.size(); k++) {
    const FileSpec &s = specs[k];
    TString raw_dir = base_path + Form("run_%d/RAW/", s.run);
    Job j;
    j.filepath = raw_dir + Form("DataR_run_%d%s.BIN", s.run, s.suffix.Data());
    j.output_name = Form("DataR_run_%d%s", s.run, s.suffix.Data());
    j.run = s.run;
    if (s.suffix == "")
      file0_jobs.push_back(j);
    else
      remaining_jobs.push_back(j);
    std::cout << "Will process: " << j.filepath << std::endl;
  }

  std::map<Int_t, UShort_t> run_headers;
  for (std::size_t k = 0; k < file0_jobs.size(); k++) {
    const Job &j = file0_jobs[k];
    UShort_t header =
        InitUtils::ConvertCoMPASSBinToROOT(j.filepath, j.output_name, 0);
    run_headers[j.run] = header;
    std::cout << "Run " << j.run << " saved global header: 0x" << std::hex
              << header << std::dec << std::endl;
  }

  if (remaining_jobs.empty()) {
    std::cout << "All conversions complete." << std::endl;
    return;
  }

  Int_t n_remaining = Int_t(remaining_jobs.size());
  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_remaining);

  std::cout << "Converting " << n_remaining << " remaining files with "
            << n_workers << " workers." << std::endl;

  for (Int_t i = 0; i < n_remaining; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_remaining);

    for (Int_t k = i; k < batch_end; ++k) {
      const Job &j = remaining_jobs[k];
      futures.push_back(std::async(std::launch::async, ConvertOneFile,
                                   j.filepath, j.output_name,
                                   run_headers[j.run]));
    }

    for (size_t k = 0; k < futures.size(); ++k) {
      Bool_t result = futures[k].get();
      if (!result) {
        std::cerr << "FAILED: " << remaining_jobs[i + k].output_name
                  << std::endl;
      }
    }
  }

  std::cout << "All conversions complete." << std::endl;
}
