#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "IOUtils.hpp"
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
#include <sstream>
#include <thread>
#include <vector>

static std::mutex bintoroot_print_mutex;

struct FlagSpec {
  UInt_t mask;
  const char *name;
};

static const FlagSpec kFlagSpecs[] = {
    {CoMPASSData::DEADTIME_OCCURRED, "Deadtime"},
    {CoMPASSData::TIMESTAMP_ROLLOVER, "Timestamp Rollover"},
    {CoMPASSData::TIMESTAMP_RESET_EXT, "Timestamp Reset (Ext)"},
    {CoMPASSData::FAKE_EVENT, "Fake Event"},
    {CoMPASSData::MEMORY_FULL, "Memory Full"},
    {CoMPASSData::TRIGGER_LOST, "Trigger Lost"},
    {CoMPASSData::N_TRIGGERS_LOST, "N Triggers Lost"},
    {CoMPASSData::SATURATION_IN_GATE, "Saturation In Gate"},
    {CoMPASSData::TRIGGERS_1024_COUNTED, "1024 Triggers"},
    {CoMPASSData::FIRST_AFTER_BUSY, "First After Busy"},
    {CoMPASSData::INPUT_SATURATING, "Input Saturating"},
    {CoMPASSData::N_TRIGGERS_COUNTED, "N Triggers Counted"},
    {CoMPASSData::NOT_MATCHED_TIMEFILTER, "Not Matched Time Filter"},
    {CoMPASSData::FINE_TIMESTAMP, "Fine Timestamp"},
    {CoMPASSData::PILEUP, "Pile-up"},
    {CoMPASSData::PLL_LOCK_LOSS, "PLL Lock Loss"},
    {CoMPASSData::OVER_TEMPERATURE, "Over Temperature"},
    {CoMPASSData::ADC_SHUTDOWN, "ADC Shutdown"}};

void TallyFlags(TString output_name, std::ostream &os) {
  TString filepath = output_name + ".root";
  TFile *file = IO::OpenForReading(filepath);
  if (!file || file->IsZombie()) {
    os << "  (could not open " << filepath << " to tally flags)" << std::endl;
    return;
  }
  TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
  if (!tree) {
    os << "  (no Data_R tree in " << filepath << ")" << std::endl;
    file->Close();
    return;
  }

  UInt_t flags;
  tree->SetBranchAddress("Flags", &flags);

  Long64_t n_entries = tree->GetEntries();
  Int_t n_specs = sizeof(kFlagSpecs) / sizeof(kFlagSpecs[0]);
  std::vector<Long64_t> counts(n_specs, 0);
  Long64_t any_flag = 0;

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    if (flags != 0)
      any_flag++;
    for (Int_t k = 0; k < n_specs; k++) {
      if (flags & kFlagSpecs[k].mask)
        counts[k]++;
    }
  }

  os << "  Flag tally (" << n_entries << " hits, " << any_flag
     << " with any flag set):" << std::endl;
  for (Int_t k = 0; k < n_specs; k++) {
    if (counts[k] == 0)
      continue;
    Double_t pct = n_entries > 0 ? 100.0 * counts[k] / n_entries : 0.0;
    os << "    " << kFlagSpecs[k].name << ": " << counts[k] << " (" << pct
       << "%)" << std::endl;
  }
  file->Close();
}

Bool_t ConvertOneFile(TString filepath, TString output_name,
                      UShort_t global_header) {
  UShort_t returned =
      InitUtils::ConvertCoMPASSBinToROOT(filepath, output_name, global_header);
  std::ostringstream summary;
  summary << "Converted: " << output_name << " (header 0x" << std::hex
          << returned << std::dec << ")" << std::endl;
  TallyFlags(output_name, summary);
  std::lock_guard<std::mutex> lock(bintoroot_print_mutex);
  std::cout << summary.str();
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
    TallyFlags(j.output_name, std::cout);
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
