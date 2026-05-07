#include "Constants.hpp"
#include "EventBuilder.cpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include "Timing.cpp"
#include "TraceCreator.cpp"
#include <TError.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <streambuf>
#include <thread>
#include <vector>

static std::mutex pipeline_log_mutex;

Bool_t RunPipelineOneFile(FileSpec spec) {
  std::vector<TString> raw_filepath = {RawRootName(spec)};
  std::vector<TString> sorted_name = {SortedName(spec)};
  std::vector<TString> events_name = {EventsName(spec)};
  std::vector<TString> file_labels = {FileLabel(spec)};

  std::vector<TimeShiftResult> results = CalcTimeShiftsBeamMethod(
      raw_filepath, file_labels, Constants::REF_BOARD,
      Constants::BOARD_CHANNELS[Constants::REF_BOARD], Constants::BOARD_CHANNELS,
      Constants::MIN_ENERGY, Constants::MAX_ENERGY, Constants::OVERLAP_MARGIN_S,
      Constants::THRESH_DT_US, kTRUE);
  ApplyTimeShift(raw_filepath, results, kTRUE);
  TimesortData(raw_filepath, sorted_name, kTRUE);

  BuildEvents(sorted_name, events_name, kTRUE);

  BuildTraces(events_name, file_labels, kTRUE);

  {
    std::lock_guard<std::mutex> lock(pipeline_log_mutex);
    std::cout << "Pipeline complete for " << file_labels[0] << std::endl;
  }
  return kTRUE;
}

void PipelineRunner() {
  ROOT::EnableThreadSafety();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");

  std::ofstream null_sink;
  std::streambuf *saved_cout = nullptr;
  std::streambuf *saved_cerr = nullptr;
  Int_t saved_error_level = gErrorIgnoreLevel;

  if (Constants::SILENT_PIPELINE) {
    null_sink.open("/dev/null");
    saved_cout = std::cout.rdbuf(null_sink.rdbuf());
    saved_cerr = std::cerr.rdbuf(null_sink.rdbuf());
    gErrorIgnoreLevel = kFatal;
  }

  std::vector<FileSpec> specs = BuildFileSpecs();
  Int_t n_specs = Int_t(specs.size());
  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);

  std::cout << "Running pipeline on " << n_specs << " files with " << n_workers
            << " workers." << std::endl;

  for (Int_t i = 0; i < n_specs; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_specs);

    for (Int_t k = i; k < batch_end; ++k) {
      futures.push_back(
          std::async(std::launch::async, RunPipelineOneFile, specs[k]));
    }

    for (size_t k = 0; k < futures.size(); ++k) {
      Bool_t result = futures[k].get();
      if (!result) {
        std::cerr << "FAILED: " << FileLabel(specs[i + k]) << std::endl;
      }
    }
  }

  std::cout << "All pipelines complete." << std::endl;

  if (Constants::SILENT_PIPELINE) {
    std::cout.rdbuf(saved_cout);
    std::cerr.rdbuf(saved_cerr);
    gErrorIgnoreLevel = saved_error_level;
    null_sink.close();
    std::cout << "Pipeline finished (output was silenced)." << std::endl;
  }
}
