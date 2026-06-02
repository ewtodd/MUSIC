#include "Pipeline.hpp"

std::mutex fused_log_mutex;

Bool_t FusedExists(const TString &subpath) {
  TString full = IO::GetRootFilesBaseDir() + "/" + subpath;
  return !gSystem->AccessPathName(full);
}

Double_t FusedSecSince(const std::chrono::steady_clock::time_point &t0) {
  return std::chrono::duration<Double_t>(std::chrono::steady_clock::now() - t0)
      .count();
}

// /proc/self/statm reports VmRSS in pages (field 2); 4 KiB/page on Linux.
// Worker-local label so concurrent log lines stay attributable.
void PrintMemUsage(const char *label) {
  Long64_t rss = 0;
  std::ifstream statm("/proc/self/statm");
  Long64_t dummy;
  statm >> dummy >> rss;
  Double_t rss_gb = rss * 4096.0 / (1024.0 * 1024.0 * 1024.0);
  std::lock_guard<std::mutex> lock(fused_log_mutex);
  std::cout << "[MEM] " << label << ": " << rss_gb << " GB RSS" << std::endl;
}

Bool_t EnsureRunHeaderFused(Int_t run, UShort_t &header) {
  if (BinaryToRoot::ReadHeaderSidecar(run, header))
    return kTRUE;
  FileSpec s0;
  s0.run = run;
  s0.suffix = "";
  TString bin_path = FileSet::CompassBinPath(s0);
  std::pair<std::vector<RawHit>, UShort_t> p =
      InitUtils::ConvertCoMPASSBinToHits(bin_path, 0);
  if (p.second == 0)
    return kFALSE;
  header = p.second;
  BinaryToRoot::WriteHeaderSidecar(run, header);
  return kTRUE;
}

Bool_t RunFusedPipelineForFile(FileSpec spec, UShort_t run_header,
                               const EventBuilder::SlotMap &slot_map,
                               const std::vector<ChannelCal> &sim_chans) {
  TString file_label = FileSet::FileLabel(spec);
  std::chrono::steady_clock::time_point t_total =
      std::chrono::steady_clock::now();
  Double_t t_parse = 0, t_timing = 0, t_apply = 0, t_events = 0, t_cal = 0,
           t_traces = 0;

  if (Constants::SKIP_EXISTING &&
      FusedExists(FileSet::EventsName(spec) + ".root")) {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cout << "[skip] " << file_label << " events file exists" << std::endl;
    return kTRUE;
  }

  TString bin_path = FileSet::CompassBinPath(spec);
  if (gSystem->AccessPathName(bin_path)) {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cerr << "[fail] " << file_label << " BIN missing: " << bin_path
              << std::endl;
    return kFALSE;
  }

  UShort_t use_header = (spec.suffix == "") ? 0 : run_header;

  PrintMemUsage((TString("before binary read ") + file_label).Data());

  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  std::pair<std::vector<RawHit>, UShort_t> parsed =
      InitUtils::ConvertCoMPASSBinToHits(bin_path, use_header);
  t_parse = FusedSecSince(t0);

  PrintMemUsage((TString("after binary read ") + file_label).Data());

  std::vector<RawHit> &hits = parsed.first;
  if (hits.empty() || parsed.second == 0) {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cerr << "[fail] " << file_label << " parse produced no hits"
              << std::endl;
    return kFALSE;
  }
  if (spec.suffix == "")
    BinaryToRoot::WriteHeaderSidecar(spec.run, parsed.second);

  t0 = std::chrono::steady_clock::now();
  TimeShiftResult shift_result = Timing::CalcTimeShiftsBeamMethodFromHits(
      hits, file_label, Constants::TIMING_REF_BOARD,
      Constants::TIMING_REF_BOARD_CHANNELS, Constants::TIMING_MIN_ENERGY,
      Constants::TIMING_MAX_ENERGY, Constants::TIMING_OVERLAP_MARGIN_S,
      Constants::TIMING_THRESH_DT_US);
  t_timing = FusedSecSince(t0);

  PrintMemUsage((TString("after timing ") + file_label).Data());

  t0 = std::chrono::steady_clock::now();
  Timing::ApplyShiftsInPlace(hits, shift_result.board_shifts);
  Timing::SortHitsByTimestamp(hits);
  t_apply = FusedSecSince(t0);

  PrintMemUsage((TString("after apply+sort ") + file_label).Data());

  t0 = std::chrono::steady_clock::now();
  Bool_t build_ok = EventBuilder::BuildEventsFromSortedHits(
      hits, slot_map, FileSet::EventsName(spec), file_label);
  t_events = FusedSecSince(t0);

  PrintMemUsage((TString("after event build ") + file_label).Data());

  std::vector<RawHit>().swap(hits);

  PrintMemUsage((TString("after hits free ") + file_label).Data());

  if (!build_ok) {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cerr << "[fail] " << file_label << " event build failed" << std::endl;
    return kFALSE;
  }

  if (!sim_chans.empty()) {
    t0 = std::chrono::steady_clock::now();
    CalibrateBeam::CalibrateBeamOneSubfile(spec, sim_chans);
    t_cal = FusedSecSince(t0);
    PrintMemUsage((TString("after calibration ") + file_label).Data());
  }

  std::vector<TString> events_name_vec = {FileSet::EventsName(spec)};
  std::vector<TString> file_labels_vec = {file_label};

  if (Constants::RUN_TRACES) {
    t0 = std::chrono::steady_clock::now();
    TraceCreator::BuildTraces(events_name_vec, file_labels_vec, kTRUE, kTRUE);
    t_traces = FusedSecSince(t0);
    PrintMemUsage((TString("after traces ") + file_label).Data());
  }

  Double_t total = FusedSecSince(t_total);
  {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cout << std::fixed << std::setprecision(1) << "[done] " << file_label
              << " total=" << total << "s  parse=" << t_parse
              << "  timing=" << t_timing << "  apply=" << t_apply
              << "  events=" << t_events << "  cal=" << t_cal
              << "  traces=" << t_traces << std::endl;
  }
  return kTRUE;
}

void Pipeline::Run() {
  ROOT::EnableThreadSafety();
  GpuAccel::Init();
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TString log_path = project_root + "/pipeline_fused.log";
  std::ofstream log_file(log_path.Data());
  std::streambuf *saved_cout = std::cout.rdbuf(log_file.rdbuf());
  std::streambuf *saved_cerr = std::cerr.rdbuf(log_file.rdbuf());
  Int_t saved_error_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  std::vector<FileSpec> specs = FileSet::BuildFileSpecs();
  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> unique_runs;
  for (Int_t k = 0; k < n_specs; k++)
    unique_runs.insert(specs[k].run);

  std::cout << "Phase A: gathering global headers for " << unique_runs.size()
            << " run(s)..." << std::endl;
  std::map<Int_t, UShort_t> run_headers;
  for (std::set<Int_t>::const_iterator it = unique_runs.begin();
       it != unique_runs.end(); ++it) {
    UShort_t h;
    if (!EnsureRunHeaderFused(*it, h)) {
      std::cerr << "Header gather FAILED for run " << *it << std::endl;
      continue;
    }
    run_headers[*it] = h;
    std::cout << "  Run " << *it << " header 0x" << std::hex << h << std::dec
              << std::endl;
  }

  EventBuilder::SlotMap slot_map = EventBuilder::BuildSlotMap();

  std::vector<ChannelCal> sim_chans;
  if (Constants::SKIP_CALIBRATION) {
    std::cout << "SKIP_CALIBRATION=true; skipping beam calibration and eres "
                 "aggregation (events kept in raw ADC)."
              << std::endl;
  } else {
    TString sim_beam_path = CalibrateBeam::DefaultSimBeamPath(project_root);
    sim_chans = CalibrateBeam::LoadSimChans(sim_beam_path);
    if (sim_chans.empty())
      std::cerr << "Sim anchors unavailable; per-subfile calibration disabled"
                << std::endl;
  }

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::MAX_FUSED_WORKERS);
  std::cout << "Phase B: fused pipeline on " << n_specs << " files with "
            << n_workers << " workers." << std::endl;

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
        FileSpec spec = specs[k];
        UShort_t header =
            run_headers.count(spec.run) ? run_headers[spec.run] : UShort_t(0);
        Bool_t ok = RunFusedPipelineForFile(spec, header, slot_map, sim_chans);
        if (!ok) {
          std::lock_guard<std::mutex> lk(fused_log_mutex);
          std::cerr << "FAILED: " << FileSet::FileLabel(spec) << std::endl;
        }
      }
    });
  }
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();

  std::cout << "All fused pipelines complete." << std::endl;

  if (!sim_chans.empty()) {
    std::cout << "Phase C: per-run eres TOML aggregation" << std::endl;
    for (std::set<Int_t>::const_iterator it = unique_runs.begin();
         it != unique_runs.end(); ++it) {
      std::vector<FileSpec> run_specs;
      for (Int_t k = 0; k < n_specs; k++)
        if (specs[k].run == *it)
          run_specs.push_back(specs[k]);
      CalibrateBeam::AggregateEresTomlForRun(*it, run_specs);
    }
  }

  std::cout.rdbuf(saved_cout);
  std::cerr.rdbuf(saved_cerr);
  gErrorIgnoreLevel = saved_error_level;
  log_file.close();
  std::cout << "Fused pipeline finished. Output logged to " << log_path
            << std::endl;
}
