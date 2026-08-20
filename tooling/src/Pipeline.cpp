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

  if (Constants::cfg.USE_SOLARIS_DATA) {
    // SOLARIS: find first available file (chunk or original) for header
    std::vector<TString> suffixes = FileSet::DiscoverSolRunSuffixes(run);
    Bool_t found = kFALSE;
    for (Int_t k = 0; k < Int_t(suffixes.size()); k++) {
      FileSpec s0;
      s0.run = run;
      s0.suffix = suffixes[k];
      TString sol_path = FileSet::SolBinPath(s0);
      if (gSystem->AccessPathName(sol_path))
        continue;

      SOLReader reader;
      if (!reader.Open(sol_path.Data()))
        continue;
      if (!reader.ReadEvent()) {
        reader.Close();
        continue;
      }
      header = reader.GetCurrentEvent().block_header;
      reader.Close();
      found = kTRUE;
      break;
    }
    if (!found) {
      std::cerr << "SOL header gather FAILED for run " << run
                << " (no accessible files)" << std::endl;
      return kFALSE;
    }
  } else {
    // CoMPASS: read global header from first .BIN file
    FileSpec s0;
    s0.run = run;
    s0.suffix = "";
    TString bin_path = FileSet::CompassBinPath(s0);
    std::pair<std::vector<RawHit>, UShort_t> p =
        InitUtils::ConvertCoMPASSBinToHits(bin_path, 0);
    if (p.second == 0)
      return kFALSE;
    header = p.second;
  }

  BinaryToRoot::WriteHeaderSidecar(run, header);
  return kTRUE;
}

Bool_t RunFusedPipelineForFile(FileSpec spec, UShort_t run_header,
                               const EventBuilder::SlotMap &slot_map,
                               const std::vector<ChannelCal> &chans) {
  TString file_label = FileSet::FileLabel(spec);
  std::chrono::steady_clock::time_point t_total =
      std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point t0;
  Double_t t_parse = 0, t_timing = 0, t_apply = 0, t_events = 0, t_cal = 0;

  // SKIP_EXISTING skips the expensive processing when the events file already
  // exists -- but the plots below are still (re)made from that file.
  const Bool_t skip_processing =
      Constants::cfg.SKIP_EXISTING &&
      FusedExists(FileSet::EventsName(spec) + ".root");

  if (skip_processing) {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cout << "[skip-build] " << file_label
              << " events exist; re-making plots only" << std::endl;
  } else {
    TString bin_path;
    if (Constants::cfg.USE_SOLARIS_DATA) {
      bin_path = FileSet::SolBinPath(spec);
    } else {
      bin_path = FileSet::CompassBinPath(spec);
    }

    if (gSystem->AccessPathName(bin_path)) {
      std::lock_guard<std::mutex> lock(fused_log_mutex);
      std::cerr << "[fail] " << file_label
                << (Constants::cfg.USE_SOLARIS_DATA ? " SOL" : " BIN")
                << " missing: " << bin_path << std::endl;
      return kFALSE;
    }

    PrintMemUsage((TString("before binary read ") + file_label).Data());

    t0 = std::chrono::steady_clock::now();
    std::vector<RawHit> hits;

    if (Constants::cfg.USE_SOLARIS_DATA) {
      // SOLARIS: stream blocks directly to RawHit (no intermediate SOLHit
      // vector)
      SOLReader sol_reader;
      sol_reader.SetSkipTraces(kTRUE);
      if (!sol_reader.Open(bin_path.Data())) {
        std::lock_guard<std::mutex> lock(fused_log_mutex);
        std::cerr << "[fail] " << file_label << " cannot open SOL file"
                  << std::endl;
        return kFALSE;
      }
      // Reserve an estimate to avoid repeated reallocation (~1M blocks typical)
      hits.reserve(1048576);
      while (sol_reader.ReadEvent()) {
        const SOLData &sol = sol_reader.GetCurrentEvent();
        RawHit raw;
        raw.board = 0;
        raw.channel = sol.channel;
        raw.energy = sol.energy;
        raw.timestamp = sol.timestamp * 1000;
        raw.flags = MapSOLFlagsToCoMPASS(sol.flags_high, sol.flags_low);
        hits.push_back(raw);
      }
      sol_reader.Close();
    } else {
      // CoMPASS: direct conversion to RawHit
      UShort_t use_header = (spec.suffix == "") ? 0 : run_header;
      std::pair<std::vector<RawHit>, UShort_t> parsed =
          InitUtils::ConvertCoMPASSBinToHits(bin_path, use_header);
      hits = parsed.first;
      if (spec.suffix == "")
        BinaryToRoot::WriteHeaderSidecar(spec.run, parsed.second);
    }
    t_parse = FusedSecSince(t0);

    PrintMemUsage((TString("after binary read ") + file_label).Data());

    if (hits.empty()) {
      std::lock_guard<std::mutex> lock(fused_log_mutex);
      std::cerr << "[fail] " << file_label << " parse produced no hits"
                << std::endl;
      return kFALSE;
    }

    if (Constants::cfg.TIMING_DO_BOARD_SYNC || Constants::cfg.TIMING_DO_SORT) {
      t0 = std::chrono::steady_clock::now();
      TimeShiftResult shift_result = Timing::CalcTimeShiftsBeamMethodFromHits(
          hits, file_label, Constants::cfg.TIMING_REF_BOARD,
          Constants::cfg.TIMING_REF_BOARD_CHANNELS,
          Constants::cfg.TIMING_MIN_ENERGY, Constants::cfg.TIMING_MAX_ENERGY,
          Constants::cfg.TIMING_OVERLAP_MARGIN_S,
          Constants::cfg.TIMING_THRESH_DT_US);
      t_timing = FusedSecSince(t0);

      PrintMemUsage((TString("after timing ") + file_label).Data());

      t0 = std::chrono::steady_clock::now();
      Timing::ApplyShiftsInPlace(hits, shift_result.board_shifts);
      if (Constants::cfg.TIMING_DO_SORT)
        Timing::SortHitsByTimestamp(hits);
      t_apply = FusedSecSince(t0);

      PrintMemUsage((TString("after apply+sort ") + file_label).Data());
    } else {
      std::lock_guard<std::mutex> lock(fused_log_mutex);
      std::cout << "[skip-timing] " << file_label
                << " board sync and sort both disabled; skipping to event "
                   "build"
                << std::endl;
    }

    t0 = std::chrono::steady_clock::now();
    Bool_t build_ok = EventBuilder::BuildEventsFromSortedHits(
        hits, slot_map, FileSet::EventsName(spec), file_label);
    t_events = FusedSecSince(t0);

    PrintMemUsage((TString("after event build ") + file_label).Data());

    std::vector<RawHit>().swap(hits);

    PrintMemUsage((TString("after hits free ") + file_label).Data());

    if (!build_ok) {
      std::lock_guard<std::mutex> lock(fused_log_mutex);
      std::cerr << "[fail] " << file_label << " event build failed"
                << std::endl;
      return kFALSE;
    }
  }

  // Calibration fits beam peaks from the (fresh or pre-existing) events file,
  // so it runs in plot-only mode too, regenerating its diagnostic hists.
  if (!chans.empty()) {
    t0 = std::chrono::steady_clock::now();
    CalibrateBeam::CalibrateBeamOneSubfile(spec, chans);
    t_cal = FusedSecSince(t0);
    PrintMemUsage((TString("after calibration ") + file_label).Data());
  }

  // Normed summary histograms (requires calibration from above).
  {
    t0 = std::chrono::steady_clock::now();
    EventsSummary::BuildNormedSummaryHistograms(FileSet::EventsName(spec),
                                                file_label);
    Double_t t_normed = FusedSecSince(t0);
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cout << "  normed summary: " << t_normed << "s" << std::endl;
  }

  Double_t total = FusedSecSince(t_total);
  {
    std::lock_guard<std::mutex> lock(fused_log_mutex);
    std::cout << std::fixed << std::setprecision(1) << "[done] " << file_label
              << " total=" << total << "s  parse=" << t_parse
              << "  timing=" << t_timing << "  apply=" << t_apply
              << "  events=" << t_events << "  cal=" << t_cal << std::endl;
  }
  return kTRUE;
}

void Pipeline::Run() {
  ROOT::EnableThreadSafety();
  GpuAccel::Init();
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");

  TString log_path = project_root + "/pipeline_fused.log";
  std::ofstream log_file(log_path.Data());
  std::streambuf *saved_cout = std::cout.rdbuf(log_file.rdbuf());
  std::streambuf *saved_cerr = std::cerr.rdbuf(log_file.rdbuf());
  Int_t saved_error_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  std::vector<FileSpec> specs = FileSet::BuildRawOrProcessedFileSpecs();
  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> unique_runs;
  for (Int_t k = 0; k < n_specs; k++)
    unique_runs.insert(specs[k].run);

  // The global header is consumed only when a subfile is (re)built from its
  // BIN; a fully-processed run never touches the raw dir, so skip it.
  std::set<Int_t> runs_needing_header;
  for (Int_t k = 0; k < n_specs; k++) {
    Bool_t will_build = !(Constants::cfg.SKIP_EXISTING &&
                          FusedExists(FileSet::EventsName(specs[k]) + ".root"));
    if (will_build)
      runs_needing_header.insert(specs[k].run);
  }

  std::cout << "Phase A: gathering global headers for "
            << runs_needing_header.size() << " run(s)..." << std::endl;
  std::map<Int_t, UShort_t> run_headers;
  for (std::set<Int_t>::const_iterator it = runs_needing_header.begin();
       it != runs_needing_header.end(); ++it) {
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

  std::vector<ChannelCal> chans;
  if (Constants::cfg.SKIP_CALIBRATION) {
    std::cout << "SKIP_CALIBRATION=true; skipping beam calibration and eres "
                 "aggregation (events kept in raw ADC)."
              << std::endl;
  } else {
    chans = CalibrateBeam::BuildChannels();
  }

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::cfg.MAX_FUSED_WORKERS);
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
        Bool_t ok = RunFusedPipelineForFile(spec, header, slot_map, chans);
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

  if (!chans.empty()) {
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
