#include "TraceCreator.hpp"

void TraceCreator::WriteH2CanvasToFile(TFile *file, TH2F *h,
                                       const TString &save_name,
                                       const TString &subdir) {
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  h->GetYaxis()->SetTitleOffset(1.3);
  c->SetLeftMargin(0.18);
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  file->cd();
  c->Write(h->GetName(), TObject::kOverwrite);
  delete c;
}

void TraceCreator::BuildMeVSummaryHistograms(const TString &input_filename,
                                             const TString &file_label,
                                             const FileSpec &spec) {
  TString input_filepath = input_filename + ".root";

  TFile *input_file = IO::OpenForWriting(input_filepath, "UPDATE");
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "[" << file_label
              << "] cannot open UPDATE for MeV summary: " << input_filepath
              << std::endl;
    if (input_file)
      delete input_file;
    return;
  }
  TTree *input_tree = static_cast<TTree *>(input_file->Get("events"));
  if (!input_tree) {
    std::cerr << "[" << file_label << "] no events tree for MeV summary"
              << std::endl;
    input_file->Close();
    delete input_file;
    return;
  }
  EnergyView ev;
  ev.Attach(input_tree);
  if (!ev.is_mev) {
    std::cout << "[" << file_label
              << "] no calibration in events file; skipping MeV summary build"
              << std::endl;
    input_file->Close();
    delete input_file;
    return;
  }

  UInt_t flags_or = 0;
  input_tree->SetBranchAddress("FlagsOR", &flags_or);

  const Double_t strip_e_min = Constants::STRIP_E_MIN_MEV;
  const Double_t strip_e_max = Constants::STRIP_E_MAX_MEV;
  const Double_t total_e_min = Constants::TOTAL_E_MIN_MEV;
  const Double_t total_e_max = Constants::TOTAL_E_MAX_MEV;
  const Double_t cathode_e_max = Constants::CATHODE_E_MAX_MEV;

  TH2F *h_music = new TH2F(
      "hMUSIC_MeV", "MUSIC strip energies (complete events);Strip;Energy [MeV]",
      18, -0.5, 17.5, 400, strip_e_min, strip_e_max);
  TH2F *h_music_clean =
      new TH2F("hMUSICClean_MeV",
               "MUSIC strip energies (complete, no flags);Strip;Energy [MeV]",
               18, -0.5, 17.5, 400, strip_e_min, strip_e_max);
  TH2F *h_music_flagged =
      new TH2F("hMUSICFlagged_MeV",
               "MUSIC strip energies (complete, flagged);Strip;Energy [MeV]",
               18, -0.5, 17.5, 400, strip_e_min, strip_e_max);
  TH2F *h2_strip[18];
  for (Int_t s = 0; s < 18; s++) {
    h2_strip[s] =
        new TH2F(Form("h2_totalE_vs_stripE_s%d_MeV", s),
                 Form(";Strip %d #DeltaE [MeV];Total #DeltaE [MeV]", s), 200,
                 strip_e_min, strip_e_max, 400, total_e_min, total_e_max);
  }
  TH2F *h2_cath = new TH2F("h2_totalE_vs_cathode_MeV",
                           ";Cathode #DeltaE [MeV];Total #DeltaE [MeV]", 200,
                           0.0, cathode_e_max, 400, total_e_min, total_e_max);

  Long64_t n_entries = input_tree->GetEntries();
  std::cout << "[" << file_label << "] rebuilding MeV summary over "
            << n_entries << " events..." << std::endl;
  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    ev.Decode();
    Double_t event_total = 0.0;
    for (Int_t s = 0; s < 18; s++)
      event_total += ev.total[s];
    Bool_t has_any_flag = (flags_or != 0);
    for (Int_t s = 0; s < 18; s++) {
      h_music->Fill(Double_t(s), ev.total[s]);
      if (has_any_flag)
        h_music_flagged->Fill(Double_t(s), ev.total[s]);
      else
        h_music_clean->Fill(Double_t(s), ev.total[s]);
      h2_strip[s]->Fill(ev.total[s], event_total);
    }
    if (ev.cathode_adc != -1)
      h2_cath->Fill(ev.cathode, event_total);
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString music_subdir = "events_nearest_mev/" + file_label;
    TString trace_subdir = "trace_summary_mev/" + file_label;
    WriteH2CanvasToFile(input_file, h_music, "music_strip_energies_mev",
                        music_subdir);
    WriteH2CanvasToFile(input_file, h_music_clean,
                        "music_strip_energies_clean_mev", music_subdir);
    WriteH2CanvasToFile(input_file, h_music_flagged,
                        "music_strip_energies_flagged_mev", music_subdir);
    for (Int_t s = 0; s < 18; s++)
      WriteH2CanvasToFile(input_file, h2_strip[s],
                          Form("totalE_vs_stripE_s%d_mev", s), trace_subdir);
    WriteH2CanvasToFile(input_file, h2_cath, "totalE_vs_cathode_mev",
                        trace_subdir);
  }

  for (Int_t s = 0; s < 18; s++)
    delete h2_strip[s];
  delete h_music;
  delete h_music_clean;
  delete h_music_flagged;
  delete h2_cath;

  input_file->Close();
  delete input_file;
}

TGraph *TraceCreator::BuildTraceFromTotals(const Double_t *total) {
  Int_t s_lo = Constants::IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::IGNORE_STRIP_17 ? 16 : 17;
  Int_t n_pts = s_hi - s_lo + 1;
  TGraph *g = new TGraph(n_pts);
  for (Int_t k = 0; k < n_pts; k++)
    g->SetPoint(k, s_lo + k, total[s_lo + k]);
  return g;
}

TGraph *TraceCreator::BuildEventTrace(const EnergyView &ev) {
  return BuildTraceFromTotals(ev.total);
}

void TraceCreator::BuildTraces(std::vector<TString> input_filenames,
                               std::vector<TString> file_labels,
                               Bool_t save_plots, Bool_t reprocess) {
  if (!reprocess || !save_plots)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = input_filename + ".root";
    TString file_label = file_labels[i];

    TFile *input_file = IO::OpenForReading(input_filepath);
    if (!input_file || input_file->IsZombie()) {
      std::cerr << "Error opening file: " << input_filepath << std::endl;
      continue;
    }

    TTree *input_tree = static_cast<TTree *>(input_file->Get("events"));
    if (!input_tree) {
      std::cerr << "No events tree in: " << input_filepath << std::endl;
      input_file->Close();
      continue;
    }

    FileSpec spec = FileSet::ResolveFileSpec(file_label);

    EnergyView ev;
    ev.Attach(input_tree);
    if (!ev.is_mev)
      std::cerr << "[" << file_label
                << "] WARNING: no MeV branches; using uncalibrated ADC values"
                << std::endl;
    const char *unit = ev.Unit();
    const Double_t strip_e_min =
        ev.is_mev ? Constants::STRIP_E_MIN_MEV : Constants::STRIP_E_MIN_ADC;
    const Double_t strip_e_max =
        ev.is_mev ? Constants::STRIP_E_MAX_MEV : Constants::STRIP_E_MAX_ADC;

    Long64_t n_entries = input_tree->GetEntries();
    Long64_t n_to_save =
        TMath::Min(Long64_t(Constants::MAX_TRACE_SAVES), n_entries);
    std::cout << "Saving " << n_to_save << " per-event traces from "
              << input_filename << "..." << std::endl;

    for (Long64_t j = 0; j < n_to_save; j++) {
      input_tree->GetEntry(j);
      ev.Decode();

      Int_t s_lo = Constants::IGNORE_STRIP_0 ? 1 : 0;
      Int_t s_hi = Constants::IGNORE_STRIP_17 ? 16 : 17;
      Int_t n_pts = s_hi - s_lo + 1;
      TGraph *TraceTotal = BuildEventTrace(ev);
      TGraph *TraceLeft = new TGraph(n_pts);
      TGraph *TraceRight = new TGraph(n_pts);
      for (Int_t k = 0; k < n_pts; k++) {
        Int_t s = s_lo + k;
        Double_t l = (s == 0 || s == 17) ? 0.0 : ev.left[s];
        Double_t r = (s == 0 || s == 17) ? 0.0 : ev.right[s];
        TraceLeft->SetPoint(k, s, l);
        TraceRight->SetPoint(k, s, r);
      }
      TString total_title =
          Form("Event %lld;Strip Index;#DeltaE [%s]", j, unit);
      PlottingUtils::ConfigureGraph(TraceTotal, kBlack, total_title);
      TraceTotal->SetMarkerColor(kBlack);
      TraceTotal->GetYaxis()->SetRangeUser(strip_e_min, strip_e_max);

      TraceLeft->SetLineColor(kBlue + 1);
      TraceLeft->SetMarkerColor(kBlue + 1);
      TraceLeft->SetLineWidth(2);

      TraceRight->SetLineColor(kRed + 1);
      TraceRight->SetMarkerColor(kRed + 1);
      TraceRight->SetLineWidth(2);

      {
        std::lock_guard<std::mutex> lock(g_plot_mutex);

        TString trace_subdir = "traces/" + file_label;
        TString trace_tag = Form("event%lld", j);

        TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
        TraceTotal->Draw("ALP");
        TraceLeft->Draw("LP SAME");
        TraceRight->Draw("LP SAME");
        if (Constants::SAVE_PLOTS)
          PlottingUtils::SaveFigure(c, "trace_" + trace_tag, trace_subdir,
                                    PlotSaveOptions::kLINEAR);
        delete c;

        std::cout << "Saved event " << j << " under " << trace_subdir
                  << std::endl;
      }

      delete TraceTotal;
      delete TraceLeft;
      delete TraceRight;
    }

    Bool_t had_cal = ev.is_mev;
    input_file->Close();
    delete input_file;
    if (had_cal)
      BuildMeVSummaryHistograms(input_filename, file_label, spec);
    std::cout << "Finished processing " << input_filename << std::endl;
  }
}
