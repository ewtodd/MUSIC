#include "EventsSummary.hpp"

void CreateSummaryHistograms(SummaryHistograms &h,
                             const SummaryHistConfig &cfg) {
  TString musicTitle =
      "MUSIC strip energies (complete events);Strip;#DeltaE [" +
      cfg.unit_label + "]";
  h.h_music =
      new TH2F(PlottingUtils::GetRandomName().Data(), musicTitle, 18, -0.5,
               17.5, cfg.music_energy_bins, cfg.strip_e_min, cfg.strip_e_max);
  h.h_mult =
      new TH1F(PlottingUtils::GetRandomName().Data(),
               "Event multiplicity (complete events);Multiplicity;Counts", 36,
               -0.5, 35.5);

  for (Int_t s = 1; s <= 16; s++) {
    Double_t lMax = cfg.strip_e_max;
    Double_t rMax = cfg.strip_e_max;
    if (cfg.odd_even_split) {
      lMax = (s % 2 == 0) ? cfg.left_even_max : cfg.left_odd_max;
      rMax = (s % 2 == 0) ? cfg.right_even_max : cfg.right_odd_max;
    }
    TString rlTitle = Form(";Strip %d L #DeltaE [", s) + cfg.unit_label + "];" +
                      Form("Strip %d R #DeltaE [", s) + cfg.unit_label + "]";
    h.h2_R_vs_L[s] =
        new TH2F(PlottingUtils::GetRandomName().Data(), rlTitle, 200,
                 cfg.strip_e_min, lMax, 200, cfg.strip_e_min, rMax);
  }

  if (Constants::cfg.HAS_CATHODE)
    h.h1_cathode = new TH1F(PlottingUtils::GetRandomName().Data(),
                            ";Cathode #DeltaE [" + cfg.unit_label + "];Counts",
                            400, 0.0, cfg.cathode_max);
  if (Constants::cfg.HAS_STRIP17)
    h.h1_strip17 = new TH1F(PlottingUtils::GetRandomName().Data(),
                            ";Strip17 #DeltaE [" + cfg.unit_label + "];Counts",
                            400, cfg.strip_e_min, cfg.strip_e_max);

  if (Constants::cfg.HAS_STRIP0 && Constants::cfg.HAS_GRID) {
    h.h2_strip0_vs_grid =
        new TH2F(PlottingUtils::GetRandomName().Data(),
                 ";Grid #DeltaE [" + cfg.unit_label + "];Strip0 #DeltaE [" +
                     cfg.unit_label + "]",
                 200, 0.0, cfg.grid_max, 200, cfg.strip_e_min, cfg.strip0_max);
  } else if (Constants::cfg.HAS_STRIP0) {
    h.h1_strip0 = new TH1F(PlottingUtils::GetRandomName().Data(),
                           ";Strip0 #DeltaE [" + cfg.unit_label + "];Counts",
                           400, cfg.strip_e_min, cfg.strip0_max);
  } else if (Constants::cfg.HAS_GRID) {
    h.h1_grid = new TH1F(PlottingUtils::GetRandomName().Data(),
                         ";Grid #DeltaE [" + cfg.unit_label + "];Counts", 400,
                         0.0, cfg.grid_max);
  }
}

void SaveAndDeleteSummaryHistograms(SummaryHistograms &h, TFile *out_file,
                                    const TString &subdir,
                                    const TString &plot_suffix) {
  TString musicName = "music_strip_energies" + plot_suffix;
  TString multName = "multiplicity" + plot_suffix;
  TString cathName = "cathode" + plot_suffix;
  TString strip17Name = "strip17" + plot_suffix;
  TString s0gName = "strip0_vs_grid" + plot_suffix;
  TString s0Name = "strip0" + plot_suffix;
  TString gName = "grid" + plot_suffix;

  TCanvas *c_music = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c_music->cd();
  PlottingUtils::ConfigureAndDraw2DHistogram(h.h_music, c_music);
  h.h_music->GetYaxis()->SetTitleOffset(1.4);
  c_music->SetLeftMargin(0.18);
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_music, musicName, subdir,
                              PlotSaveOptions::kLINEAR);
  out_file->cd();
  c_music->Write(h.h_music->GetName(), TObject::kOverwrite);
  delete c_music;

  TCanvas *c_mult = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c_mult->cd();
  PlottingUtils::ConfigureAndDrawHistogram(h.h_mult, kBlue + 1);
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_mult, multName, subdir, PlotSaveOptions::kLOG);
  delete c_mult;

  for (Int_t s = 1; s <= 16; s++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDraw2DHistogram(h.h2_R_vs_L[s], c);
    h.h2_R_vs_L[s]->GetYaxis()->SetTitleOffset(1.3);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, TString("R_vs_L_s") + s + plot_suffix,
                                subdir, PlotSaveOptions::kLINEAR);
    out_file->cd();
    c->Write(h.h2_R_vs_L[s]->GetName(), TObject::kOverwrite);
    delete c;
  }

  if (h.h1_cathode) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDrawHistogram(h.h1_cathode, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, cathName, subdir, PlotSaveOptions::kLOG);
    out_file->cd();
    h.h1_cathode->Write(h.h1_cathode->GetName(), TObject::kOverwrite);
    delete c;
  }

  if (h.h1_strip17) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDrawHistogram(h.h1_strip17, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, strip17Name, subdir, PlotSaveOptions::kLOG);
    out_file->cd();
    h.h1_strip17->Write(h.h1_strip17->GetName(), TObject::kOverwrite);
    delete c;
  }

  if (h.h2_strip0_vs_grid) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDraw2DHistogram(h.h2_strip0_vs_grid, c);
    h.h2_strip0_vs_grid->GetYaxis()->SetTitleOffset(1.3);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, s0gName, subdir, PlotSaveOptions::kLINEAR);
    out_file->cd();
    c->Write(h.h2_strip0_vs_grid->GetName(), TObject::kOverwrite);
    delete c;
  }
  if (h.h1_strip0) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDrawHistogram(h.h1_strip0, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, s0Name, subdir, PlotSaveOptions::kLOG);
    out_file->cd();
    h.h1_strip0->Write(h.h1_strip0->GetName(), TObject::kOverwrite);
    delete c;
  }
  if (h.h1_grid) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->cd();
    PlottingUtils::ConfigureAndDrawHistogram(h.h1_grid, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c, gName, subdir, PlotSaveOptions::kLOG);
    out_file->cd();
    h.h1_grid->Write(h.h1_grid->GetName(), TObject::kOverwrite);
    delete c;
  }

  for (Int_t s = 1; s <= 16; s++)
    delete h.h2_R_vs_L[s];
  delete h.h_music;
  delete h.h_mult;
  delete h.h1_cathode;
  delete h.h1_strip17;
  delete h.h2_strip0_vs_grid;
  delete h.h1_strip0;
  delete h.h1_grid;
  h.h_music = nullptr;
  h.h_mult = nullptr;
  h.h1_cathode = nullptr;
  h.h1_strip17 = nullptr;
  h.h2_strip0_vs_grid = nullptr;
  h.h1_strip0 = nullptr;
  h.h1_grid = nullptr;
}

void EventsSummary::SaveSampleTraces(const std::vector<TGraph *> &traces,
                                     const TString &save_name,
                                     const TString &subdir, Double_t y_min,
                                     Double_t y_max, const char *y_title) {
  if (traces.empty())
    return;

  Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  TH2F *frame = new TH2F(PlottingUtils::GetRandomName().Data(),
                         Form(";Strip;%s", y_title), s_hi - s_lo + 1,
                         s_lo - 0.5, s_hi + 0.5, 100, y_min, y_max);
  frame->SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c->cd();
  frame->Draw();

  for (Int_t i = 0; i < Int_t(traces.size()); i++) {
    traces[i]->SetLineColor(kBlack);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }

  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  delete c;
  delete frame;
}

void EventsSummary::BuildNormedSummaryHistograms(const TString &input_filename,
                                                 const TString &file_label) {
  TString input_filepath = input_filename + ".root";

  TFile *input_file = IO::OpenForWriting(input_filepath, "UPDATE");
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "[" << file_label
              << "] cannot open UPDATE for normed summary: " << input_filepath
              << std::endl;
    if (input_file)
      delete input_file;
    return;
  }
  TTree *input_tree = static_cast<TTree *>(input_file->Get("events"));
  if (!input_tree) {
    std::cerr << "[" << file_label << "] no events tree for normed summary"
              << std::endl;
    input_file->Close();
    delete input_file;
    return;
  }
  EnergyView ev;
  ev.Attach(input_tree);
  if (!ev.is_normed) {
    std::cout
        << "[" << file_label
        << "] no calibration in events file; skipping normed summary build"
        << std::endl;
    input_file->Close();
    delete input_file;
    return;
  }

  const Double_t strip_e_min = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t strip_e_max = Constants::cfg.STRIP_DE_MAX_NORMED;

  SummaryHistConfig cfg;
  cfg.unit_label = "a.u.";
  cfg.strip_e_min = strip_e_min;
  cfg.strip_e_max = strip_e_max;
  cfg.odd_even_split = kFALSE;
  cfg.left_odd_max = strip_e_max;
  cfg.left_even_max = strip_e_max;
  cfg.right_odd_max = strip_e_max;
  cfg.right_even_max = strip_e_max;
  cfg.cathode_max = 1.0;
  cfg.strip17_max = strip_e_max;
  cfg.grid_max = 1.0;
  cfg.strip0_max = strip_e_max;
  cfg.music_energy_bins = 400;

  SummaryHistograms h;
  CreateSummaryHistograms(h, cfg);

  Long64_t n_entries = input_tree->GetEntries();
  std::cout << "[" << file_label << "] building normed summary over "
            << n_entries << " events..." << std::endl;

  // Sample traces for overlay plot
  std::vector<TGraph *> sample_traces;
  Long64_t sample_stride = 0;
  if (Constants::cfg.SAVE_SAMPLE_TRACES > 0) {
    sample_stride = n_entries / Long64_t(Constants::cfg.SAVE_SAMPLE_TRACES);
    if (sample_stride < 1)
      sample_stride = 1;
  }

  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    ev.Decode();

    for (Int_t s = 0; s < 18; s++)
      h.h_music->Fill(Double_t(s), ev.total[s]);

    for (Int_t s = 1; s <= 16; s++)
      h.h2_R_vs_L[s]->Fill(ev.left[s], ev.right[s]);

    if (h.h1_cathode && ev.cathode > 0.0)
      h.h1_cathode->Fill(ev.cathode);

    if (h.h1_strip17)
      h.h1_strip17->Fill(ev.total[17]);

    if (h.h2_strip0_vs_grid)
      h.h2_strip0_vs_grid->Fill(ev.grid, ev.total[0]);

    if (h.h1_strip0)
      h.h1_strip0->Fill(ev.total[0]);

    if (h.h1_grid)
      h.h1_grid->Fill(ev.grid);

    Int_t mult = 0;
    for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
      mult += ev.hits_adc[k];
    h.h_mult->Fill(Double_t(mult));

    // Collect sample traces
    if (sample_stride > 0 && j % sample_stride == 0 &&
        Int_t(sample_traces.size()) < Constants::cfg.SAVE_SAMPLE_TRACES) {
      sample_traces.push_back(EventsSummary::BuildTraceFromTotals(ev.total));
    }
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString subdir = "events_summary_normed/" + file_label;
    SaveAndDeleteSummaryHistograms(h, input_file, subdir, "_normed");

    // Sample traces overlay
    if (!sample_traces.empty()) {
      SaveSampleTraces(sample_traces, "sample_traces_normed", subdir,
                       strip_e_min, strip_e_max, "#DeltaE [a.u.]");
    }
  }

  for (Int_t i = 0; i < Int_t(sample_traces.size()); i++)
    delete sample_traces[i];

  input_file->Close();
  delete input_file;
}

TGraph *EventsSummary::BuildTraceFromTotals(const Double_t *total) {
  Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  Int_t n_pts = s_hi - s_lo + 1;
  TGraph *g = new TGraph(n_pts);
  for (Int_t k = 0; k < n_pts; k++)
    g->SetPoint(k, s_lo + k, total[s_lo + k]);
  return g;
}
