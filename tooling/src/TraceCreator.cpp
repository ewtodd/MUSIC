#include "TraceCreator.hpp"

void TraceCreator::BuildNormedSummaryHistograms(const TString &input_filename,
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
  const Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  const Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;

  // Same histogram set as the raw ADC event builder, but in a.u.
  TH2F *h_music =
      new TH2F("hMUSIC_Normed",
               "MUSIC strip energies (complete events);Strip;#DeltaE [a.u.]",
               18, -0.5, 17.5, 400, strip_e_min, strip_e_max);
  TH1F *h_mult = new TH1F("hMult_Normed",
                          "Event multiplicity (complete events);"
                          "Multiplicity;Counts",
                          36, -0.5, 35.5);
  TH2F *h2_R_vs_L[18] = {nullptr};
  for (Int_t s = 1; s <= 16; s++) {
    h2_R_vs_L[s] = new TH2F(
        Form("h2_R_vs_L_s%d_Normed", s),
        Form(";Strip %d L #DeltaE [a.u.];Strip %d R #DeltaE [a.u.]", s, s), 200,
        strip_e_min, strip_e_max, 200, strip_e_min, strip_e_max);
  }
  // Optional 1D histograms based on hardware flags.
  TH1F *h1_cathode = nullptr;
  TH1F *h1_strip17 = nullptr;
  if (Constants::cfg.HAS_CATHODE)
    h1_cathode = new TH1F("h1_cathode_Normed", ";Cathode #DeltaE [a.u.];Counts",
                          400, 0.0, 1.0);
  if (Constants::cfg.HAS_STRIP17)
    h1_strip17 = new TH1F("h1_strip17_Normed", ";Strip17 #DeltaE [a.u.];Counts",
                          400, strip_e_min, strip_e_max);
  TH2F *h2_strip0_vs_grid = nullptr;
  TH1F *h1_strip0 = nullptr;
  TH1F *h1_grid = nullptr;
  if (Constants::cfg.HAS_STRIP0 && Constants::cfg.HAS_GRID) {
    h2_strip0_vs_grid = new TH2F("h2_strip0_vs_grid_Normed",
                                 ";Grid #DeltaE [a.u.];Strip0 #DeltaE [a.u.]",
                                 200, 0.0, 1.0, 200, strip_e_min, strip_e_max);
  } else if (Constants::cfg.HAS_STRIP0) {
    h1_strip0 = new TH1F("h1_strip0_Normed", ";Strip0 #DeltaE [a.u.];Counts",
                         400, strip_e_min, strip_e_max);
  } else if (Constants::cfg.HAS_GRID) {
    h1_grid = new TH1F("h1_grid_Normed", ";Grid #DeltaE [a.u.];Counts", 400,
                       0.0, 1.0);
  }

  Long64_t n_entries = input_tree->GetEntries();
  std::cout << "[" << file_label << "] building normed summary over "
            << n_entries << " events..." << std::endl;
  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    ev.Decode();

    for (Int_t s = 0; s < 18; s++)
      h_music->Fill(Double_t(s), ev.total[s]);

    for (Int_t s = 1; s <= 16; s++)
      h2_R_vs_L[s]->Fill(ev.right[s], ev.left[s]);

    if (h1_cathode && ev.cathode > 0.0)
      h1_cathode->Fill(ev.cathode);

    if (h1_strip17)
      h1_strip17->Fill(ev.total[17]);

    if (h2_strip0_vs_grid)
      h2_strip0_vs_grid->Fill(ev.grid, ev.total[0]);

    if (h1_strip0)
      h1_strip0->Fill(ev.total[0]);

    if (h1_grid)
      h1_grid->Fill(ev.grid);

    Int_t mult = 0;
    for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
      mult += ev.hits_adc[k];
    h_mult->Fill(Double_t(mult));
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString subdir = "events_summary_normed/" + file_label;

    // hMUSIC
    TCanvas *c_music = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music, c_music);
    h_music->GetYaxis()->SetTitleOffset(1.4);
    c_music->SetLeftMargin(0.18);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music, "music_strip_energies_normed", subdir,
                                PlotSaveOptions::kLINEAR);
    input_file->cd();
    c_music->Write(h_music->GetName(), TObject::kOverwrite);
    delete c_music;

    // hMult
    TCanvas *c_mult = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h_mult, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_mult, "multiplicity_normed", subdir,
                                PlotSaveOptions::kLOG);
    delete c_mult;

    // R vs L for each split strip
    for (Int_t s = 1; s <= 16; s++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_R_vs_L[s], c);
      h2_R_vs_L[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, Form("R_vs_L_s%d_normed", s), subdir,
                                  PlotSaveOptions::kLINEAR);
      input_file->cd();
      c->Write(h2_R_vs_L[s]->GetName(), TObject::kOverwrite);
      delete c;
    }

    // Optional cathode histogram
    if (h1_cathode) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_cathode, kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, "cathode_normed", subdir,
                                  PlotSaveOptions::kLOG);
      input_file->cd();
      h1_cathode->Write(h1_cathode->GetName(), TObject::kOverwrite);
      delete c;
    }

    // Optional strip17 histogram
    if (h1_strip17) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_strip17, kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, "strip17_normed", subdir,
                                  PlotSaveOptions::kLOG);
      input_file->cd();
      h1_strip17->Write(h1_strip17->GetName(), TObject::kOverwrite);
      delete c;
    }

    // Strip0 vs grid (2D) or strip0/grid (1D)
    if (h2_strip0_vs_grid) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_strip0_vs_grid, c);
      h2_strip0_vs_grid->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, "strip0_vs_grid_normed", subdir,
                                  PlotSaveOptions::kLINEAR);
      input_file->cd();
      c->Write(h2_strip0_vs_grid->GetName(), TObject::kOverwrite);
      delete c;
    }
    if (h1_strip0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_strip0, kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, "strip0_normed", subdir,
                                  PlotSaveOptions::kLOG);
      input_file->cd();
      h1_strip0->Write(h1_strip0->GetName(), TObject::kOverwrite);
      delete c;
    }
    if (h1_grid) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_grid, kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, "grid_normed", subdir,
                                  PlotSaveOptions::kLOG);
      input_file->cd();
      h1_grid->Write(h1_grid->GetName(), TObject::kOverwrite);
      delete c;
    }
  }

  for (Int_t s = 1; s <= 16; s++)
    delete h2_R_vs_L[s];
  delete h_music;
  delete h_mult;
  delete h1_cathode;
  delete h1_strip17;
  delete h2_strip0_vs_grid;
  delete h1_strip0;
  delete h1_grid;

  input_file->Close();
  delete input_file;
}

TGraph *TraceCreator::BuildTraceFromTotals(const Double_t *total) {
  Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  Int_t n_pts = s_hi - s_lo + 1;
  TGraph *g = new TGraph(n_pts);
  for (Int_t k = 0; k < n_pts; k++)
    g->SetPoint(k, s_lo + k, total[s_lo + k]);
  return g;
}
