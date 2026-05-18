#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <vector>

namespace CathodeAnodeRatioNS {

const Int_t kAnodeBins = 400;
const Int_t kRatioBins = 400;
const Double_t kRatioMin = 0.0;
const Double_t kRatioMax = 2.0;

void CathodeAnodeRatioOneFile(const FileSpec &spec) {
  TString events_path = EventsName(spec) + ".root";
  TString file_label = FileLabel(spec);

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *file = IO::OpenForReading(events_path);
  if (!file || file->IsZombie()) {
    std::cerr << "[" << file_label << "] Cannot open " << events_path
              << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("events"));
  if (!tree) {
    std::cerr << "[" << file_label << "] No events tree in " << events_path
              << std::endl;
    file->Close();
    return;
  }

  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_mev)
    std::cerr << "[" << file_label
              << "] WARNING: no MeV branches; using uncalibrated ADC values"
              << std::endl;
  const char *unit = ev.Unit();
  const Double_t anode_min =
      ev.is_mev ? Constants::STRIP_E_MIN_MEV : Constants::STRIP_E_MIN_ADC;
  const Double_t anode_max =
      ev.is_mev ? Constants::STRIP_E_MAX_MEV : Constants::STRIP_E_MAX_ADC;

  std::vector<TH2F *> h_per_strip(18, static_cast<TH2F *>(nullptr));
  for (Int_t s = 0; s < 18; s++) {
    TString name = Form("h2_cathRatio_vs_strip%d_%s", s, file_label.Data());
    TString title =
        Form("Strip %d: Cathode / #Sigma_{strips} vs Strip %d signal (%s);"
             "Strip %d #DeltaE [%s];Cathode / #Sigma anode",
             s, s, file_label.Data(), s, unit);
    h_per_strip[s] = new TH2F(name, title, kAnodeBins, anode_min, anode_max,
                              kRatioBins, kRatioMin, kRatioMax);
  }

  TH2F *h_cath_over_strip_vs_strip =
      new TH2F(Form("h2_cathOverStrip_vs_stripNum_%s", file_label.Data()),
               Form("Cathode / single-strip anode vs Strip number (%s);"
                    "Strip number;Cathode / Strip #DeltaE",
                    file_label.Data()),
               18, -0.5, 17.5, kRatioBins, kRatioMin, kRatioMax);

  Long64_t n_entries = tree->GetEntries();
  Long64_t n_no_cathode = 0;
  Long64_t n_zero_anode = 0;
  Long64_t n_used = 0;

  std::cout << "[" << file_label << "] Cathode/anode ratio over " << n_entries
            << " events..." << std::endl;

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    ev.Decode();

    if (!ev.is_mev && ev.cathode_adc == -1) {
      n_no_cathode++;
      continue;
    }

    Double_t anode_sum = 0;
    for (Int_t s = 0; s < 18; s++)
      anode_sum += ev.total[s];

    if (anode_sum <= 0) {
      n_zero_anode++;
      continue;
    }

    Double_t ratio = ev.cathode / anode_sum;
    for (Int_t s = 0; s < 18; s++)
      h_per_strip[s]->Fill(ev.total[s], ratio);

    for (Int_t s = 0; s < 18; s++) {
      if (ev.total[s] <= 0)
        continue;
      h_cath_over_strip_vs_strip->Fill(Double_t(s), ev.cathode / ev.total[s]);
    }

    n_used++;

    if ((j + 1) % 100000 == 0)
      std::cout << "  " << (j + 1) << "/" << n_entries << std::endl;
  }

  Double_t inv_total = n_entries > 0 ? 1.0 / Double_t(n_entries) : 0.0;
  std::cout << "[" << file_label << "] entries=" << n_entries
            << "  no-cathode=" << n_no_cathode << " ("
            << TString::Format("%.4f", n_no_cathode * inv_total) << ")"
            << "  zero-anode-skipped=" << n_zero_anode << "  used=" << n_used
            << std::endl;

  TString out_subdir = "cathode_ratio/" + file_label;

  for (Int_t s = 0; s < 18; s++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_per_strip[s], c);
    PlottingUtils::SaveFigure(c, Form("cathRatio_vs_strip%02d", s), out_subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
    delete h_per_strip[s];
  }

  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_cath_over_strip_vs_strip, c);
    PlottingUtils::SaveFigure(c, "cathOverStrip_vs_stripNum", out_subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  delete h_cath_over_strip_vs_strip;

  file->Close();
}

} // namespace CathodeAnodeRatioNS

void CathodeAnodeRatio(TString file_label = "") {
  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = BuildFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from BuildFileSpecs()" << std::endl;
      return;
    }
  } else {
    FileSpec s = ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  for (Int_t k = 0; k < Int_t(specs.size()); k++)
    CathodeAnodeRatioNS::CathodeAnodeRatioOneFile(specs[k]);
}
