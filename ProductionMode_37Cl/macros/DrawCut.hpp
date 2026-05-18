#ifndef DRAW_CUT_HPP
#define DRAW_CUT_HPP

#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include <TCanvas.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <vector>

namespace DrawCutNS {

inline void DrawCutOneFile(const TString &hist_name, const TString &cut_name,
                           const FileSpec &spec) {
  TString events_path = EventsName(spec) + ".root";
  TString file_label = FileLabel(spec);

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *file = IO::OpenForWriting(events_path, "UPDATE");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << events_path << std::endl;
    return;
  }

  TH2 *hist = static_cast<TH2 *>(file->Get(hist_name));
  if (!hist) {
    std::cerr << "No histogram '" << hist_name << "' in " << events_path
              << std::endl;
    file->Close();
    return;
  }

  gROOT->SetBatch(kFALSE);

  TCanvas *c = new TCanvas("c_drawcut", hist->GetTitle(), 900, 700);
  hist->Draw("COLZ");
  c->SetLogz();
  c->Update();

  std::cout << "[" << file_label << "] Drawing region '" << cut_name << "' on '"
            << hist_name << "'" << std::endl;
  std::cout << "Click points around the region; double-click to close."
            << std::endl;

  TCutG *cut = static_cast<TCutG *>(gPad->WaitPrimitive("CUTG", "CutG"));
  if (!cut) {
    std::cerr << "No cut drawn" << std::endl;
    file->Close();
    return;
  }

  cut->SetName(cut_name);
  cut->SetTitle(
      Form("region %s drawn on %s", cut_name.Data(), hist_name.Data()));
  // VarX is reused as the source-hist tag; ComputeCutXY reads it back.
  cut->SetVarX(hist_name);
  cut->SetVarY("");
  cut->SetLineColor(kRed + 1);
  cut->SetLineWidth(3);

  // Reapply the histogram + cut on top so the canvas leaves CutG editor
  // mode and the data is visible behind the polygon. Drag points to adjust.
  hist->Draw("COLZ");
  cut->Draw("L SAME");
  gPad->Modified();
  gPad->Update();

  std::cout << "Adjust polygon points if needed, then press ENTER to save..."
            << std::flush;
  std::cin.get();

  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir)
    cuts_dir = file->mkdir("cuts");
  cuts_dir->cd();
  cut->Write(cut_name, TObject::kOverwrite);

  std::cout << "Saved cut '" << cut_name << "' to " << events_path << ":/cuts/"
            << std::endl;

  delete c;
  file->Close();
}

} // namespace DrawCutNS

inline void DrawCut(TString hist_name, TString cut_name,
                    TString file_label = "") {
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
    DrawCutNS::DrawCutOneFile(hist_name, cut_name, specs[k]);
}

#endif
