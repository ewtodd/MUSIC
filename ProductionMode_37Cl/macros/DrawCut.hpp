#ifndef DRAW_CUT_HPP
#define DRAW_CUT_HPP

#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <sys/select.h>
#include <unistd.h>
#include <vector>

namespace DrawCutNS {

// Block until ENTER is pressed on stdin, while still pumping ROOT GUI
// events every 50ms so the canvas remains interactive (right-click menus,
// vertex drags, etc. all need ProcessEvents to fire continuously).
inline void WaitForEnterPumpingEvents() {
  fd_set rfds;
  struct timeval tv;
  while (true) {
    FD_ZERO(&rfds);
    FD_SET(0, &rfds); // stdin
    tv.tv_sec = 0;
    tv.tv_usec = 50000; // 50 ms
    int rv = select(1, &rfds, nullptr, nullptr, &tv);
    if (rv > 0 && FD_ISSET(0, &rfds)) {
      std::cin.get();
      return;
    }
    gSystem->ProcessEvents();
  }
}

} // namespace DrawCutNS

namespace DrawCutNS {

// Draw an interactive TCutG on a histogram already present in `file`, save
// the cut to `file:/cuts/`. The file must be open for writing (UPDATE).
// `label` is used only for the console banner.
inline Bool_t DrawCutOnFile(TFile *file, const TString &hist_name,
                            const TString &cut_name, const TString &label) {
  TH2 *hist = static_cast<TH2 *>(file->Get(hist_name));
  if (!hist) {
    std::cerr << "No histogram '" << hist_name << "' in " << file->GetName()
              << std::endl;
    return kFALSE;
  }

  // Force-disable batch and ensure a TApplication exists so the GUI event
  // loop actually runs WaitPrimitive's clicks reach the canvas.
  gROOT->SetBatch(kFALSE);
  if (!gApplication) {
    static Int_t app_argc = 1;
    static char app_arg0[] = "drawcut";
    static char *app_argv[] = {app_arg0};
    new TApplication("drawcut", &app_argc, app_argv);
  }

  // Unique canvas name per call so a stale canvas from a previous iteration
  // can't shadow the new one in gROOT's list.
  static Int_t s_canvas_counter = 0;
  TString canvas_name = Form("c_drawcut_%d", ++s_canvas_counter);
  TCanvas *c = new TCanvas(canvas_name, hist->GetTitle(), 900, 700);
  c->cd();
  hist->Draw("COLZ");
  c->SetLogz();
  c->Update();
  gSystem->ProcessEvents();

  std::cout << "[" << label << "] Drawing region '" << cut_name << "' on '"
            << hist_name << "'" << std::endl;
  std::cout << "Click points around the region; double-click to close."
            << std::endl;

  // WaitPrimitive on the explicit canvas rather than gPad — gPad can be
  // stale across iterations.
  TCutG *cut = static_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
  if (!cut) {
    std::cerr << "No cut drawn" << std::endl;
    delete c;
    return kFALSE;
  }

  cut->SetName(cut_name);
  cut->SetTitle(
      Form("region %s drawn on %s", cut_name.Data(), hist_name.Data()));
  // VarX is reused as the source-hist tag; ComputeCutXY reads it back.
  cut->SetVarX(hist_name);
  cut->SetVarY("");
  cut->SetLineColor(kRed + 1);
  cut->SetLineWidth(3);
  cut->SetMarkerColor(kRed + 1);
  cut->SetMarkerStyle(20);
  cut->SetMarkerSize(1.6);

  // Reapply the histogram + cut on top so the canvas leaves CutG editor
  // mode and the data is visible behind the polygon.
  c->cd();
  hist->Draw("COLZ");
  cut->Draw("LP SAME");
  c->SetEditable(kTRUE);
  c->Modified();
  c->Update();
  gSystem->ProcessEvents();

  std::cout << "Right-click on the polygon for the editor menu (InsertPoint,"
               " RemovePoint, SetPoint via x/y dialog)."
            << std::endl;
  std::cout << "Press ENTER when ready to save..." << std::flush;
  WaitForEnterPumpingEvents();

  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir)
    cuts_dir = file->mkdir("cuts");
  cuts_dir->cd();
  cut->Write(cut_name, TObject::kOverwrite);
  // Flush the directory key index so the cut survives an interrupted run.
  cuts_dir->SaveSelf(kTRUE);

  std::cout << "Saved cut '" << cut_name << "' to " << file->GetName()
            << ":/cuts/" << std::endl;

  delete c;
  return kTRUE;
}

// Like DrawCutOnFile, but pre-seeds the canvas with a clone of `template_cut`
// (e.g. from a related channel) for the user to drag-adjust. Press ENTER
// to accept the (possibly edited) polygon and save under `cut_name`.
inline Bool_t DrawCutOnFileFromTemplate(TFile *file, const TString &hist_name,
                                        const TString &cut_name,
                                        TCutG *template_cut,
                                        const TString &label) {
  TH2 *hist = static_cast<TH2 *>(file->Get(hist_name));
  if (!hist) {
    std::cerr << "No histogram '" << hist_name << "' in " << file->GetName()
              << std::endl;
    return kFALSE;
  }
  if (!template_cut) {
    std::cerr << "No template cut for '" << cut_name << "'" << std::endl;
    return kFALSE;
  }

  gROOT->SetBatch(kFALSE);
  if (!gApplication) {
    static Int_t app_argc = 1;
    static char app_arg0[] = "drawcut";
    static char *app_argv[] = {app_arg0};
    new TApplication("drawcut", &app_argc, app_argv);
  }

  static Int_t s_tmpl_counter = 0;
  TString canvas_name = Form("c_drawcut_tmpl_%d", ++s_tmpl_counter);
  TCanvas *c = new TCanvas(canvas_name, hist->GetTitle(), 900, 700);
  c->cd();
  hist->Draw("COLZ");
  c->SetLogz();

  TCutG *cut = static_cast<TCutG *>(template_cut->Clone(cut_name));
  cut->SetName(cut_name);
  cut->SetTitle(Form("region %s (from %s) on %s", cut_name.Data(),
                     template_cut->GetName(), hist_name.Data()));
  cut->SetVarX(hist_name);
  cut->SetVarY("");
  cut->SetLineColor(kRed + 1);
  cut->SetLineWidth(3);
  cut->SetMarkerColor(kRed + 1);
  cut->SetMarkerStyle(20);
  cut->SetMarkerSize(1.6);
  cut->Draw("LP SAME");
  c->SetEditable(kTRUE);
  c->Modified();
  c->Update();
  gSystem->ProcessEvents();

  std::cout << "[" << label << "] '" << cut_name << "' on '" << hist_name
            << "' (cloned from " << template_cut->GetName() << ")" << std::endl;
  std::cout << "Right-click on the polygon for the editor menu (InsertPoint,"
               " RemovePoint, SetPoint via x/y dialog)."
            << std::endl;
  std::cout << "Press ENTER when ready to save..." << std::flush;
  WaitForEnterPumpingEvents();

  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir)
    cuts_dir = file->mkdir("cuts");
  cuts_dir->cd();
  cut->Write(cut_name, TObject::kOverwrite);
  cuts_dir->SaveSelf(kTRUE);
  std::cout << "Saved cut '" << cut_name << "' to " << file->GetName()
            << ":/cuts/" << std::endl;
  delete c;
  return kTRUE;
}

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

  DrawCutOnFile(file, hist_name, cut_name, file_label);
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
