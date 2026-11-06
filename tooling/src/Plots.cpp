#include "Plots.hpp"

#include <TSystem.h>
#include <iostream>

void Plots::SaveFigure(TCanvas *canvas, TString output_name,
                       TString output_subdirectory,
                       PlotSaveOptions save_options) {
  PlottingUtils::SaveFigure(canvas, output_name, output_subdirectory,
                            save_options);

  // The path PlottingUtils wrote under: plots base, then the optional
  // subdirectory, then the file.
  const TString base = PlottingUtils::GetPlotsBaseDir();
  const TString full_dir =
      output_subdirectory == "" ? base : base + "/" + output_subdirectory;

  FilePdf(full_dir, output_name);
  if (save_options != PlotSaveOptions::kLINEAR)
    FilePdf(full_dir, "log_" + output_name);
}

void Plots::FilePdf(const TString &full_dir, const TString &output_name) {
  const TString from = full_dir + "/" + output_name + ".pdf";
  // AccessPathName answers 0 when the path *is* reachable, so a false return
  // means the file is there; a kPNG-only format writes no PDF and falls out.
  if (gSystem->AccessPathName(from))
    return;
  const TString pdf_dir = full_dir + "/pdf";
  if (gSystem->AccessPathName(pdf_dir))
    gSystem->mkdir(pdf_dir, kTRUE);
  if (gSystem->Rename(from, pdf_dir + "/" + output_name + ".pdf"))
    std::cerr << "Plots: could not move " << from << " into " << pdf_dir
              << std::endl;
}
