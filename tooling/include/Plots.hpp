#ifndef PLOTS_HPP
#define PLOTS_HPP

#include "PlottingUtils.hpp"
#include <TString.h>

/**
 * @file Plots.hpp
 * @brief The one figure-save call, with the PDFs filed under the plot
 *        directory's `pdf/` subfolder.
 *
 * Every binary sets the style once (InitUtils::SetROOTPreferences, save
 * format kPNGAndPDF), so PlottingUtils::SaveFigure() writes each figure two
 * ways, as `.png` and as `.pdf`, into the same directory. SaveFigure() here
 * is that same call, then moves the PDFs it wrote into `<plot dir>/pdf/` —
 * the `log_` variants included — so a plot directory keeps only PNGs and
 * every vector copy sits in one predictable place. Under a single-format
 * style the move either files the one PDF or, for kPNG, does nothing.
 */
class Plots {
public:
  /**
   * @brief PlottingUtils::SaveFigure(), then file the PDFs it wrote.
   *
   * The figure lands at `<plots base>/<subdirectory>/<name>.{png,pdf}`, and
   * the log-y variant at `<plots base>/<subdirectory>/log_<name>.{png,pdf}`
   * when @p save_options asks for one; afterwards each written PDF is moved
   * to the same name under `<plots base>/<subdirectory>/pdf/`. Missing
   * directories are created. Everything else — style, batch graphics, the
   * Portrait PDF pages — is PlottingUtils's, unchanged.
   *
   * @param canvas              Canvas to save. Must not be null.
   * @param output_name         Base filename, without extension.
   * @param output_subdirectory Subdirectory under the plots base; empty means
   *                            the base itself.
   * @param save_options        Which axis scalings to write, kBOTH by
   *                            default.
   */
  static void SaveFigure(TCanvas *canvas, TString output_name,
                         TString output_subdirectory = "",
                         PlotSaveOptions save_options = PlotSaveOptions::kBOTH);

private:
  /// @brief Move `<full_dir>/<output_name>.pdf` into `<full_dir>/pdf/`,
  ///        creating that subfolder, when the file exists.
  static void FilePdf(const TString &full_dir, const TString &output_name);
};

#endif
