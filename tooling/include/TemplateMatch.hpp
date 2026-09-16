#ifndef TEMPLATE_MATCH_HPP
#define TEMPLATE_MATCH_HPP
/**
 * @file TemplateMatch.hpp
 * @brief The template-match second pass: an (a,n) count per strip that does
 * not go through the tag, to set against the tagged count.
 *
 * The tag is a chain of per-strip conditions, each with an efficiency that is
 * only estimated. This pass counts the same reaction a different way. From
 * the tagged events inside each strip's (a,n) region it takes the mean trace,
 * the **template** for a reaction at that strip; then it goes back over every
 * event in the data and asks which template, if any, fits the trace better
 * than the flat beam does, by chi-square over the strips with the measured
 * per-strip beam widths. An event is counted at the strip whose template fits
 * best when that fit beats the beam by more than `TEMPLATE_DELTA_CHI2`.
 *
 * Both methods then get the same treatment: a bootstrap that resamples the
 * template with the beam widths and pushes it through the classifier gives
 * this method's efficiency and its migration to the next strip, and the same
 * resampling of a flat trace gives its false-positive rate on beam. The
 * report sets the two corrected counts side by side, per strip.
 *
 * Only a comparison is written; the cross section keeps using the tag.
 */
#include <Rtypes.h>
#include <vector>

/// @brief One strip's outcome, both methods.
struct TemplateStripResult {
  Int_t reac = -1;
  Double_t n_template_events = 0.0; ///< Events the template was averaged from.
  Long64_t n_tagged = 0;            ///< Tagged at this strip (cache).
  Double_t eff_tag = 0.0;       ///< Bootstrap tag efficiency, 0 when no store.
  Double_t corrected_tag = 0.0; ///< Region count over eff_tag.
  Long64_t n_matched = 0;       ///< Classified at this strip by the templates.
  Long64_t n_both = 0;          ///< Of those, also tagged at this strip.
  Double_t eff_template = 0.0;  ///< Bootstrap: template resamples classed here.
  Double_t migrate_template = 0.0;   ///< Classed one strip late instead.
  Double_t false_positive = 0.0;     ///< Flat-beam resamples classed here.
  Double_t background = 0.0;         ///< false_positive x beam count at strip.
  Double_t corrected_template = 0.0; ///< (n_matched - background) / eff.
};

/**
 * @brief Build the templates, bootstrap them, run the data pass, and report.
 */
class TemplateMatch {
public:
  /// @return `kFALSE` when the scatter cache, its noise sigmas or the
  ///         reservoir are missing, or no strip yields a template.
  Bool_t Run();

private:
  std::vector<std::vector<Double_t>>
      templates_; ///< By ReacIndex; empty = none.
  std::vector<TemplateStripResult> results_;
  void DrawTemplate(Int_t reac, const std::vector<Double_t> &tpl,
                    Double_t n_events);
  void WriteReport(Long64_t seen, Long64_t considered, const char *tag_part);
};

#endif
