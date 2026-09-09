#ifndef TAG_EFFICIENCY_HPP
#define TAG_EFFICIENCY_HPP
/**
 * @file TagEfficiency.hpp
 * @brief The tag-efficiency correction, and its contract with the cross
 * section.
 *
 * Whatever method produces it, the efficiency side hands the cross section one
 * record per channel and reaction strip: the count it applies to, the fraction
 * of true reactions at that strip that end up in that count, and the fraction
 * counted one strip late instead.
 *
 * The cross section then applies one formula,
 *
 *     N_true(r) = (n_counted(r) - migrate(r-1) * N_true(r-1)) / eff(r)
 *
 * and never needs to know how the numbers were made. Change the method, keep
 * the record, and the cross section keeps working.
 */
#include "RegionCuts.hpp"
#include <Rtypes.h>
#include <TCutG.h>
#include <TString.h>
#include <map>
#include <vector>

class TFile;
class TH2F;
class TTree;
struct CrossSectionChannel;

/**
 * @brief One channel and strip's efficiency, as the cross section consumes it.
 *
 * The whole contract between the two sides. See the file overview for the
 * formula these feed.
 */
struct TagEfficiencyRecord {
  Int_t reac = -1; ///< Reaction strip index; `-1` for an unset record.
  /// The numerator this efficiency belongs to: the events counted at this strip
  /// by the same selection the efficiency was measured for. Pairing an
  /// efficiency with a differently-selected count is the error this field
  /// exists to prevent.
  Double_t n_counted = 0.0;
  /// Fraction of true reactions at this strip that end up in #n_counted.
  Double_t eff = 0.0;
  /// Uncertainty on #eff, or `0` when the method provides none.
  Double_t eff_err = 0.0;
  /// Fraction of true reactions at this strip counted at strip `reac + 1`
  /// instead, which the cross section unfolds.
  Double_t migrate = 0.0;
  /// Of #eff, the part the tag alone passes; the remainder is the region cut.
  Double_t tag_eff = 0.0;
};

/// @brief Persistence for the efficiency records.
namespace TagEfficiencyStore {

/// @brief Path to the store.
/// @return `root_files/tag_efficiency.root` under the dataset's results.
TString Path();
/**
 * @brief Replace one channel's records, keeping every other channel's.
 * @param channel Channel whose records to replace.
 * @param records The new records, one per reaction strip.
 * @param method  Label identifying how they were produced, stamped into the
 *                store so a later reader knows what it is holding.
 */
void Write(const TString &channel,
           const std::vector<TagEfficiencyRecord> &records,
           const TString &method);
/**
 * @brief Read one record.
 * @param channel Channel to look up.
 * @param reac    Reaction strip index.
 * @param[out] record Set on success; untouched otherwise.
 * @return `kFALSE` when the store holds no such record.
 */
Bool_t Load(const TString &channel, Int_t reac, TagEfficiencyRecord &record);
/// @brief The method label stamped into the store.
/// @return The label, or an empty string when no store exists.
TString Method();
} // namespace TagEfficiencyStore

/**
 * @brief The bootstrap method for measuring the tag efficiency.
 *
 * Takes the strip's mean trace of the channel's events **from the data**,
 * resamples it with the measured per-strip beam widths, and pushes the result
 * through the real tag and the channel's real hand-drawn region cut. Running
 * the actual selection on realistic pseudo-events is what makes the resulting
 * efficiency apply to the count it is paired with.
 *
 * Runs once per configured channel.
 */
class TagEfficiency {
public:
  /**
   * @brief Measure every configured channel and write the store and figures.
   * @return `kFALSE` when a prerequisite is missing — the scatter cache, the
   *         noise sigmas, or the hand-drawn cuts.
   */
  Bool_t Run();

  /// @brief The label stamped into the store for this method.
  static const char *MethodLabel() {
    return "bootstrap: mean trace inside the hand-drawn cut, resampled with "
           "the measured per-strip beam widths";
  }

private:
  struct Strip {
    TCutG *cut = nullptr;       // the channel's hand-drawn cut at this strip
    std::vector<Double_t> mean; // mean trace of the events inside it
    Double_t count = 0.0;       // how many
  };
  struct Outcome {
    Long64_t n = 0, tagged = 0, inside = 0, next_inside = 0;
  };

  Bool_t LoadSigmas(TFile &cache);
  // One channel end to end; kTRUE when it produced at least one record.
  Bool_t RunChannel(const CrossSectionChannel &ch, TTree *reservoir);
  void LoadDrawnCuts(const TString &region);
  void MeanTraces(TTree *reservoir);
  Outcome Bootstrap(Int_t reac, TH2F *plane);
  void DrawMeanTrace(const TString &subdir, Int_t reac, const Strip &st);
  void DrawPlane(const TString &subdir, Int_t reac, TH2F *plane, TCutG *cut);

  std::map<Int_t, Strip> strips_; // the channel being run
};

#endif
