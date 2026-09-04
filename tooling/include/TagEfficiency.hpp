#ifndef TAG_EFFICIENCY_HPP
#define TAG_EFFICIENCY_HPP
// The tag-efficiency correction and the contract between it and the cross
// section.
//
// Whatever method produces it, the efficiency side hands the cross section one
// record per channel and reaction strip: the count it applies to, the fraction
// of true reactions at that strip which end up in that count, and the fraction
// counted one strip late instead. The cross section then applies one formula,
//
//   N_true(r) = (n_counted(r) - migrate(r-1) * N_true(r-1)) / eff(r)
//
// and never needs to know how the numbers were made. Change the method, keep
// the record, and the cross section keeps working.
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

struct TagEfficiencyRecord {
  Int_t reac = -1;
  // The numerator this efficiency belongs to: the events counted at this strip
  // by the same selection the efficiency was measured for.
  Double_t n_counted = 0.0;
  // Fraction of true reactions at this strip that end up in n_counted, and its
  // uncertainty (0 when the method gives none).
  Double_t eff = 0.0;
  Double_t eff_err = 0.0;
  // Fraction of true reactions at this strip counted at strip reac+1 instead.
  Double_t migrate = 0.0;
  // Of eff, the part the tag alone passes (the rest is the region cut).
  Double_t tag_eff = 0.0;
};

namespace TagEfficiencyStore {
// root_files/tag_efficiency.root in the dataset's results.
TString Path();
// Replace one channel's records, stamped with the method's label. Other
// channels' records in the store are kept.
void Write(const TString &channel,
           const std::vector<TagEfficiencyRecord> &records,
           const TString &method);
// The record for a channel and strip; kFALSE (record untouched) when the store
// has none.
Bool_t Load(const TString &channel, Int_t reac, TagEfficiencyRecord &record);
// The method label stamped in the store, or "" when there is no store.
TString Method();
} // namespace TagEfficiencyStore

// The bootstrap method: the strip's mean trace of the channel's events from
// the data, resampled with the measured per-strip beam widths, through the
// real tag and the channel's real region cut. Runs once per configured
// channel.
class TagEfficiency {
public:
  // Runs every channel, writes the store and the figures. kFALSE when a
  // prerequisite is missing (cache, noise sigmas, drawn cuts).
  Bool_t Run();

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
