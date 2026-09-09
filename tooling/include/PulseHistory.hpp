#ifndef PULSE_HISTORY_HPP
#define PULSE_HISTORY_HPP

// Pole-zero pulse-history correction on the raw hit stream
#include "BinaryUtils.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <vector>

class TH1D;
class TH2D;
class TProfile;

namespace PulseHistory {

const UInt_t kFlagClamped = 0x2000;
const Int_t kNBins = 12;
const Double_t kLogLo = -6.0; // log10(dt = 1us in seconds)
const Double_t kLogHi = -4.0; // log10(dt = 100us in seconds)
const Int_t kMaxAmpBins = 6;

// Channel groups, one kernel each: the long and short end of the split strips
// on either chain, and the two single-pad guards. A short end fires in only a
// fraction of beam events and its height follows the track position, so its
// kernel is fitted on the fired hits alone and describes a smaller share of
// its variance.
enum Group {
  kNone = 0,
  kLongLeft = 1,
  kLongRight = 2,
  kShortLeft = 3,
  kShortRight = 4,
  kGuard0 = 5,
  kGuard17 = 6,
  kNGroups = 7
};
const char *GroupName(Int_t g);
// Short tag for file names: L, R, Ls, Rs, S0, S17.
const char *GroupTag(Int_t g);
// 0 for the L chain, 1 for the R chain, -1 for the guards.
Int_t ChainOf(Int_t g);
Bool_t IsLongGroup(Int_t g);

struct Kernel {
  Bool_t ok = kFALSE;
  Int_t n_amp = 1;
  Double_t k[kMaxAmpBins][kNBins]; // [amplitude band][dt bin]
  Double_t intercept = 0.0;
  Double_t r2 = 0.0;
  Double_t rms_before = 0.0, rms_after = 0.0;
  Long64_t n = 0;      // (event, channel) pairs in the fit
  Long64_t n_beam = 0; // beam-like events for this group's selection
  Kernel();
};

struct Result {
  Kernel kernel[kNGroups];
  Long64_t n_hits = 0, n_seeds = 0, n_beam_events = 0;
  Long64_t n_corrected = 0, n_clamped = 0;
  Long64_t n_clamped_group[kNGroups]; // clamps per group
  Double_t mean_shift[kNGroups];      // mean of the subtracted term, ADC
  Bool_t sorted_input = kTRUE;
  // Beam peak per (board, channel) index, ADC. Zero where not a long end.
  std::vector<Double_t> mode;
  // Diagnostics, filled when SAVE_PLOTS is on and drawn by SavePlots, which
  // also frees them. Per group: the channel deviation against the predicted
  // shift, its distribution before and after, the deviation against the time
  // to the immediately previous pulse before and after, and the shift
  // applied per hit.
  TH2D *dev_vs_pred[kNGroups];
  TH1D *dev_before[kNGroups], *dev_after[kNGroups], *shift[kNGroups];
  TProfile *dtprev_before[kNGroups], *dtprev_after[kNGroups];
  Result();
  ~Result();
};

// Bin index of a dt in seconds, or -1 outside the kernel's range.
Int_t BinOf(Double_t dt_s);
// Centre of a bin in microseconds, for reports.
Double_t BinCentreUs(Int_t b);
// Amplitude band of a previous pulse, given the channel's beam peak.
Int_t AmpBinOf(Double_t e_prev, Double_t mode, Int_t n_amp);

// Group of a (board, channel) under the active channel map.
std::vector<Int_t> BuildGroupMap();

// Measure the kernels on this subfile's beam-like events. Hits must be time
// ordered (FOR NOW they are sorted in place if not, which the event builder
// needs anyway). Returns kFALSE when no group could be fitted. TO DO: always
// pass through pulse history, sort here prior to event building by default so
// that the CUDA sort can be used
Bool_t Measure(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
               Result &res, const TString &file_label);

// Apply the measured kernels in place. Only bins up to
// Constants::cfg.PULSE_HISTORY_APPLY_MAX_US are used; beyond that the fitted
// coefficients are degenerate with the intercept and physically zero.
void Apply(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
           Result &res);

TString Report(const Result &res, const TString &file_label);
void SavePlots(Result &res, const TString &file_label);
void WriteToEventsFile(const TString &events_subpath, const Result &res);

} // namespace PulseHistory

#endif
