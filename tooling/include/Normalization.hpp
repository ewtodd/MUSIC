#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TTree.h>

/**
 * @brief A view over one event's energies, decoding raw ADC into calibrated
 * units.
 *
 * Binds to an events tree or chain, then decodes each entry: per-channel gains
 * are applied, then per-strip alignment factors. Whether the result is in
 * arbitrary calibrated units or still raw ADC depends on whether a calibration
 * tree was found — see #is_normed.
 *
 * The layout follows the detector: strips 1-16 are segmented and read at a
 * left and a right end, held in arrays of 16 indexed by `strip - 1`; strips 0
 * and 17 are unsegmented and each a single value. A strip's total deposit is
 * never stored; Total() adds the two ends when asked.
 *
 * @note Holds a non-owning pointer to the tree it is attached to, which must
 *       outlive it. Gains are reloaded automatically when a `TChain` rolls over
 *       to a new file, tracked by #loaded_tree_.
 */
struct EnergyView {
  /// @name Raw ADC, straight from the events tree
  /// @{
  UShort_t leftdE_adc[16];  ///< Left end of strips 1-16, index `strip - 1`.
  UShort_t rightdE_adc[16]; ///< Right end of strips 1-16, index `strip - 1`.
  UShort_t strip0_adc;      ///< Strip 0, unsegmented.
  UShort_t strip17_adc;     ///< Strip 17, unsegmented.
  UShort_t hits_adc[36];    ///< Raw per-slot hit ADC values.
  Short_t cathode_adc;      ///< Raw cathode ADC value.
  Short_t grid_adc;         ///< Raw Frisch grid ADC value.
  /// @}

  /// @name Calibration, from the file's calibration tree
  /// @{
  Float_t gain_left[16];  ///< Left-end gain of strips 1-16, index `strip - 1`.
  Float_t gain_right[16]; ///< Right-end gain of strips 1-16.
  /// ADC subtracted from an end before its gain, only when that end fired: a
  /// short end's energy carries a fixed amount over the charge it collected
  /// whenever it triggers (see ChannelCal::offset_adc). Zero on the long
  /// ends and when no calibration is loaded.
  Float_t offset_left[16];
  Float_t offset_right[16];
  Float_t gain_strip0;  ///< Strip 0 gain.
  Float_t gain_strip17; ///< Strip 17 gain.
  Float_t gain_cathode; ///< Cathode gain.
  /// Grid gain: the reciprocal of its modal beam peak, so the grid reads 1.0
  /// for a beam event. Zero on files calibrated before it was written, where
  /// Decode() falls back to the full-scale normalisation.
  Float_t gain_grid;
  /// Per-strip multiplicative alignment, indexed by strip 0-17 and applied to
  /// every end of the strip after its gain. Identity (1.0) when no alignment
  /// is loaded, so an unaligned dataset decodes unchanged.
  Float_t strip_factor[18];
  /// Per-strip additive alignment in calibrated units, indexed by strip 0-17.
  /// Together with #strip_factor it puts the strip total on the scale with
  /// the beam at 1.0 and the two-particle pile-up at 2.0. It belongs to the
  /// strip, not an end, and is carried by the long end (L on odd strips, R on
  /// even) or the unsegmented value, so Total() stays the plain sum of the
  /// ends. Zero when no alignment is loaded.
  Float_t strip_offset[18];
  /// Whether a calibration tree was found. When true, the decoded values are in
  /// arbitrary calibrated units; when false they remain raw ADC.
  Bool_t is_normed;
  /// @}

  /// @name Decoded per-event values
  /// In arbitrary units when #is_normed, otherwise raw ADC. Each end already
  /// carries its strip's alignment factor, and the long end its offset, so a
  /// strip's deposit is the plain sum of its ends, which is what Total()
  /// returns.
  /// @{
  Double_t left[16];  ///< Left end of strips 1-16, index `strip - 1`.
  Double_t right[16]; ///< Right end of strips 1-16, index `strip - 1`.
  Double_t strip0;    ///< Strip 0.
  Double_t strip17;   ///< Strip 17.
  Double_t cathode;   ///< Cathode energy.
  Double_t grid;      ///< Grid energy.
  /// @}

  TTree *tree_;       ///< Events tree or chain bound by Attach(). Not owned.
  Int_t loaded_tree_; ///< Tree number whose gains are loaded; `-1` for none.

  /// @brief Construct detached, with zeroed values and no gains loaded.
  EnergyView()
      : strip0_adc(0), strip17_adc(0), cathode_adc(0), grid_adc(0),
        gain_strip0(0.0f), gain_strip17(0.0f), gain_cathode(0.0f),
        gain_grid(0.0f), is_normed(kFALSE), strip0(0.0), strip17(0.0),
        cathode(0.0), grid(0.0), tree_(nullptr), loaded_tree_(-1) {
    for (Int_t k = 0; k < 16; k++) {
      leftdE_adc[k] = 0;
      rightdE_adc[k] = 0;
      gain_left[k] = 0.0f;
      gain_right[k] = 0.0f;
      offset_left[k] = 0.0f;
      offset_right[k] = 0.0f;
      left[k] = 0.0;
      right[k] = 0.0;
    }
    for (Int_t s = 0; s < 18; s++) {
      strip_factor[s] = 1.0f;
      strip_offset[s] = 0.0f;
    }
  }

  /**
   * @brief Bind to an events tree and set up the branch addresses.
   * @param t Tree or chain to read. Borrowed; must outlive this view.
   * @return `kTRUE` if the expected branches were found.
   */
  Bool_t Attach(TTree *t);

  /**
   * @brief Decode the currently loaded entry into the value members.
   *
   * Applies the per-channel offsets (to ends that fired) and gains, drops the
   * short end when `IGNORE_SHORT_STRIPS` is set, then multiplies every end by
   * its strip's #strip_factor and adds #strip_offset to the long end when it
   * fired (a reading of 0 stays 0, as for the per-channel offsets).
   *
   * @note Reads whichever entry the bound tree has loaded, so call
   *       `GetEntry()` first.
   */
  void Decode();

  /// @brief Load the calibration gains for the current file.
  /// @note Called automatically when a chain rolls over to a new tree; only
  ///       needed directly when driving decoding by hand.
  void LoadGains();

  /// @brief A strip's deposit: the sum of its ends, or the unsegmented value.
  /// @param strip Anode strip, 0 to 17.
  Double_t Total(Int_t strip) const {
    if (strip <= 0)
      return strip0;
    if (strip >= 17)
      return strip17;
    return left[strip - 1] + right[strip - 1];
  }
  /// @brief The value on a beam-gate axis: a strip's deposit, or the grid
  ///        (GATE_AXIS_GRID) over its anchor.
  Double_t Axis(Int_t axis) const {
    return axis == GATE_AXIS_GRID ? grid : Total(axis);
  }
  /// @brief Every strip's deposit at once, for code that wants an array.
  /// @param[out] out 18 values, one per strip.
  void Totals(Double_t *out) const {
    for (Int_t s = 0; s < 18; s++)
      out[s] = Total(s);
  }
};

#endif
