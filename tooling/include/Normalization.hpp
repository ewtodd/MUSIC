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
 * tree was found — see #is_normed and Unit().
 *
 * @note Holds a non-owning pointer to the tree it is attached to, which must
 *       outlive it. Gains are reloaded automatically when a `TChain` rolls over
 *       to a new file, tracked by #loaded_tree_.
 */
struct EnergyView {
  UShort_t left_0_17_adc[18]; ///< Raw left-side strip ADC values.
  UShort_t rightdE_adc[18];   ///< Raw right-side strip ADC values.
  UShort_t hits_adc[36];      ///< Raw per-slot hit ADC values.
  Short_t cathode_adc;        ///< Raw cathode ADC value.
  Short_t grid_adc;           ///< Raw Frisch grid ADC value.

  Float_t gain_left[18];  ///< Per-strip left-side gain.
  Float_t gain_right[18]; ///< Per-strip right-side gain.
  Float_t gain_cathode;   ///< Cathode gain.
  /// Whether a calibration tree was found. When true, the decoded values are in
  /// arbitrary calibrated units; when false they remain raw ADC.
  Bool_t is_normed;

  /// Per-strip multiplicative alignment, `pol3_reference / centroid`, applied
  /// to #total after the per-channel gain. Identity (1.0) when no alignment is
  /// loaded, so an unaligned dataset decodes unchanged.
  Float_t strip_factor[18];

  /// @name Decoded per-event values
  /// In arbitrary units when #is_normed, otherwise raw ADC.
  /// @{
  Double_t left[18];  ///< Left-side energy per strip.
  Double_t right[18]; ///< Right-side energy per strip.
  Double_t total[18]; ///< Summed energy per strip, after #strip_factor.
  Double_t cathode;   ///< Cathode energy.
  Double_t grid;      ///< Grid energy.
  /// @}

  TTree *tree_;       ///< Events tree or chain bound by Attach(). Not owned.
  Int_t loaded_tree_; ///< Tree number whose gains are loaded; `-1` for none.

  /// @brief Construct detached, with zeroed values and no gains loaded.
  EnergyView()
      : cathode_adc(0), grid_adc(0), gain_cathode(0.0f), is_normed(kFALSE),
        cathode(0.0), grid(0.0), tree_(nullptr), loaded_tree_(-1) {}

  /**
   * @brief Bind to an events tree and set up the branch addresses.
   * @param t Tree or chain to read. Borrowed; must outlive this view.
   * @return `kTRUE` if the expected branches were found.
   */
  Bool_t Attach(TTree *t);

  /**
   * @brief Decode the currently loaded entry into the value members.
   *
   * Applies the per-channel gains, then #strip_factor to #total.
   *
   * @note Reads whichever entry the bound tree has loaded, so call
   *       `GetEntry()` first.
   */
  void Decode();

  /// @brief Load the calibration gains for the current file.
  /// @note Called automatically when a chain rolls over to a new tree; only
  ///       needed directly when driving decoding by hand.
  void LoadGains();

  /// @brief Unit label for the decoded values, for axis titles.
  /// @return Arbitrary units when #is_normed, otherwise an ADC label.
  const char *Unit() const;
};

#endif
