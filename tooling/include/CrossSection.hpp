#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP
/**
 * @file CrossSection.hpp
 * @brief Absolute cross section per reaction strip, for every declared channel.
 *
 * A read-only pass over the scatter cache and the region cuts, covering every
 * channel in `CrossSectionConfig::CHANNELS`. For each channel and strip:
 *
 *     sigma = N_reac / (N_beam * n_gas * L_strip)
 *
 * where `N_reac` is the number of reactions of that channel at that strip, and
 * `N_beam` the number of beam particles that reached it **under every cut a
 * reaction there also had to pass** — so those cuts cancel in the ratio rather
 * than needing their own efficiency.
 *
 * `N_reac` comes from one of two places:
 *
 * - When the tag-efficiency store holds a record for the channel and strip, it
 *   is that record's count, unfolded by its efficiency and by the migration
 *   from the strip before. See TagEfficiency.hpp for the contract.
 * - Otherwise it is the count the mixture fit attributed to the reaction
 *   component, uncorrected, with the fit-versus-core disagreement taken as the
 *   region systematic.
 *
 * Beam energies come from the simulated unreacted beam, calibrated against the
 * measured per-strip energy loss, so each strip's centre-of-mass energy is the
 * simulation's rather than a nominal dE/dx table's — see BeamEnergies. Each
 * channel's TALYS curve is the sum of the residual channels its exit list
 * names, and its shape sets that channel's effective energies.
 */
#include <Rtypes.h>
#include <TString.h>
#include <map>
#include <utility>
#include <vector>

class TCutG;
class TFile;
class TGraph;
class TH2F;
struct CrossSectionChannel;

/**
 * @brief Computes and reports the absolute cross sections.
 *
 * Read-only with respect to the analysis products it consumes: it evaluates the
 * current cache, cuts and efficiency store, and changes none of them.
 */
class CrossSection {
public:
  /**
   * @brief Run every channel and write the tables and figures.
   *
   * Produces the table, the comparison against published values, one figure per
   * channel, and — when more than one channel is configured — the combined
   * figure.
   *
   * @return `kFALSE` when a prerequisite is missing, with the reason printed.
   */
  Bool_t Run();

  /**
   * @brief Residual nucleus left by a named exit channel.
   *
   * For an alpha on a beam nucleus, given an exit named as in
   * `CrossSectionChannel::talys_exits`.
   *
   * @param exit   Exit channel name.
   * @param z_beam Beam proton number.
   * @param a_beam Beam mass number.
   * @param[out] z Residue proton number.
   * @param[out] a Residue mass number.
   *
   * @return `kFALSE` when the name does not parse, leaving the outputs unset.
   */
  static Bool_t ExitResidue(const TString &exit, Int_t z_beam, Int_t a_beam,
                            Int_t &z, Int_t &a);
  /// @brief Display label for a channel.
  /// @param ch Channel to label.
  /// @return Its configured label, or one derived from its exit list when it
  ///         has none.
  static TString Label(const CrossSectionChannel &ch);

private:
  struct TalysCurve {
    TString label;
    TGraph *axn;
  };
  struct Point {
    Int_t reac;
    Double_t e_in, e_out, e_eff, e_eff_lo, e_eff_hi; // E_cm [MeV]
    Double_t n_reac, n_denom;
    Double_t sigma, stat, sys; // [mb]
    Double_t Err() const;
  };
  struct ChannelResult {
    const CrossSectionChannel *ch;
    std::vector<TalysCurve> talys;
    std::vector<Point> points;
  };

  Bool_t LoadCache();
  Bool_t LoadBeam();
  // Every residual-production graph of every model in
  // root_files/talys/talys_xs.root, keyed by model then (Z, A).
  void LoadTalys();
  // The channel's curves: per model, the sum over its exits' residues.
  std::vector<TalysCurve> ChannelCurves(const CrossSectionChannel &ch) const;
  Bool_t RunChannel(const CrossSectionChannel &ch, ChannelResult &out);
  // One strip's point, or kFALSE (with a printed reason) when it has none.
  Bool_t Strip(const CrossSectionChannel &ch,
               const std::vector<TalysCurve> &talys, Int_t reac, Point &pt);
  void CompareReference(const ChannelResult &r) const;
  // The figure for these channels; name is the file's basename.
  void Draw(const std::vector<const ChannelResult *> &rs,
            const TString &name) const;

  static Long64_t ReadCount(TFile &f, const char *name, Bool_t &ok);
  static TGraph *Clipped(TGraph *g, Double_t e_lo, Double_t e_hi);
  static Double_t EffectiveEnergy(TGraph *axn, Double_t e_out, Double_t e_in);
  static TCutG *ScaledCut(TCutG *cut, Double_t scale);
  static Double_t CountInCut(TH2F *scatter, TCutG *cut, Double_t scale);
  static Double_t Enclosed(Double_t nsigma);

  TFile *cache_ = nullptr;
  Long64_t n_seen_ = 0, n_beam_ = 0;
  Double_t n_gas_ = 0.0, areal_ = 0.0;
  Double_t dE_[18], e_strip0_ = 0.0, cm_frac_ = 0.0, e_mid_[17];
  std::vector<TString> talys_labels_;
  std::vector<std::map<std::pair<Int_t, Int_t>, TGraph *>> talys_raw_;
  // Unfolding state carried from one strip to the next within a channel.
  Int_t prev_reac_ = -1;
  Double_t prev_true_ = 0.0, prev_migrate_ = 0.0;
};

#endif
