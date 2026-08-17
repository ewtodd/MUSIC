# L/R gain-matching fix + automatic pileup removal — August 2026

## Summary

Two problems were fixed for 37Cl:

1. **Gain mismatch ("sawtooth")**: beam events looked fine, but the more
   non-beam-like an event was, the worse its per-strip totals alternated
   odd/even. The short side of every strip was calibrated at its *beam-mode*
   ADC reading instead of its *full-charge* response, so any event that put a
   different-than-beam fraction of its charge on the short end read the wrong
   total (up to 3–12× too big).
2. **Automatic noise/pileup removal**: pileup was only handled by CoMPASS
   flags (useless here — 0 flagged events) and hand-tuned energy thresholds.
   Two changes make it automatic: complete events with a multi-hit anode
   channel (the unambiguous signature of two ions in one coincidence window)
   are now rejected at the event-builder level, and event completeness now
   respects `IGNORE_STRIP_0` (which was silently discarding ~2/3 of otherwise
   good events for 37Cl).

---

## 1. L/R gain mismatch — the sawtooth

### Symptom

Per-strip totals of reaction events (and any event whose charge deposition
differs from the beam) alternate high/low between odd and even strips, and
the amplitude grows the further the event is from the beam. Beam events are
fine.

### Root cause

The two ends of each strip have different ADC scales: the "long" side preamp
reads ~1.5× the "short" side preamp, and the long side receives ~83% of the
charge (η ≈ 0.83 for this detector/geometry). A correct calibration must
satisfy charge conservation:

```
L / F_L + R / F_R = Q / Q_beam          (partition-independent total)
```

where `F_L`, `F_R` are the *full-charge* ADC responses of each end (what each
side reads when it carries the entire deposit). The old pass 1 of
`ComputeLRGainMatch` instead anchored the short side at its **beam-mode**
reading — the tiny value it reads at the beam partition — via a 2D
correlation ridge that locked onto the beam blob's top edge (odd strips) or
pileup (even strips):

| strip | old short anchor | true full-charge response |
|-------|------------------|---------------------------|
| 1     | 251 ADC          | ~1260 ADC (median fallback) |
| 13    | 307 ADC          | 1400 ADC (measured)        |
| 14    | 227 ADC          | 1310 ADC (measured)        |

The strip total with an over-gained short side scales as

```
total = (Q/Q_beam) · [η′/η + (1−η′)·(anchor error)]
```

so an event whose partition η′ differs from the beam's η reads 3–12× too big
on the affected strip. Pass 2 (eSum alignment) and the per-strip strip
factors absorb the error for *beam* events (consistent partitions), which is
exactly why beam looked okay and everything else looked worse the further it
got from beam. Odd/even alternation appeared because the short anchors
differed systematically between parities (odd strips ~250–330 ADC, even
~145–315 ADC).

Measured on run 100 (Events_Run100_1_c000.root, 337k events): mean strip
total vs short-side charge fraction rose from **1.88 → 4.56 a.u. (+140%)**
with the old calibration.

### The fix (`CalibrateBeamLRGain.cpp`, `CalibrateBeamInternal.hpp`)

Pass 1 now anchors the two sides as follows:

- **Long side** — unchanged: `gain_long = 1 / long_beam_peak`. Anchoring the
  long side at its beam peak keeps long-only events at 1.0 a.u., matching
  every downstream threshold (noise 0.85, beam gates) that assumes a flat
  beam.
- **Short side** — anchored at its **full-charge response**, measured from a
  *relative shoulder slice*: the mode of the short side among events where
  the long side reads low (`long < 0.35 × long_beam_peak`). Those events
  carry most of the charge on the short end, so the short side reads the same
  charge the long side reads at full charge. The mode is searched in
  `[0.40, 1.25] × long_beam_peak`: the lower edge sits above the short side's
  beam-mode pile, the upper edge below its pileup reading (2× full charge).
  A mode hugging either band edge is treated as unmeasured.
- **Fallback** — upstream strips (1–9 for 37Cl) have so few extreme-partition
  events that the slice is empty; they use the **same-parity median** of the
  measured full-charge anchors (short side = R on odd strips, L on even
  strips; each parity's preamps share one electronics scale), with a
  cross-parity last resort. Pass 2's eSum alignment then fine-tunes each
  strip individually with full beam statistics.

The old 2D-correlation-ridge anchor and the fixed-window `[270,430]` shoulder
slice were removed (both fail on this data; the fixed window was empty for
most upstream strips).

### Verification

After re-running `calibrate-beam` on the same subfile:

- Measured shoulder anchors land at the full-charge scale: strips 10–16 give
  **1014–1400 ADC**; upstream strips take the median fallbacks
  (odd 1263, even 1205 ADC).
- The L/R gain-match eSum reference drops from 1.59 → **1.11 a.u.** — the
  beam eSum now sits at ≈1.0–1.17 as designed, so the pileup/noise thresholds
  (1.4 / 0.85 / 2.0) are physical again.
- The partition test is now **flat: 1.77 → 1.85 a.u. (+5%)** instead of the
  old +140% rise. The remaining ±10–15% U-shape is the residual η-position
  dependence that the notebook addresses with a separate poly2 ridge
  correction (not implemented here — see Follow-ups).

---

## 2. Automatic noise/pileup removal

### Multi-hit rejection (`EventBuilder.cpp`, new `REJECT_MULTI_HIT_EVENTS` knob)

The event builder always counted per-channel hits in the window (the `Hits`
branch) but only used the per-channel *dedup* value — with
`DEDUP_STRATEGY=kLARGEST_ENERGY` a second ion in the 8 µs window made one
strip read 2× instead of flagging the event. At the 45 kHz rate that is
~2.5% of events.

`REJECT_MULTI_HIT_EVENTS` (Bool_t, default kFALSE, kTRUE for 37Cl) now rejects
complete events where any **anode** channel fired more than once in the
window. This is the definition of pileup — no energy thresholds to tune. The
hit counts remain in the tree regardless, so the selection is fully
reversible.

**What it catches, and what it can't.** The check catches *electronics-resolved*
double hits — the second ion's grid hit lost in grid dead-time, so its anode
hits join the previous event. On the new runs that is ~5% of events (visible
in the pipeline log as `Events with multi-hit on any anode: ...`; a dedicated
`Stored events (REJECT_MULTI_HIT_EVENTS=true): ...; N rejected as multi-hit
pileup` line reports the rejection rate). It does **not** catch the case the
shaper merges two pulses on a channel into one hit with summed amplitude
(count 1, no hit-count signature). Merged pileup is handled by the
energy-based filters downstream instead:

- two ions within one shaping time sit at nearly the same chamber position,
  so their tracks overlap on *every* strip and the event reads ~2×
  everywhere — caught by `IsPileup` (≥2 strips ≥ 1.4 a.u.) in
  strip-sum-scatter;
- a single strip reading 1.6–2.4× beam with no other high strip is ~0.11% of
  events and is *physically ambiguous* with a reaction jump (both raise one
  strip's charge); no automatic filter can separate them without killing
  reaction events — the analysis's reaction gates handle this band.
  `IsHighStrip` (> 2.0 a.u.) catches its upper part.
  An L/R-ratio test (single-channel pileup breaks the strip's partition
  ratio) was evaluated and rejected: the short side's beam level is so small
  that normal partition spread trips it (~12% false-positive rate).

### Completeness respects `IGNORE_STRIP_0` (`EventBuilder.cpp`)

`CheckEventComplete` required Strip0 unconditionally, but 37Cl sets
`IGNORE_STRIP_0=kTRUE` (its Strip0 rarely fires: 64–66% miss rate). The
result: only ~33% of events were stored as "complete", and the other ~2/3
were discarded even though every physics strip had fired. Completeness now
mirrors the downstream `AllStripsFired` predicate: Strip0 (and Strip17, when
ignored) only counts when the analysis actually uses it. This multiplies the
usable statistics by ~3 for 37Cl.

---

## 3. Supporting change — scatter cache invalidation

`strip-sum-scatter` caches its scatters/reservoir in
`StripSumScatter_cache.root` under a fingerprint of every knob that changes
its contents — but the calibration content was not part of it, so a
re-calibration with identical filter settings silently reused the old
cache. `BuildFingerprint` now folds per-run sums of `GainLeft`, `GainRight`
and `StripFactor` into the fingerprint (version bumped v12 → v13), so any
re-calibration automatically invalidates and rebuilds the cache.

---

## 4. Cluster-histogram skip flag + peak-like sanity gate (run-84 follow-up)

After re-running on run 84, two more problems surfaced and were fixed:

### 4a. `SKIP_CLUSTER_HISTS` flag

The interactive region-trace overlay always drew the nine per-class
cluster-variable histograms (`ClusterVarHists`: energy, peak3, plateau,
tail, reacstrip, mult, trigtaildev, reacslope3, beamdev — raw and
Savitzky-Golay passes, ~54 PNGs), which are only useful for the Python-side
blind clustering. New `StripSumScatterConfig::SKIP_CLUSTER_HISTS` (default
kFALSE, kTRUE for 37Cl) skips them entirely; the region-trace overlays are
still produced.

### 4b. Peak-like sanity gate in the beam-peak fitter

Run 84's Strip17 spectrum is a smear (mode ≈ 245 ADC, mean ≈ 1106, tail to
5000+), and the beam-peak fitter's narrow refit locked onto an arbitrary
sub-peak (fit μ = 81 ADC with σ = 932 — an 11σ-wide "peak"). The resulting
gain made Strip17 totals read 0–38 a.u. (70% of events above 2 a.u.), which
inflated the cluster trigger's pooled beam sigma to ~1.6 a.u.; the 5σ
reaction-onset gate then sat at ~8 a.u., so **no event ever triggered** and
the trigger-centered cluster variables (`trigtaildev`, `reacslope3`) came
out empty, tripping ROOT's canvas-range errors.

Fix: `FitBeamPeakGaussian` now rejects any fit whose width exceeds half its
centroid (`σ > 0.5·μ` → uncalibrated), and `ReduceToAnchors` applies the
same gate to the robust-mode fallback. A no-peak channel reads 0 a.u. — the
same net effect as an ignored/no-fire guard (37Cl ignores Strip17 anyway).
After the fix on run 84: beam sigma back to 0.091 a.u., trigger gate 0.46
a.u., and reaction events trigger again ((α,α′): 71%, (α,n): 54%).

## Files changed

| File | Change |
|------|--------|
| `tooling/include/CalibrateBeamInternal.hpp` | `StripPairSamples::shoulder` → `slice_short`/`slice_long` (relative shoulder slice) |
| `tooling/src/CalibrateBeamLRGain.cpp` | Short-side anchor = full-charge response from the relative shoulder slice + same-parity median fallback; removed 2D-correlation ridge and fixed `[270,430]` window |
| `tooling/src/EventBuilder.cpp` | Completeness respects `IGNORE_STRIP_0/17`; `REJECT_MULTI_HIT_EVENTS` rejection of multi-hit anode events (+ visible rejection log line) |
| `tooling/include/Constants.hpp` | New `REJECT_MULTI_HIT_EVENTS` knob |
| `tooling/src/Constants.cpp` | Default `REJECT_MULTI_HIT_EVENTS = kFALSE`; new `SKIP_CLUSTER_HISTS` default |
| `tooling/include/Constants.hpp` | New `SKIP_CLUSTER_HISTS` knob |
| `tooling/src/StripSumScatter.cpp` | Fingerprint v13 → v14 includes per-run calibration sums (GainLeft/GainRight/StripSlope/StripIntercept) |
| `tooling/src/StripSumScatterScatter.cpp` | `SKIP_CLUSTER_HISTS` skips `ClusterVarHists` |
| `tooling/src/CalibrateBeamFits.cpp` | Peak-like sanity gate in `FitBeamPeakGaussian` (σ ≤ 0.5·μ, else uncalibrated) |
| `tooling/include/CalibrateBeam.hpp` | `StripAlignmentResult`: `factors` → `slopes`/`intercepts` + `pileup_centroids` |
| `tooling/src/CalibrateBeam.cpp` | Two-point linear normalization: double-beam-gated pileup centroid per strip, slope/intercept derivation; calibration tree writes `StripSlope`/`StripIntercept` |
| `tooling/include/Normalization.hpp`, `tooling/src/Normalization.cpp` | `EnergyView` applies `total = slope·total + intercept` (legacy `StripFactor` read path kept) |
| `analysis/37Cl/config/Constants.cpp` | `REJECT_MULTI_HIT_EVENTS = kTRUE`, `SKIP_CLUSTER_HISTS = kTRUE` (plus the user's test-run edits: RUNS 84/85, candidate strip, post-trigger strips) |

## How to re-run and validate

```sh
nix build .#37Cl          # or: make DATASET=37Cl
./result/bin/pipeline     # rebuilds events + calibration for the configured runs
./result/bin/strip-sum-scatter   # cache auto-rebuilds (fingerprint changed)
```

Watch for in the pipeline log:

- `short_anchor=... ADC (shoulder, n=...)` values in the ~1000–1400 ADC range
  (previously ~150–700).
- `L/R gain match: eSum reference (median) = ~1.1` (previously ~1.6).
- `[uncalibrated no-peak] ... gain 0` lines for channels with no real beam
  peak (e.g. a noisy guard strip) instead of a garbage anchor.
- `strip N ... pileup=... slope=... intercept=...` lines in the alignment
  pass: pileup centroids ~2.0–2.2 a.u., slopes ~0.87–0.95 (run 84).
- `Complete events:` fraction jumps from ~33% to ~90%+ (Strip0 no longer
  required).
- `Events with multi-hit on any anode:` still counted, but those events are
  no longer stored — the new `Stored events (REJECT_MULTI_HIT_EVENTS=true):
  ...; N rejected as multi-hit pileup` line gives the rejection count.

## 5. Two-point linear normalization (sublinearity correction)

The ADC response is slightly sublinear: measured on run 84, the double-beam
(merged-pileup) population reads **2.04–2.17 a.u. where 2× the beam level
predicts 2.27** — ~5–10% compression at 2× charge, larger than the <2% the
notebook assumed. A single multiplicative factor cannot fix this (it maps
the beam onto 1.0 and leaves the pileup wherever it lands).

`FindStripCentroidAlignment` now derives a per-strip **two-point linear
correction** instead of the multiplicative factor:

- beam centroid B per strip — the peak nearest 1.0 in the eSum projection
  (unchanged machinery; the long-only population defines the 1.0 level, so
  every existing threshold stays valid);
- pileup centroid P per strip — the mode of the **double-beam-gated**
  population: a second 2D histogram filled only from events whose total
  energy (sum over strips 1–16) sits in [1.7, 2.4]× the beam-total mode.
  The gate excludes reaction events (they add only a fraction to the total),
  and only both-ends events feed the mode (the long-only pileup population
  is missing the short side's charge and would bias the anchor low). The
  mode is smoothed and refined with a narrow Gaussian fit for sub-bin
  precision;
- slope = 1/(P − B), intercept = 1 − slope·B, so the line passes through
  (B, 1.0) and (P, 2.0). Strips whose pileup peak is unmeasurable (upstream
  strips with few both-ends pileup events, or an unphysical slope outside
  [0.8, 1.25]) fall back to the **median slope** of the measured strips,
  anchored on their own beam centroid.

**Why the corrected pileup sits ~1–2% below 2.0, and why that is
deliberate.** The Gaussian centroid is the anchor rather than the raw
max-bin mode because it is *stable*: the raw mode jumps by a histogram bin
with per-subfile statistics, which made the per-subfile calibrations
inconsistent across a run (tried and rejected). The Gaussian centroid sits
systematically above the skewed peak's mode (the partial-pileup tail on the
high side), so the corrected peak lands consistently ~1.5% below 2.0 — the
same small offset in every subfile, which is what the downstream analysis
sees. The line itself maps the measured pileup centroid onto exactly 2.0 by
construction.

The calibration tree now stores `StripSlope[18]/F` and
`StripIntercept[18]/F` (replacing `StripFactor`); `EnergyView::Decode`
applies `total = slope·total + intercept` after the per-channel gains (with
a legacy read path for old files carrying `StripFactor`). The scatter-cache
fingerprint (v14) folds the new branch sums in.

Verified on run 84 after re-calibration: long-only beam modes 0.985–1.005,
double-beam-gated pileup modes 1.965–2.005 (strips 2–16) — consistently
~1.5% below 2.0, the same offset in every subfile (see above for why that
is deliberate). The ADC's raw pileup reading (1.83–1.93× beam, i.e. ~5–10%
sublinear at 2× charge) is physics; the fit maps it onto ~2.0.

## Known limitations / follow-ups

- **Residual η-position dependence**: with the long side anchored at its beam
  peak, the total still carries ±10–15% partition dependence for extreme
  partitions (the notebook's "linear" scheme has the same). A full η
  correction (poly2 ridge on the (L,R) correlation, then per-event total
  reweighting) would remove it; it is deliberately out of scope here.
- **Upstream short anchors are imputed**: strips 1–9 rarely show the short
  side at full charge, so their short gain comes from the same-parity median.
  If the per-strip preamp spread ever grows, pass 2's eSum alignment is the
  safety net that re-levels those strips.
- **Sublinearity beyond 2×**: the two-point line is exact at beam (1.0) and
  pileup (2.0) and interpolates/extrapolates linearly between. The
  triple-charge region (>3×) may still be off by a few %; the existing
  `PILEUP_THRESHOLD`/`HIGH_STRIP_THRESHOLD` filters reject it anyway.
- **Upstream short anchors are imputed**: strips 1–9 rarely show the short
  side at full charge, so their short gain comes from the same-parity median.
  If the per-strip preamp spread ever grows, pass 2's eSum alignment is the
  safety net that re-levels those strips.
