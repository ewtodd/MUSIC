# VLM fine-tuning: current result and next steps

## Current result

The 37Cl simulator corpus now contains 10,000 events for each reaction strip
2--8 for `(alpha,n)`, `(alpha,alpha')`, and `(alpha,p)`, plus 10,000 beam
events. `(alpha,p)` is currently the fourth `other` class.

Gemma 4 E2B was trained with 4-bit QLoRA. The final 2,000-example-per-class
run used 6,745 training events and 1,255 validation events, with reaction strip
8 held out. Its final validation loss was 0.626. Recall was 0.181 for
`(alpha,n)`, 0.761 for `(alpha,alpha')`, 0.993 for beam, and 0.705 for other.
The preceding 500-example-per-class run had a worse total loss of 0.795 but
substantially better `(alpha,n)` recall of 0.651. The final-epoch checkpoint
therefore is not the right model-selection rule for the physics objective.

On the 6,578-event experimental compute-regions reservoir, the final adapter
ranked every existing tagged event above every beam-flat event (AUC 1.000) and
produced no beam-flat false positives at `p(an) >= 0.5`. However, only 3.4% of
existing tagged events passed that threshold. Acceptance fell from 6.2% at
reaction strip 2 to zero at strips 7--9. The model has learned useful
beam-versus-reaction ranking, but not a calibrated or strip-transferable
`(alpha,n)` tag.

The current adapter is generated and ignored by Git at:

```text
analysis/37Cl/models/vlm-gemma-e2b-lora
```

## Recommended order of work

1. Add threshold sweeps and ROC/precision-recall reporting on a complete
   beam-gated experimental subfile. Do not choose 0.5 merely because it is the
   conventional probability threshold: this is a softmax over four selected
   letter tokens and is not calibrated. The 400 beam-flat cache events are
   useful but unusually easy negatives, so they are not sufficient alone.

2. Save an adapter after every epoch and choose the checkpoint by `(alpha,n)`
   recall or F1 subject to a maximum beam false-positive rate. Also report
   macro-F1 and the full confusion matrix. Do not select by total validation
   cross-entropy: the completed runs show that it favors the wrong model.

3. Add horizontal translation augmentation. The held-out-strip failure shows
   that the model is using absolute reaction position. Shift complete simulated
   traces left and right during training, then perform leave-one-strip-out
   validation for every reaction strip. Consider a reaction-relative crop if
   translation augmentation is insufficient.

4. Improve the `other` class. `(alpha,p)` is a useful start but does not cover
   pileup, incomplete tracks, off-beam particles, detector excursions, and
   other experimental backgrounds. Add simulated examples where possible and
   carefully selected real negatives where simulation is inadequate.

5. Quantify simulation-to-data agreement per strip before generating more
   events: beam mean and width, plateau height, downstream collapse, transition
   shape, and end-strip behavior. More samples from a mismatched simulator will
   reinforce the mismatch rather than fix it.

6. Train XGBoost, a small MLP, and a compact 1D CNN on exactly the same traces
   and leave-one-strip-out splits. The input is only 18 numerical values. Gemma
   should be retained only if it transfers better to experimental data or
   achieves comparable performance with materially less supervision.

## Suggested next experiment

- Use 500, 1,000, and 2,000 examples per class.
- Add horizontal shifts and per-epoch checkpoints.
- Select by `(alpha,n)` recall with a bounded beam false-positive rate.
- Evaluate every held-out reaction strip, not only strip 8.
- Compare the selected Gemma adapter directly with XGBoost on identical splits.
- Apply both to one complete beam-gated experimental subfile and inspect the
  threshold sweep before attempting a cross section.

Do not spend time generating a larger simulator corpus until these checks show
that data volume, rather than strip shortcuts or simulation-domain mismatch, is
the limiting factor.

## Follow-up completed after the first study

The branch was rebased onto `origin/main` and the first three recommendations
were implemented. Gemma training now translates deviations from the beam trace
by up to three strips, saves every epoch, and selects by `(alpha,n)` acceptance
under a 1% aggregate simulated-background leakage constraint.

The augmented 500-example-per-class Gemma run selected epoch 3. On held-out
strip 8 its argmax recalls were 0.683 `(alpha,n)`, 0.562 `(alpha,alpha')`, 0.927
beam, and 0.338 other. At the constrained threshold `p(an) >= 0.164`, its
`(alpha,n)` acceptance was 0.952 with zero simulated-beam leakage under the old
beam-only threshold criterion. Re-evaluate this checkpoint with the corrected
aggregate-background criterion before using that threshold.

Two identical-split comparisons were added:

- XGBoost, trained conventionally on 2,000 examples per class, had held-out
  `(alpha,n)` argmax recall of 0.987, 0.966, 0.915, 0.922, 0.839, 0.741, and
  0.422 for strips 2 through 8. Beam recall remained 0.993 or better.
- Google TabFM 1.0.0, used zero-shot through in-context examples rather than
  weight training, had held-out `(alpha,n)` recall of 0.878, 0.922, 0.924,
  0.920, 0.937, 0.663, and 0.365 for strips 2 through 8 with 500 examples per
  class, 100 rows per context, and eight estimators. Beam recall was 0.991 or
  better. On strip 8, the corrected 1% aggregate-background operating point
  gave 0.492 `(alpha,n)` acceptance and 0.4% background leakage.

TabFM is the best match to the original goal of avoiding task-specific weight
training. It nearly matches XGBoost through strip 6 and clearly beats the
fine-tuned VLM as a no-training method. Both tabular methods still degrade at
strips 7--8, confirming that downstream simulation/domain transfer is the next
physics problem rather than model scale.

Next, apply TabFM and XGBoost to the complete beam-gated experimental subfile.
Use simulated labels only to choose an operating threshold, then inspect the
experimental selected traces and per-strip rates. The stale-Pandas-cache segmentation fault in the experimental loader was fixed
by replacing pickle-backed scalar loading with direct ROOT scalar buffers and a
NumPy archive. Do not pass the complete 106,110-event reservoir to TabFM in one
call: that exhausted GPU resources badly enough to reset the NVIDIA kernel
module. `vlm_tabfm.py --experimental` now defaults to a 1,000-row probe,
calls TabFM in independent 32-row chunks, and checkpoints each requested row range so an
interrupted run can resume; increase those limits cautiously.

A safe 6,000-event experimental test used four estimators, 100 context rows,
and 16-row inference chunks. Two events (0.033%) passed the strip-8 simulated
operating point `p(an) >= 0.5116`; 59 (0.98%) exceeded 0.01. The median score
was about `1e-6`, the 99th percentile 0.0091, and the maximum 0.557. This is
consistent with a rare-reaction reservoir and confirms that bounded TabFM
inference works, but these sequential events did not overlap the existing
compute-regions cache by `SeedTs`, so this sample does not yet measure event
identity agreement.

In parallel, compare simulation and data distributions at strips 7--8 before
trusting either model for a cross section. The TabFM pretrained weights are
non-commercial and may only be used under their
`tabfm-non-commercial-v1.0` license.
