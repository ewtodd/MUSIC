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
