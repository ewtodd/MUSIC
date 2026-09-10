"""Staged checks on the VLM path, cheapest first.

Run this before kicking off anything long. The stages are ordered by what
they cost, and the script stops at the first failure, so a broken chat
template surfaces in seconds instead of an hour into a classification run:

    stage 1  rasterization on real events       no network, no model
    stage 2  processor and chat template        processor files only (MBs)
    stage 3  one real forward pass on 8 events  the checkpoint (GBs)

Stage 2 is the one that matters most. The batched multimodal chat template
-- a list of conversations, left padding, one image apiece -- is the seam
where the API is most likely to differ from what this code assumes, and it
needs no weights to exercise. Stage 3 then confirms only that the answer
slot really is the last position and that the four class letters come back
as a usable distribution.

    python vlm_selftest.py            # all stages
    python vlm_selftest.py 1 2        # only these

Stage 1 also writes a contact sheet to the plots dir so the images can be
looked at rather than trusted, next to the region_traces_reac<r> figure the
prompt uses as its worked example. If a person cannot tell (a,n) from beam
in that PNG, no model will either, and that is a rendering bug rather than a
model result.

Everything here runs on real beam-gated events, capped to
VLM_SELFTEST_FILES subfiles so it stays quick. There is deliberately no
labelled truth set: the comparison this whole path is built for is against
the published analysis, not against labels invented locally.
"""

import sys

import numpy as np

import blind_an
import config
import vlm
import vlm_apply

_PROBE_EVENTS = 8
_BEAM_RMS = [0.0]  # filled by stage 1


def _fail(msg):
    raise AssertionError(msg)


def _reservoir():
    """A small beam-gated slice of real events, with its beam reference.

    No reaction-strip selection anywhere in this file. blind_an's CFD gates
    on a peak excess above BLIND_REAC_ONSET_NSIGMA * beam RMS -- 5 * 0.0346
    = 0.173 on run 16 -- while the (a,n) plateau sits about 0.10-0.15 above
    beam, so the threshold is above the signal and what survives it is
    pileup. Sampling probe events through it would show the model, and the
    contact sheet, the wrong events.
    """
    X, _both, _ts = vlm_apply.beam_gated_reservoir(
        max_files=config.VLM_SELFTEST_FILES)
    beam_ref, _sigma = blind_an._beam_reference(X)
    if beam_ref is None:
        _fail("no pure-beam events in the slice; cannot build a reference")
    return X, beam_ref


def stage1():
    """Rasterization on real events: the window, jitter, determinism, and a
    contact sheet to put beside region_traces_reac<r>.png."""
    print("stage 1: rasterization on real events")
    X, beam_ref = _reservoir()
    print("  beam reference " + " ".join(f"{v:.2f}" for v in beam_ref))

    # The window has to contain the traces. These are the values
    # StripSumScatter::DrawRegionTraces frames region_traces_reac<r> with, so
    # anything leaving it is leaving the view a person classifies from too.
    lo, hi = config.VLM_DE_MIN, config.VLM_DE_MAX
    inside = float(((X >= lo) & (X <= hi)).mean())
    print(f"  {100.0 * inside:.1f}% of strip samples fall inside "
          f"[{lo}, {hi}]")
    if inside < 0.90:
        _fail(f"only {100.0 * inside:.0f}% of the data is on canvas; the "
              "render window is wrong")

    # A plain random draw from the gated reservoir -- no selection at all, so
    # the sheet shows what the model will actually be handed.
    rng = np.random.default_rng(config.SEED)
    per_row = config.VLM_SHEET_PER_ROW
    n_show = min(per_row * config.VLM_SHEET_ROWS, X.shape[0])
    pick = rng.choice(X.shape[0], n_show, replace=False)
    img = vlm.render_traces(X[pick], beam_ref)
    if img.dtype != np.uint8 or img.shape[1:] != (config.VLM_IMG_H,
                                                  config.VLM_IMG_W, 3):
        _fail(f"render_traces returned {img.shape} {img.dtype}")
    ink = (img != np.asarray(config.VLM_BG_RGB, np.uint8)).any(axis=3)
    frac = ink.mean(axis=(1, 2))
    if frac.min() <= 0.0:
        _fail("an event rendered as a blank canvas")
    print(f"  {n_show} random gated events, ink fraction "
          f"{frac.min():.3f}..{frac.max():.3f}")
    rows = [np.concatenate(list(img[i:i + per_row]), axis=1)
            for i in range(0, n_show - per_row + 1, per_row)]

    probe = X[pick[:_PROBE_EVENTS]]
    a = vlm.render_traces(probe, beam_ref)
    if not np.array_equal(a, vlm.render_traces(probe, beam_ref)):
        _fail("render_traces is not deterministic")
    jit = vlm.render_traces(probe, beam_ref, dy=1.0, dhalf=0.5)
    if np.array_equal(a, jit):
        _fail("dy/dhalf changed nothing; the seed ensemble would report a "
              "systematic of exactly zero")
    print("  jitter moves "
          f"{100.0 * (a != jit).any(axis=(1, 2, 3)).mean():.0f}% of probes")

    out = config.PLOTS_DIR / config.VLM_PLOT_SUBDIR
    out.mkdir(parents=True, exist_ok=True)
    path = out / "vlm_selftest_sheet.png"
    from PIL import Image
    sheet = np.concatenate(rows, axis=0)
    Image.fromarray(sheet).resize((sheet.shape[1] * 3, sheet.shape[0] * 3),
                                  Image.NEAREST).save(path)
    print(f"  contact sheet -> {path}")
    print("  stage 1 OK")


def stage2():
    """Processor and chat template -- no weights downloaded."""
    print("stage 2: processor + chat template")
    builder = vlm.PromptBuilder()
    builder.report()
    print(f"  prompt ({len(builder.prompt)} chars):")
    for line in builder.prompt.splitlines():
        print(f"    | {line}")

    # The budget argument is undocumented, and setting attributes on the
    # image processor demonstrably does nothing. Dump what the processor
    # actually exposes so the name can be read off rather than guessed, then
    # put it in config.VLM_BUDGET_KWARG.
    if builder.n_image_tokens != builder.budget:
        import inspect
        ip = getattr(builder.processor, "image_processor", None)
        print(f"  image processor: {type(ip).__name__}")
        keys = sorted(getattr(ip, "to_dict", dict)().keys()) if ip else []
        hits = [k for k in keys
                if any(w in k.lower()
                       for w in ("token", "budget", "patch", "size", "pixel"))]
        print(f"    config keys mentioning token/budget/patch/size/pixel:")
        for k in hits:
            print(f"      {k} = {getattr(ip, k, '?')}")
        for obj, what in ((ip, "image_processor.__call__"),
                          (builder.processor, "processor.__call__")):
            if obj is None:
                continue
            try:
                params = list(inspect.signature(obj.__call__).parameters)
            except (TypeError, ValueError):
                continue
            print(f"    {what} params: {params}")
        print("    -> set config.VLM_BUDGET_KWARG to whichever of these "
              "selects the soft-token count")

    tok = builder.processor.tokenizer
    print(f"  class tokens {config.VLM_CLASS_TOKENS} -> "
          f"{builder.class_token_ids}")
    side = getattr(tok, "padding_side", None)
    if side != "left":
        _fail(f"tokenizer padding_side is {side!r}, not 'left'; with right "
              "padding logits[:, -1] is a pad token for the short rows")

    images = np.zeros((3, config.VLM_IMG_H, config.VLM_IMG_W, 3), np.uint8)
    images[1] = 255
    inputs = builder.build_inputs(images)
    keys = sorted(inputs.keys())
    print(f"  build_inputs keys: {keys}")
    if "input_ids" not in inputs:
        _fail(f"no input_ids in processor output (got {keys})")
    ids = inputs["input_ids"]
    if ids.shape[0] != 3:
        _fail(f"batched template collapsed {3} conversations to "
              f"{ids.shape[0]} rows -- apply_chat_template is not batching")
    print(f"  input_ids {tuple(ids.shape)}, "
          f"{[k for k in keys if 'pixel' in k or 'image' in k]}")

    tid = builder.image_token_id()
    if tid is None:
        print("  NOTE: no image token id on the processor; the budget could "
              "not be verified without a model (stage 3 retries it)")
    else:
        per_row = [int((row == tid).sum()) for row in ids]
        print(f"  image tokens per row: {per_row}")
        if len(set(per_row)) != 1:
            _fail(f"rows disagree on image token count: {per_row}")
        if per_row[0] == 0:
            _fail("no image tokens in the prompt; the image never reached "
                  "the template")

    # Every row must end on the same token: that is the answer slot the
    # readout indexes with logits[:, -1].
    last = {int(row[-1]) for row in ids}
    if len(last) != 1:
        _fail(f"rows end on different tokens {last}; logits[:, -1] would be "
              "reading a different position per row")
    print(f"  all rows end on token {last.pop()} (the answer slot)")

    # A permuted builder must produce a different prompt but the same letters.
    n_cls = len(config.VLM_CLASSES)
    other = vlm.PromptBuilder(letter_order=tuple(reversed(range(n_cls))))
    if other.prompt == builder.prompt:
        _fail("letter permutation did not change the prompt")
    if other.class_token_ids != builder.class_token_ids:
        _fail("class token ids must stay in LETTER order under permutation")
    print("  letter permutation changes the prompt, not the token ids")
    print("  stage 2 OK")


def stage3():
    """One real forward pass: the answer slot and the class distribution."""
    print("stage 3: forward pass")
    X, beam_ref = _reservoir()
    rng = np.random.default_rng(config.SEED)
    pick = rng.choice(X.shape[0], _PROBE_EVENTS, replace=False)

    clf = vlm.EventClassifier()
    print(f"  logits kwarg: {clf._logits_kwarg}")
    print(f"  reference figure: {clf.builder.reference is not None} "
          f"(strip {clf.builder.reference_strip})")
    probs, labels = clf.classify(X[pick], beam_ref, batch=_PROBE_EVENTS)
    if probs.shape != (pick.size, len(config.VLM_CLASSES)):
        _fail(f"class_probs returned {probs.shape}")
    if not np.allclose(probs.sum(axis=1), 1.0, atol=1e-4):
        _fail("posteriors do not sum to 1")
    if not np.isfinite(probs).all():
        _fail("non-finite posteriors")
    for i, name in enumerate(config.VLM_CLASSES):
        print(f"    {name:>6}: mean p = {probs[:, i].mean():.3f}")
    # Only p(an) decides anything. These are a RANDOM draw from the gated
    # reservoir, which is overwhelmingly unreacted beam, so p(an) should be
    # LOW here -- how the rest of the mass splits between beam and (a,a') is
    # not a result. The (a,n) rate at threshold is what a yield is built on,
    # so a high rate on a random draw is the failure that matters.
    pan = vlm.an_probability(probs)
    tagged = vlm.is_an(probs)
    print(f"  p(an) on a random gated draw: mean {pan.mean():.3f}, "
          f"max {pan.max():.3f}")
    print(f"  tagged (a,n) at p>={config.VLM_AN_THRESHOLD}: "
          f"{int(tagged.sum())}/{tagged.size}")
    if pan.mean() > 0.5:
        _fail(f"mean p(an) is {pan.mean():.2f} on a draw that is almost all "
              "unreacted beam; every yield built on this would be beam")
    if float(probs.std(axis=0).max()) < 0.01:
        print("  WARNING: posteriors barely vary across events -- the answer "
              "may not depend on the image.")
    print("  stage 3 OK")


_STAGES = {1: stage1, 2: stage2, 3: stage3}


def main():
    want = [int(a) for a in sys.argv[1:]] or sorted(_STAGES)
    print(f"music-ml: dataset {config.DATASET}, VLM self-test "
          f"({config.VLM_MODEL}), stages {want}")
    for n in want:
        if n not in _STAGES:
            raise SystemExit(f"no stage {n}; have {sorted(_STAGES)}")
        _STAGES[n]()
    print("all requested stages passed")


if __name__ == "__main__":
    main()
