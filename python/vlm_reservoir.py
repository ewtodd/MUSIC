"""Tier A of the event-identity check: the VLM on the compute-regions
reservoir alone.

compute-regions keeps every event its tag fired on (plus a sample of
beam-flat events it took its baseline from) in root_files/
StripSumScatter_cache.root, tree `traces`: calibrated per-end deposits in
the same beam-equals-one units the Python trace view uses, the reaction-strip
bit mask, the beam-flat flag and SeedTs. So the cheapest possible question
-- does the model say (a,n) on the events the tag says (a,n) on, and beam on
the ones the tag says are beam -- needs no events-file loading, no beam
gate and a few minutes of GPU: the reservoir is ~1500 events.

What this does NOT test is the other half of "the same events": whether
the model tags events the tag did not. That is Tier B (vlm_apply on a
subfile, joined on SeedTs against this reservoir), and only that pass can
see false positives. Read the numbers here as agreement on the tag's
positives, not as a classifier score; on 37Cl the tag itself is mostly
background, so disagreement is not automatically the model's error.

The tag carries no per-event label beyond its own bit, so the beam-flat
events are the only negative reference available here, and the AUC quoted
is tagged-vs-beam-flat: how well p(an) ranks the two populations, threshold
free. The seeds in config.VLM_SEEDS are all run, as in vlm_apply, so the
spread of the agreement fraction across them is on record.

Run inside the dataset dev shell, from python/:

    python vlm_reservoir.py                       # config.VLM_MODEL
    python vlm_reservoir.py --model google/gemma-4-E2B-it --load-in none
"""

import argparse
import re

import numpy as np

import config
import vlm

_CACHE_BASENAME = "StripSumScatter_cache.root"


def load_reservoir():
    """The `traces` tree as arrays.

    Returns (X, meta): X is (n, 18) float32 per-strip totals, strips 0 and
    17 from their scalars and the rest left+right, which is TraceEvt::Total
    on the C++ side. meta is a dict of seed_ts (uint64), reac_mask (uint32),
    reac (int, the LOWEST tagged strip or -1), beam_flat (bool), both_mult
    (int) and reac_min (the bit base, parsed from the cache fingerprint).
    """
    import ROOT

    path = config.ROOT_FILES_DIR / _CACHE_BASENAME
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise FileNotFoundError(f"no compute-regions cache at {path}")
    fp = f.Get("fingerprint")
    m = re.search(r"reac\[(\d+),(\d+)\]", fp.GetTitle() if fp else "")
    if not m:
        raise RuntimeError("cache fingerprint does not name the reaction "
                           "strip range; cannot decode reac_mask")
    reac_min = int(m.group(1))
    t = f.Get("traces")
    if not t:
        raise RuntimeError(f"no traces tree in {path}")
    n = t.GetEntries()
    X = np.empty((n, config.N_STRIPS), dtype=np.float32)
    seed_ts = np.empty(n, dtype=np.uint64)
    reac_mask = np.empty(n, dtype=np.uint32)
    beam_flat = np.empty(n, dtype=bool)
    both_mult = np.empty(n, dtype=np.int32)
    for i in range(n):
        t.GetEntry(i)
        X[i, 0] = t.strip0dE
        X[i, 17] = t.strip17dE
        for k in range(16):
            X[i, k + 1] = t.leftdE[k] + t.rightdE[k]
        seed_ts[i] = np.uint64(t.seed_ts)
        reac_mask[i] = np.uint32(t.reac_mask)
        beam_flat[i] = bool(t.beam_flat)
        both_mult[i] = int(t.both_mult)
    f.Close()
    # Lowest set bit -> strip. The C++ sets one bit per tagged strip; an
    # event tagged at two strips is reported, not silently split.
    reac = np.full(n, -1, dtype=np.int32)
    tagged = reac_mask > 0
    low = np.zeros(n, dtype=np.int32)
    mm = reac_mask.astype(np.int64)
    low[tagged] = np.array([int(v & -v).bit_length() - 1 for v in mm[tagged]])
    reac[tagged] = low[tagged] + reac_min
    multi = int((np.array([bin(int(v)).count("1") for v in mm]) > 1).sum())
    print(f"  {path.name}: {n} reservoir events, {int(tagged.sum())} tagged "
          f"(bit base reac {reac_min}), {int(beam_flat.sum())} beam-flat"
          + (f", {multi} tagged at more than one strip" if multi else ""))
    return X, {
        "seed_ts": seed_ts,
        "reac_mask": reac_mask,
        "reac": reac,
        "beam_flat": beam_flat,
        "both_mult": both_mult,
        "reac_min": reac_min,
    }


def rank_auc(pos, neg):
    """P(p_pos > p_neg) by ranks (Mann-Whitney), ties counted half."""
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(allv.size, dtype=np.float64)
    # Average ranks over ties.
    sorted_v = allv[order]
    i = 0
    while i < sorted_v.size:
        j = i
        while j + 1 < sorted_v.size and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    r_pos = ranks[:pos.size].sum()
    return float((r_pos - pos.size * (pos.size + 1) / 2.0) /
                 (pos.size * neg.size))


def report(probs, labels, meta, thr):
    """Per-group table: n, mean/median p(an), tagged fractions, argmax mix."""
    pan = vlm.an_probability(probs)
    classes = list(config.VLM_CLASSES)
    groups = [("beam-flat", meta["beam_flat"])]
    for r in sorted(set(int(v) for v in meta["reac"] if v >= 0)):
        groups.append((f"tag reac {r:2d}", meta["reac"] == r))
    groups.append(("tag (all)", meta["reac"] >= 0))
    head = (f"  {'group':<12} {'n':>5} {'<p(an)>':>8} {'median':>7} "
            f"{'p>=0.3':>7} {f'p>={thr}':>7} {'p>=0.7':>7}  argmax "
            + " ".join(f"{c:>5}" for c in classes))
    print(head)
    out = {}
    for name, sel in groups:
        n = int(sel.sum())
        if n == 0:
            continue
        p = pan[sel]
        lab = labels[sel]
        mix = " ".join(f"{100.0 * (lab == i).mean():4.0f}%"
                       for i in range(len(classes)))
        print(f"  {name:<12} {n:>5} {p.mean():8.3f} {np.median(p):7.3f} "
              f"{(p >= 0.3).mean():7.3f} {(p >= thr).mean():7.3f} "
              f"{(p >= 0.7).mean():7.3f}         {mix}")
        out[name] = (n, float(p.mean()), float((p >= thr).mean()))
    auc = rank_auc(pan[meta["reac"] >= 0], pan[meta["beam_flat"]])
    print(f"  AUC p(an) tagged vs beam-flat: {auc:.3f}  "
          "(0.5 = no ranking, 1 = tagged always above beam-flat)")
    return out, auc


def contact_sheet(X, beam_ref, pan, meta, tag):
    """Rows of beam-flat, then the tagged events with the highest and the
    lowest p(an): the two ends of the model's ranking, next to what it was
    told beam looks like. p(an) values are printed in the same order."""
    from PIL import Image

    per_row = config.VLM_SHEET_PER_ROW
    rng = np.random.default_rng(config.SEED)
    tagged = np.flatnonzero(meta["reac"] >= 0)
    flat = np.flatnonzero(meta["beam_flat"])
    order = tagged[np.argsort(pan[tagged])]
    rows = [
        ("beam-flat (random)", rng.choice(flat, min(per_row, flat.size),
                                          replace=False)),
        ("tagged, highest p(an)", order[::-1][:per_row]),
        ("tagged, lowest p(an)", order[:per_row]),
    ]
    imgs = []
    for name, idx in rows:
        img = vlm.render_traces(X[idx], beam_ref)
        pad = per_row - img.shape[0]
        if pad > 0:
            blank = np.empty((pad, ) + img.shape[1:], dtype=np.uint8)
            blank[:] = np.asarray(config.VLM_BG_RGB, dtype=np.uint8)
            img = np.concatenate([img, blank])
        imgs.append(np.concatenate(list(img), axis=1))
        print(f"    {name:<24} p(an): "
              + " ".join(f"{pan[i]:.2f}" for i in idx)
              + "   reac: " + " ".join(f"{meta['reac'][i]:d}" for i in idx))
    out = config.PLOTS_DIR / config.VLM_PLOT_SUBDIR
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"vlm_reservoir_sheet_{tag}.png"
    sheet = np.concatenate(imgs, axis=0)
    Image.fromarray(sheet).resize((sheet.shape[1] * 3, sheet.shape[0] * 3),
                                  Image.NEAREST).save(path)
    print(f"  contact sheet -> {path}")


def _model_tag(model_id):
    return re.sub(r"[^A-Za-z0-9]+", "_", model_id.split("/")[-1]).strip("_")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", default=None,
                    help="Hub id; default config.VLM_MODEL")
    ap.add_argument("--load-in", default=None, choices=("none", "8bit",
                                                         "4bit"),
                    help="override config.VLM_LOAD_IN ('none' = bf16)")
    ap.add_argument("--budget", type=int, default=None,
                    help="override config.VLM_TOKEN_BUDGET")
    ap.add_argument("--seeds", type=int, nargs="*", default=None,
                    help="override config.VLM_SEEDS")
    ap.add_argument("--adapter", default=None,
                    help="LoRA adapter directory for the configured base model")
    args = ap.parse_args()
    if args.model:
        config.VLM_MODEL = args.model
    if args.load_in:
        config.VLM_LOAD_IN = None if args.load_in == "none" else args.load_in
    if args.budget:
        config.VLM_TOKEN_BUDGET = args.budget
    if args.adapter:
        config.VLM_ADAPTER_PATH = args.adapter
    seeds = tuple(args.seeds) if args.seeds else config.VLM_SEEDS
    # seed_perturbation treats the first configured seed as the unperturbed
    # one; keep that meaning when the list is overridden.
    config.VLM_SEEDS = seeds
    # The shared-prefix cache is a throughput trick that moves p(an) by up
    # to ~0.06 on gemma4 and disables itself when it flips a tag. Fifteen
    # hundred events do not need the throughput; run them plain so the
    # numbers carry no cache systematic.
    config.VLM_PREFIX_CACHE = False

    np.random.seed(config.SEED)
    print(f"music-ml: dataset {config.DATASET}, VLM on the compute-regions "
          f"reservoir ({config.VLM_MODEL}, "
          f"{getattr(config, 'VLM_LOAD_IN', None) or config.VLM_DTYPE}, "
          f"budget {config.VLM_TOKEN_BUDGET})")
    X, meta = load_reservoir()
    long_w = config.block_widths()[0]
    X = X[:, :long_w]
    # The beam reference the prompt draws: the mean of the beam-flat events,
    # which is the same population compute-regions took its baseline from.
    flat = meta["beam_flat"]
    if not flat.any():
        raise RuntimeError("no beam-flat events in the reservoir; nothing "
                           "to draw as the beam reference")
    beam_ref = X[flat].astype(np.float64).mean(axis=0)
    # A strip the C++ ignores (IGNORE_STRIP_0/17) is zero on every event and
    # absent from the region_traces figure the prompt shows as its worked
    # example. Drawn, it would be a vertical drop off the canvas on every
    # trace, reference included -- a feature the reference figure does not
    # have. Trim it so the two views match.
    live = beam_ref > 0.0
    if not live.all():
        dead = [int(s) for s in np.flatnonzero(~live)]
        print(f"  dropping strip(s) {dead}: zero on every beam-flat event "
              "(ignored on the C++ side, not on the reference figure)")
        X = X[:, live]
        beam_ref = beam_ref[live]
    print("  beam reference (beam-flat mean) "
          + " ".join(f"{v:.2f}" for v in beam_ref))

    thr = config.VLM_AN_THRESHOLD
    tag = _model_tag(config.VLM_MODEL)
    summary = []
    for seed in seeds:
        order, dy, dhalf = vlm.seed_perturbation(seed)
        letters = ", ".join(f"{config.VLM_CLASS_TOKENS[i]}="
                            f"{config.VLM_CLASSES[c]}"
                            for i, c in enumerate(order))
        print(f"\nseed {seed}: {letters}; dy={dy:+.2f}px dhalf={dhalf:+.2f}px")
        clf = vlm.EventClassifier(letter_order=order)
        probs, labels = clf.classify(X, beam_ref, dy=dy, dhalf=dhalf)
        del clf
        groups, auc = report(probs, labels, meta, thr)
        summary.append((seed, groups, auc))
        pan = vlm.an_probability(probs)
        if seed == seeds[0]:
            contact_sheet(X, beam_ref, pan, meta, tag)
        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        out = config.CACHE_DIR / f"vlm_reservoir_{tag}_seed{seed}.npz"
        np.savez(out,
                 seed_ts=meta["seed_ts"],
                 reac_mask=meta["reac_mask"],
                 reac=meta["reac"],
                 beam_flat=meta["beam_flat"],
                 both_mult=meta["both_mult"],
                 probs=probs.astype(np.float32),
                 labels=labels.astype(np.int8),
                 classes=np.array(config.VLM_CLASSES),
                 model=np.array(config.VLM_MODEL),
                 budget=np.array(config.VLM_TOKEN_BUDGET))
        print(f"  per-event table -> {out}")

    print("\n=== summary across seeds ===")
    print(f"  {'seed':>5} {'tag agree':>10} {'flat false+':>12} {'AUC':>6}")
    for seed, groups, auc in summary:
        agree = groups.get("tag (all)", (0, 0.0, float("nan")))[2]
        fp = groups.get("beam-flat", (0, 0.0, float("nan")))[2]
        print(f"  {seed:>5} {agree:10.3f} {fp:12.3f} {auc:6.3f}")
    agrees = np.array([g.get("tag (all)", (0, 0, np.nan))[2]
                       for _, g, _ in summary])
    print(f"  tag agreement at p>={thr}: {agrees.mean():.3f} "
          f"+/- {agrees.std(ddof=0):.3f} across {agrees.size} seed(s)")
    print("done")


if __name__ == "__main__":
    main()
