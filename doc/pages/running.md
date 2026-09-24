Running {#running}
=======

[TOC]

## The pipeline

The main entry point is the **pipeline**. For each CoMPASS binary subfile it
runs: binary → raw hits → multi-board timing alignment → GPU timestamp sort →
event building → beam-energy calibration → trace creation.

It processes the dataset's files — selected by `RUNS` and `N_CHUNKS` in
`analysis/<dataset>/config/Constants.cpp` — in parallel, logging to
`analysis/<dataset>/pipeline_fused.log`.

```sh
nix build .#87Rb
nix develop .#87Rb
./result/bin/pipeline
```

@note The binaries take no command-line arguments. Every knob lives in
`Constants.hpp` / `Constants.cpp` and in the `control/` TOMLs. Each binary
resolves its dataset directory at build time and cross-checks the dataset name
if `MUSIC_DATASET_DIR` is overridden.

## The other binaries

These operate on the pipeline's output ROOT files, and are all in
`./result/bin/` after `nix build .#<dataset>`:

- **`calibrate-beam`** — per-channel beam-peak calibration, left/right gain
  matching, strip alignment, energy-resolution aggregation.
- **`strip-sum-scatter`** — reaction-strip scatter and dE plots, trace region
  overlays, comparisons against simulation. Also writes the event selection as
  configured, cut by cut in the order applied, to
  `plots/strip_sum_scatter/selection.{dot,mmd}` (and `.png`/`.pdf` when
  Graphviz's `dot` is on the PATH, as it is in the dev shells); `compute-regions`
  writes the same in the all-tagged mode.
- **`preprocess-sol`** — the SOL-processing tool: splits `.sol` files into
  time chunks and decodes them to raw-hit ROOT files.
- **`cross-section`** — the absolute cross section per reaction strip from the
  tagged counts and the region cut, with the cross-section figures.
- **`talys-xs`** — runs TALYS for each configured model over the strips'
  energies and writes the residual channels to `root_files/talys/`.

## Where output goes

Generated data — event ROOT files and plots — is written to `MUSIC_RESULTS_DIR`
when that is set, and otherwise to the baked-in dataset directory. The pipeline
log always goes to `<dataset dir>/pipeline_fused.log`.

Running from the dev shell points both `MUSIC_DATASET_DIR` and
`MUSIC_RESULTS_DIR` at the in-repo `analysis/<dataset>`, which is writable.
Override `MUSIC_RESULTS_DIR` to send output to a scratch drive instead.

## Caching

`strip-sum-scatter` caches its scatters and trace reservoir to
`analysis/<dataset>/root_files/StripSumScatter_cache.root`.

The cache carries a fingerprint of every configuration knob that changes its
contents — filters, gates, strip spans, and so on — so editing any of them
rebuilds the cache automatically on the next run. There is no cache-clearing
step to remember.
