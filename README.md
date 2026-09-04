# MUSIC — Multi-Sampling Ionization Chamber Analysis

```
 ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ 
  ███╗   ███╗ ██╗   ██╗ ███████╗ ██╗  ██████╗
  ████╗ ████║ ██║   ██║ ██╔════╝ ██║ ██╔════╝
  ██╔████╔██║ ██║   ██║ ███████╗ ██║ ██║
  ██║╚██╔╝██║ ██║   ██║ ╚════██║ ██║ ██║
  ██║ ╚═╝ ██║ ╚██████╔╝ ███████║ ██║ ╚██████╗
  ╚═╝     ╚═╝  ╚═════╝  ╚══════╝ ╚═╝  ╚═════╝
 ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ 
```

<!---->
Analysis code for data from the MUSIC active-target ionization chamber at
Argonne National Lab.
This is a **monorepo**: a single shared tooling tree is compiled separately
against each experiment's configuration.
Each experiment ("dataset", e.g.
`37Cl`, `87Rb`) lives under `analysis/<dataset>/` and is built from the one copy
of the code in `tooling/`.
<!---->
## Layout
<!---->
```
python/                 python library for machine learning
tooling/                shared source — built once per dataset
    include/  src/        C++ analysis library + headers
    mains/                one main_*.cpp per binary
    gpu/                  CUDA timestamp-sort kernel (libgpuaccel.so), dlopen'd at runtime
analysis/<dataset>/
    config/Constants.cpp  the per-dataset config overrides (run numbers, thresholds, gates, ...)
    control/              per-reaction TOML control files (gas, beam, target, detector)
    bin/                  built binaries (git-ignored; produced only by the internal Makefile)
    build/                object files + libmusic.a (git-ignored; internal Makefile)
    root_files/           event ROOT files produced by the pipeline
    sim_root_files/       simulated beam/trace ROOT files
    plots/                output figures
flake.nix               Nix dev shells + package derivation, one per dataset
tooling/Makefile        internal build recipe, driven by the flake (never run directly)
```
<!---->
The config is split in two: the struct definition and tooling-wide defaults live
in `tooling/include/Constants.hpp` / `tooling/src/Constants.cpp`, and each dataset
overrides the fields it cares about in `analysis/<dataset>/config/Constants.cpp`
(compiled into that dataset's binaries). The `control/` TOMLs are the only other
per-dataset input. The tooling itself is identical.
<!---->
## Dependencies
<!---->
Everything is pinned with [Nix](https://nixos.org/) (flakes must be
[enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently)). The flake
provides ROOT (CUDA build), a C++ toolchain, `tomlplusplus`, and the shared
[Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities) library. A
CUDA-capable GPU is used for the timestamp sort; without the CUDA toolchain
(the build skips it, or `isLaptop = true` in `flake.nix`) it falls back to CPU.
<!---->
## Building
<!---->
Binaries are made with `nix build`; `nix develop` is only for editing code.
There is no `make` step you run yourself — `tooling/Makefile` is an internal
recipe that the flake drives (as `make -C tooling`). It is no longer at the repo
root, so `make` cannot be run from there; invoking it by hand would drop binaries
into `analysis/<dataset>/bin/` instead of keeping them at the repo root.
<!---->
To edit / develop — this gives the compiler toolchain and include paths so the
LSP (clangd) can index the code (it does **not** build anything):
<!---->
```sh
nix develop .#37Cl      # default shell is 37Cl; use .#87Rb for the other dataset
```
<!---->
To build the binaries:
<!---->
```sh
nix build               # builds music-tooling-37Cl (the default package)
nix build .#87Rb        # the other dataset
./result/bin/...        # binaries land under ./result/bin (repo root)
```
<!---->
`./result` is a symlink to the built package's store path. Each binary resolves
its dataset directory at build time (baked in) and cross-checks the dataset name
if `MUSIC_DATASET_DIR` is overridden. The binary still needs the dev-shell
environment to find its runtime libraries (e.g. `libanalysis-utils.so`), so run
it from inside `nix develop .#<dataset>` rather than a bare shell. The binaries
take no command-line arguments — all knobs live in `Constants.hpp`/`Constants.cpp`
and the `control/` TOMLs.
<!---->
## Running
<!---->
The main entry point is the **pipeline**, which for each CoMPASS binary subfile
does: binary → raw hits → multi-board timing alignment → GPU timestamp sort →
event building → beam-energy calibration → trace creation. It runs the dataset's
files (set by `RUNS` / `N_CHUNKS` in `analysis/<dataset>/config/Constants.cpp`)
in parallel and logs to `analysis/<dataset>/pipeline_fused.log`:
<!---->
```sh
nix build .#87Rb
nix develop .#87Rb
./result/bin/pipeline
```
<!---->
The remaining binaries operate on the pipeline's output ROOT files:
`calibrate-beam` (per-channel beam-peak calibration, L/R gain matching, strip
alignment, eres aggregation), `strip-sum-scatter` (reaction-strip scatter/dE
plots, trace region overlays, sim comparisons), and the SOL-processing tools
(`split-sol`, `preprocess-sol`) — all in `./result/bin/` after `nix build
.#<dataset>`.
<!---->
Generated data (event ROOT files, plots) is written to `MUSIC_RESULTS_DIR` when
set, else to the baked-in dataset dir; the pipeline log goes to
`<dataset dir>/pipeline_fused.log`. You run from the dev shell, so both
`MUSIC_DATASET_DIR` and `MUSIC_RESULTS_DIR` point at the in-repo
`analysis/<dataset>` (writable) — override `MUSIC_RESULTS_DIR` to send outputs
to a scratch drive.

`strip-sum-scatter` caches its scatters and trace reservoir to
`analysis/<dataset>/root_files/StripSumScatter_cache.root`; the cache carries a
fingerprint of every config knob that changes its contents (filters, gates,
strip spans, ...), so editing any of them rebuilds the cache automatically on
the next run.

The TOML control files in `analysis/<dataset>/control/` (gas fill, beam
species/energy, target, detector response) are used to generate simulated data
for the relevant reaction using [Remix-MUSIC-Sim](https://github.com/ewtodd/Remix-MUSIC-Sim). <!---->
## Git hooks
<!---->
`.githooks/pre-commit` runs `clang-format` (LLVM style) on staged C/C++ files.
Enable it once with `git config core.hooksPath .githooks`.
