# MUSIC — Multi-Sampling Ionization Chamber Analysis
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
    bin/                  built binaries (git-ignored)
    build/                object files + libmusic.a (git-ignored)
    root_files/           event ROOT files produced by the pipeline
    sim_root_files/       simulated beam/trace ROOT files
    plots/                output figures
flake.nix               Nix dev shells + package derivation, one per dataset
Makefile                dataset-aware build
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
CUDA-capable GPU is used for the timestamp sort; toggle `isCUDA` in `flake.nix`
to disable it.
<!---->
## Building
<!---->
Enter a dataset's dev shell (this exports `MUSIC_DATASET` / `MUSIC_DATASET_DIR`,
which the Makefile and the tooling read), then build:
<!---->
```sh
nix develop .#37Cl      # or .#87Rb  (default shell is 37Cl)
make                    # build every binary for this dataset
make pipeline           # build a single binary
make clean              # remove this dataset's build artifacts
```
<!---->
You can also select the dataset directly without the shell:
<!---->
```sh
make DATASET=37Cl pipeline
```
<!---->
The packaged build (what `./result` points at) is driven by the flake's
derivation, which compiles the tooling in a clean source tree:
<!---->
```sh
nix build               # builds music-tooling-37Cl (the default package)
nix build .#87Rb        # the other dataset
./result/bin/...        # binaries land under ./result/bin
```
<!---->
With `make`, binaries land in `analysis/<dataset>/bin/`. Either way, each binary
resolves its dataset directory from `MUSIC_DATASET_DIR` at runtime, so run them
from inside the dev shell (or with that variable exported). They take no
command-line arguments — all knobs live in `Constants.hpp`/`Constants.cpp` and
the `control/` TOMLs.
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
nix develop .#87Rb
make pipeline
./analysis/87Rb/bin/pipeline
```
<!---->
The remaining binaries operate on the pipeline's output ROOT files:
`calibrate-beam` (per-channel beam-peak calibration, L/R gain matching, strip
alignment, eres aggregation), `strip-sum-scatter` (reaction-strip scatter/dE
plots, trace region overlays, sim comparisons), and the SOL-processing tools
(`split-sol`, `preprocess-sol`).

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
