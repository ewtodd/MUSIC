# MUSIC — Multi-Sampling Ionization Chamber Analysis

Analysis code for data from the MUSIC active-target ionization chamber at
Argonne National Lab. The detector measures the energy loss (dE/dx) of a heavy
ion beam strip-by-strip as it traverses a gas-filled chamber, which is used to
identify reaction products and reconstruct reactions such as `(α,α)` and `(α,n)`.

This is a **monorepo**: a single shared tooling tree is compiled separately
against each experiment's configuration. Each experiment ("dataset", e.g.
`37Cl`, `87Rb`) lives under `analysis/<dataset>/` and is built from the one copy
of the code in `tooling/`.

## Layout

```
tooling/                shared source — built once per dataset
  include/  src/        C++ analysis library + headers
  mains/                one main_*.cpp per binary
  gpu/                  CUDA timestamp-sort kernel (libgpuaccel.so), dlopen'd at runtime
analysis/<dataset>/
  config/Constants.hpp  the per-dataset config (run numbers, channel map, gates, ...)
  control/              per-reaction TOML control files (gas, beam, target, detector)
  bin/                  built binaries (git-ignored)
  build/                object files + libmusic.a (git-ignored)
  root_files/           event ROOT files produced by the pipeline
  sim_root_files/       simulated beam/trace ROOT files
  plots/                output figures
flake.nix               Nix dev shells, one per dataset
Makefile                dataset-aware build
```

The only thing that differs between datasets is `analysis/<dataset>/config/Constants.hpp`
and the `control/` TOMLs. The tooling is identical.

## Dependencies

Everything is pinned with [Nix](https://nixos.org/) (flakes must be
[enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently)). The flake
provides ROOT (CUDA build), a C++ toolchain, `tomlplusplus`, and the shared
[Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities) library. A
CUDA-capable GPU is used for the timestamp sort; toggle `isCUDA` in `flake.nix`
to disable it.

## Building

Enter a dataset's dev shell (this exports `MUSIC_DATASET` / `MUSIC_DATASET_DIR`,
which the Makefile and the tooling read), then build:

```sh
nix develop .#87Rb      # or .#37Cl  (default shell is 37Cl)
make                    # build every binary for this dataset
make pipeline           # build a single binary
make clean              # remove this dataset's build artifacts
```

You can also select the dataset directly without the shell:

```sh
make DATASET=37Cl pipeline
```

Binaries land in `analysis/<dataset>/bin/`. Each one resolves its dataset
directory from `MUSIC_DATASET_DIR` at runtime, so run them from inside the dev
shell (or with that variable exported). They take no command-line arguments —
all knobs live in `config/Constants.hpp` and the `control/` TOMLs.

## Running

The main entry point is the **pipeline**, which for each CoMPASS binary subfile
does: binary → raw hits → multi-board timing alignment → GPU timestamp sort →
event building → beam-energy calibration → trace creation. It runs the dataset's
files (set by `RUN_NUMBERS` / `N_FILES` in `Constants.hpp`) in parallel and logs
to `analysis/<dataset>/pipeline_fused.log`:

```sh
nix develop .#87Rb
make pipeline
./analysis/87Rb/bin/pipeline
```

The remaining binaries operate on the pipeline's output ROOT files: scatter/dE
plots, beam and silicon calibrations, stopping-power tables, simulation, and
timing/event diagnostics. Reaction-specific tools read the TOML control files in
`analysis/<dataset>/control/` (gas fill, beam species/energy, target, detector
response). Build them with `make` and run the corresponding binary from
`analysis/<dataset>/bin/`.

## Git hooks

`.githooks/pre-commit` runs `clang-format` (LLVM style) on staged C/C++ files.
Enable it once with `git config core.hooksPath .githooks`.
