Repository layout {#layout}
=================

[TOC]

## The monorepo

```
python/                 python library for machine learning
tooling/                shared source — built once per dataset
    include/  src/        C++ analysis library + headers
    mains/                one main_*.cpp per binary
    gpu/                  CUDA timestamp-sort kernel (libgpuaccel.so), dlopen'd at runtime
analysis/<dataset>/
    config/Constants.cpp  per-dataset config overrides (run numbers, thresholds, gates, ...)
    control/              per-reaction TOML control files (gas, beam, target, detector)
    bin/                  built binaries (git-ignored; produced only by the internal Makefile)
    build/                object files + libmusic.a (git-ignored; internal Makefile)
    root_files/           event ROOT files produced by the pipeline
    sim_root_files/       simulated beam/trace ROOT files
    plots/                output figures
flake.nix               Nix dev shells + package derivation, one per dataset
tooling/Makefile        internal build recipe, driven by the flake (never run directly)
```

## How configuration is split

The configuration is split in two. The struct definition and the tooling-wide
defaults live in `tooling/include/Constants.hpp` and `tooling/src/Constants.cpp`.
Each dataset then overrides only the fields it cares about, in
`analysis/<dataset>/config/Constants.cpp`, which is compiled into that dataset's
binaries.

The `control/` TOML files are the only other per-dataset input. **The tooling
itself is identical across datasets** — which is what makes a single API
reference meaningful for all of them.

The TOMLs (gas fill, beam species and energy, target, detector response) are
also what generate simulated data for the relevant reaction, via
[Remix-MUSIC-Sim](https://github.com/ewtodd/Remix-MUSIC-Sim).

## Git hooks

`.githooks/pre-commit` runs `clang-format` (LLVM style) on staged C/C++ files.
Enable it once with:

```sh
git config core.hooksPath .githooks
```
