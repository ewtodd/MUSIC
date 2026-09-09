Building {#building}
========

[TOC]

## Dependencies

Everything is pinned with [Nix](https://nixos.org/); flakes must be
[enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently). The flake
provides ROOT (CUDA build), a C++ toolchain, `tomlplusplus`, and the shared
[Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities) library.

A CUDA-capable GPU is used for the timestamp sort. Without the CUDA toolchain —
either because the build skipped it, or because `isLaptop = true` in
`flake.nix` — it falls back to CPU.

## Building the binaries

```sh
nix build               # builds music-tooling-37Cl (the default package)
nix build .#87Rb        # the other dataset
./result/bin/...        # binaries land under ./result/bin
```

`./result` is a symlink to the built package's store path.

@warning There is no `make` step you run yourself. `tooling/Makefile` is an
internal recipe that the flake drives as `make -C tooling`. It deliberately does
not sit at the repository root, so `make` cannot be run from there; invoking it
by hand would drop binaries into `analysis/<dataset>/bin/` instead of keeping
them at the repository root.

## Editing

```sh
nix develop .#37Cl      # default shell is 37Cl; use .#87Rb for the other
```

This gives the compiler toolchain and include paths so clangd can index the
code. It does **not** build anything.

@note A built binary still needs the dev-shell environment to find its runtime
libraries — `libanalysis-utils.so` among them — so run it from inside
`nix develop .#<dataset>` rather than a bare shell.

## Building this documentation

`doxygen` and `graphviz` are already in the dev shells:

```sh
nix develop .#37Cl
./scripts/build_docs.sh
```

Output lands in `docs/html/index.html` (gitignored), with any parse warnings in
`docs/doxygen-warnings.log`. `nix build .#docs` produces the same site as a
store path, which is what the published site is served from.

The reference is dataset-independent: it documents `tooling/`, which every
dataset shares, so there is one docs output rather than one per dataset.
