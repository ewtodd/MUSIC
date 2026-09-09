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
## Documentation
<!---->
**→ [docs.ethanwtodd.com/music](https://docs.ethanwtodd.com/music/)**
<!---->
Guides and the full API reference live there. Start with
[Repository layout](https://docs.ethanwtodd.com/music/layout.html) for the monorepo and the config split,
then [Building](https://docs.ethanwtodd.com/music/building.html) and [Running](https://docs.ethanwtodd.com/music/running.html).
<!---->
## Quick start
<!---->
Requires [Nix](https://nixos.org/) with
[flakes enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently).
<!---->
```sh
nix build               # builds music-tooling-37Cl (the default package)
nix build .#87Rb        # the other dataset
nix develop .#37Cl      # editing only; gives the toolchain and include paths
```
<!---->
Binaries land in `./result/bin`, take no command-line arguments, and must be run
from inside `nix develop .#<dataset>` so they can find their runtime libraries.
The main entry point is `./result/bin/pipeline`; every knob lives in
`Constants.hpp` / `Constants.cpp` and the `control/` TOMLs.
<!---->
There is no `make` step you run yourself — `tooling/Makefile` is an internal
recipe the flake drives.
<!---->
## Layout
<!---->
| Path | Contents |
|------|----------|
| `tooling/include/`, `tooling/src/` | the shared C++ analysis library |
| `tooling/mains/` | one `main_*.cpp` per binary |
| `tooling/gpu/` | CUDA timestamp-sort kernel, `dlopen`'d at runtime |
| `analysis/<dataset>/` | per-dataset config, control TOMLs, and output |
| `doc/pages/` | narrative documentation, rendered by Doxygen |
| `doc/theme/` | Kanagawa palette for the generated site |
| `python/` | Python library for machine learning |
<!---->
The tooling is identical across datasets; only
`analysis/<dataset>/config/Constants.cpp` and the `control/` TOMLs differ. See
[Repository layout](https://docs.ethanwtodd.com/music/layout.html) for how that split works.
<!---->
## Building the docs locally
<!---->
`doxygen` and `graphviz` are in the dev shells already:
<!---->
```sh
nix develop .#37Cl
./scripts/build_docs.sh
```
<!---->
Output lands in `docs/html/index.html` (gitignored), warnings in
`docs/doxygen-warnings.log`. `nix build .#docs` produces the same site as a
store path, and fails if doxygen emits any warning — including an undocumented
public declaration.
<!---->
## Dependencies
<!---->
Everything is pinned with Nix. The flake provides ROOT (CUDA build), a C++
toolchain, `tomlplusplus`, and the shared
[Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities) library. A
CUDA-capable GPU is used for the timestamp sort; without the CUDA toolchain it
falls back to CPU.
<!---->
## AI assisted development disclosure
<!---->
The API documentation — the Doxygen comments on the headers and the narrative
guides — was written with the help of
[Claude Code](https://claude.ai/claude-code) and/or
[son-of-anton](https://github.com/ewtodd/son-of-anton).
<!---->
Those comments describe behaviour derived by reading the implementations rather
than from a specification, so where the documentation and the code disagree, the
code is what runs. The analysis logic itself — event building, calibration,
tagging, cross-section extraction — is human-written and human-reviewed.
<!---->
## Git hooks
<!---->
`.githooks/pre-commit` runs `clang-format` (LLVM style) on staged C/C++ files.
Enable it once with `git config core.hooksPath .githooks`.
<!---->
## License
<!---->
[MIT](LICENSE)
