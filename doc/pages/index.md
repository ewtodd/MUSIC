MUSIC
=====

Analysis code for data from the MUSIC active-target ionization chamber at
Argonne National Laboratory.

This is a **monorepo**: a single shared tooling tree under `tooling/` is
compiled separately against each experiment's configuration. Each experiment
("dataset", e.g. `37Cl`, `87Rb`) lives under `analysis/<dataset>/` and is built
from that one copy of the code.

The analysis library is built on [ROOT](https://root.cern/) and on
[Analysis-Utilities](https://docs.ethanwtodd.com/analysis-utilities/), which
supplies the binary readers, waveform processing, and photopeak fitting this
code sits on top of.

## Guides

- @subpage layout — the monorepo, and where the per-dataset configuration lives
- @subpage building — Nix dev shells, building the binaries
- @subpage running — the pipeline and the binaries that consume its output

## API reference

Every class and function in `tooling/` is indexed under
<a href="annotated.html">Classes</a> and <a href="files.html">Files</a>. The
search box covers the whole reference.

@note This site documents `tooling/`, the shared analysis library, which is
identical across datasets. Per-dataset configuration and experiment inputs are
not part of the API and are not published here.
