# MUSIC — notes for agents

- `analysis/<dataset>/sim_control/*.toml` are Remix-MUSIC-Sim controls
  (gitignored). They follow that repo's strict schema (its README, "Strict TOML
  control schema"); validate with `musicsim --check <file>`. When the
  simulator's schema changes, migrate them by hand.
- `RemixSim` reads `sim_control/` (`RemixSim::ControlDir`). It uses only each
  control's `output = ".../traces_<iso>_<tag>.root"` line, loads
  `analysis/<dataset>/sim_root_files/traces_<iso>_<tag>.root`, and stops if two
  controls cover the same reaction and strip.
- Excitation, angular, and branch-weight choices in those controls are physics
  decisions for ewt.

## Commits

- A sentence-case summary line, a blank line, then a body that explains what
  changed and why. Never cram several changes into one long subject line.
- Agent commits: `--author="son-of-anton-bot <307402699+son-of-anton-bot@users.noreply.github.com>"`.
