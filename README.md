# MUSIC (Multi-Sampling Ionization Chamber at Argonne National Lab) Analysis
<!---->
## Projects
<!---->
- **ProductionMode_37Cl** — Event building and timing analysis for the MUSIC detector array, including binary data conversion, multi-board timing synchronization, and event reconstruction.
- **SiCalibration_37Cl** — Silicon diode energy loss measurement used for validating the stopping power tables for 37Cl in helium gas.
<!---->
## Dependencies
<!---->
All projects depend on [Analysis-Utilities](https://github.com/ewtodd/Analysis-Utilities), a shared library providing common ROOT-based analysis tools.
<!---->
## Reproducibility
<!---->
This project uses [Nix](https://nixos.org/) to manage dependencies.
Install it by following the instructions [here](https://nixos.org/download/).
Ensure [flakes are enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently) in your Nix configuration.
All directories contain their own Nix flake, which can be used to re-run code with the exact same environment as I did.
If projects reference a local version of Analysis-Utilities in the flake, they are not yet finished.
Contact me for access to data.
