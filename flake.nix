{
  description = "MUSIC analysis monorepo (one tooling, per-dataset config)";
  inputs = {
    # Follow Analysis-Utilities' nixpkgs so both flakes share one interpreter
    # and toolchain: its pythonPackage is then importable from this flake's
    # python3.withPackages (a foreign-python module would be silently
    # dropped), and the C++ here compiles with the same stdenv that built
    # ROOT.
    nixpkgs.follows = "utils/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    utils = {
      url = "/home/e-work/Analysis-Utilities";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaCapabilities = [ "12.0" ];
            cudaForwardCompat = false;
          };
        };
        isCUDA = true;
        analysis-utils =
          if !isCUDA then utils.packages.${system}.default else utils.packages.${system}.cuda;
        analysis-utils-py = utils.packages.${system}.pythonPackage;
        root = if !isCUDA then pkgs.root else utils.packages.${system}.rootCuda;
        clangdConfigFile = (pkgs.formats.yaml { }).generate "dot-clangd" {
          CompileFlags.Add = [
            "--cuda-gpu-arch=sm_120"
            "--no-cuda-version-check"
          ];
          Diagnostics.Suppress = [
            "no_member"
            "nested_name_spec_non_tag"
            "typename_nested_not_found"
            "template_instantiate_undefined"
          ];
        };

        # One dev shell per dataset. The only difference is the exported
        # MUSIC_DATASET / MUSIC_DATASET_DIR (which the Makefile and the tooling
        # banner / path resolution read). All datasets share the one tooling.
        mkDatasetShell =
          dataset:
          pkgs.mkShell {
            nativeBuildInputs = with pkgs; [
              pkg-config
              gnumake
              clang-tools
            ];
            buildInputs = [
              analysis-utils
              root
              pkgs.bash
              pkgs.tomlplusplus
              (pkgs.python3.withPackages (
                python-pkgs: with python-pkgs; [
                  numpy
                  pandas
                  scikit-learn
                  scipy
                  shap
                  packaging
                  torch-bin
                  xgboost
                  analysis-utils-py
                ]
              ))
            ]
            ++ pkgs.lib.optionals isCUDA [
              pkgs.cudaPackages.cuda_nvcc
              pkgs.cudaPackages.cuda_cudart
              pkgs.cudaPackages.cuda_cccl
            ];
            shellHook = ''
               echo "Analysis-Utilities version: ${analysis-utils.version}${pkgs.lib.optionalString isCUDA " (CUDA)"}"
               flake_root="$PWD"
               git_root="$(git -C "$flake_root" rev-parse --show-toplevel)"

              # --- dataset selection ---
              export MUSIC_DATASET="${dataset}"
              export MUSIC_DATASET_DIR="$git_root/analysis/${dataset}"
              echo "MUSIC dataset: ${dataset}  ($MUSIC_DATASET_DIR)"

              # Where GENERATED outputs (root_files, plots) land. Defaults to the
              # in-repo dataset dir (on /home); override this one var to redirect
              # processed output to a scratch drive without rebuilding.
              export MUSIC_RESULTS_DIR="''${MUSIC_RESULTS_DIR:-$git_root/analysis/${dataset}}"
              echo "MUSIC results: $MUSIC_RESULTS_DIR"

              ${pkgs.lib.optionalString isCUDA ''
                export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
                export LD_LIBRARY_PATH="/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
                install -m 644 ${clangdConfigFile} "$git_root/tooling/gpu/.clangd"
              ''}

              # Editor/clangd include resolution: tooling headers + this dataset's config.
              export CPLUS_INCLUDE_PATH="$git_root/tooling/include:$MUSIC_DATASET_DIR/config''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
              export ROOT_INCLUDE_PATH="$git_root/tooling/include:$MUSIC_DATASET_DIR/config''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
              export LD_LIBRARY_PATH="$git_root/tooling/gpu''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

              alias wipe-analysis='rm -r analysis/${dataset}/plots analysis/${dataset}/root_files'
              alias wipe-plots='rm -r analysis/${dataset}/plots'
              alias wipe-root='rm -r analysis/${dataset}/root_files'
              alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm *_cxx.so *_cxx.d *_cxx_ACLiC_dict_rdict.pcm AutoDict_*'
            '';
          };
      in
      {
        devShells = {
          "87Rb" = mkDatasetShell "87Rb";
          "37Cl" = mkDatasetShell "37Cl";
          default = mkDatasetShell "87Rb";
        };
      }
    );
}
