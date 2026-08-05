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
            cudaCapabilities = [ "8.9" ];
            cudaForwardCompat = false;
          };
          overlays = [
            # torch-bin 2.12 requires cuda-bindings >= 13.0.3; the default
            # cudaPackages (12.9) is too old, so pin the 13.0 set instead.
            (final: prev: {
              cudaPackages = prev.cudaPackages_13_0;
            })
            # slicer (a dep of shap) fails its test suite on python 3.14 in
            # the updated nixpkgs; skip the checks.
            (final: prev: {
              pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
                (python-final: python-prev: {
                  # slicer (a dep of shap) is broken in the updated nixpkgs:
                  # its pyproject build lacks setuptools, and its test suite
                  # fails on python 3.14. Fix the build and skip the checks.
                  slicer = python-prev.slicer.overridePythonAttrs (old: {
                    doCheck = false;
                    nativeBuildInputs = (old.nativeBuildInputs or []) ++ [
                      python-final.setuptools
                    ];
                  });
                  # shap is missing typing-extensions in its nixpkgs deps; the
                  # runtime-deps check fails without it.
                  shap = python-prev.shap.overridePythonAttrs (old: {
                    dependencies = (old.dependencies or []) ++ [
                      python-final.typing-extensions
                    ];
                  });
                  # torch-bin pins setuptools<82 but nixpkgs ships 83.0.0;
                  # relax the constraint so the runtime-deps check passes.
                  torch-bin = python-prev.torch-bin.overridePythonAttrs (old: {
                    pythonRelaxDeps = (old.pythonRelaxDeps or []) ++ [
                      "setuptools"
                    ];
                  });
                })
              ];
            })
          ];
        };
        isCUDA = true;
        analysis-utils =
          if !isCUDA then utils.packages.${system}.default else utils.packages.${system}.cuda;
        analysis-utils-py = utils.packages.${system}.pythonPackage;
        root = if !isCUDA then pkgs.root else utils.packages.${system}.rootCuda;
        clangdConfigFile = (pkgs.formats.yaml { }).generate "dot-clangd" {
          CompileFlags.Add = [
            "--cuda-gpu-arch=sm_89"
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
              pkgs.cudaPackages.cccl
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
      let
        mkPackage =
          dataset:
          pkgs.stdenv.mkDerivation {
            name = "music-tooling-${dataset}";
            src = ./.;
            nativeBuildInputs = with pkgs; [
              pkg-config
              gnumake
            ];
            buildInputs = [
              analysis-utils
              root
              pkgs.bash
              pkgs.tomlplusplus
            ]
            ++ pkgs.lib.optionals isCUDA [
              pkgs.cudaPackages.cuda_nvcc
              pkgs.cudaPackages.cuda_cudart
              pkgs.cudaPackages.cccl
            ];

            buildPhase = ''
              export MUSIC_DATASET="${dataset}"
              export MUSIC_DATASET_DIR="$sourceRoot/analysis/${dataset}"
              make -j DATASET_DIR_OUT="$out/analysis/${dataset}" \
                   GPU_LIB_OUT="$out/lib/libgpuaccel.so"
            '';

            installPhase = ''
              mkdir -p $out/bin $out/lib $out/analysis/${dataset}/config
              cp analysis/${dataset}/bin/* $out/bin/
              cp tooling/gpu/libgpuaccel.so $out/lib/
              cp -r analysis/${dataset}/config/* $out/analysis/${dataset}/config/
            '';
          };
      in
      {
        devShells = {
          "87Rb" = mkDatasetShell "87Rb";
          "37Cl" = mkDatasetShell "37Cl";
          default = mkDatasetShell "37Cl";
        };

        packages = {
          "87Rb" = mkPackage "87Rb";
          "37Cl" = mkPackage "37Cl";
          default = mkPackage "37Cl";
        };
      }
    );
}
