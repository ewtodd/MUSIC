{
  description = "MUSIC analysis monorepo (unified analysis framework, per-dataset config)";
  inputs = {
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
        isLaptop = false;
        includePython = false;
        lib = nixpkgs.lib;
        pkgs = import nixpkgs {
          inherit system;
          config = lib.mkIf (!isLaptop) {
            allowUnfree = true;
            cudaCapabilities = [ "8.9" ];
            cudaForwardCompat = false;
          };
          overlays =
            if (!isLaptop && includePython == true) then
              [
                (final: prev: {
                  cudaPackages = prev.cudaPackages_13_0;
                })
                (final: prev: {
                  pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
                    (python-final: python-prev: {
                      slicer = python-prev.slicer.overridePythonAttrs (old: {
                        doCheck = false;
                        nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [
                          python-final.setuptools
                        ];
                      });
                      shap = python-prev.shap.overridePythonAttrs (old: {
                        dependencies = (old.dependencies or [ ]) ++ [
                          python-final.typing-extensions
                        ];
                      });
                      torch-bin = python-prev.torch-bin.overridePythonAttrs (old: {
                        pythonRelaxDeps = (old.pythonRelaxDeps or [ ]) ++ [
                          "setuptools"
                        ];
                      });
                    })
                  ];
                })
              ]
            else
              [ ];
        };
        analysis-utils =
          if isLaptop then utils.packages.${system}.default else utils.packages.${system}.cuda;
        analysis-utils-py = utils.packages.${system}.pythonPackage;
        root = if isLaptop then pkgs.root else utils.packages.${system}.rootCuda;
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
            ]
            ++ pkgs.lib.optionals includePython [
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
            ++ pkgs.lib.optionals (!isLaptop) [
              pkgs.cudaPackages.cuda_nvcc
              pkgs.cudaPackages.cuda_cudart
              pkgs.cudaPackages.cccl
            ];
            shellHook = ''
               echo "Analysis-Utilities version: ${analysis-utils.version}${
                 pkgs.lib.optionalString (!isLaptop) " (CUDA)"
               }"
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

              ${pkgs.lib.optionalString (!isLaptop) ''
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
          let
            extraBuild = if (!isLaptop) then ''GPU_LIB_OUT="$out/lib/libgpuaccel.so"'' else "";
            extraInstall = if (!isLaptop) then "cp tooling/gpu/libgpuaccel.so $out/lib/" else "";
          in
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
            ++ pkgs.lib.optionals (!isLaptop) [
              pkgs.cudaPackages.cuda_nvcc
              pkgs.cudaPackages.cuda_cudart
              pkgs.cudaPackages.cccl
            ];
            buildPhase = ''
              export MUSIC_DATASET="${dataset}"
              export MUSIC_DATASET_DIR="$sourceRoot/analysis/${dataset}"
              make -j DATASET_DIR_OUT="$out/analysis/${dataset}" ${extraBuild}
            '';

            installPhase = ''
              mkdir -p $out/bin $out/lib $out/analysis/${dataset}/config
              cp analysis/${dataset}/bin/* $out/bin/ 
              ${extraInstall}
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
