{
  description = "MUSIC analysis monorepo (unified analysis framework, per-dataset config)";
  inputs = {
    nixpkgs.follows = "utils/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    # Deliberately does NOT follow our nixpkgs: overriding it changes the
    # derivation hash and forces a local rebuild of the CUDA-overlaid ROOT that
    # cache.ethanwtodd.com cannot then satisfy. Our nixpkgs follows this input's
    # instead, so the whole tree stays on one pin.
    utils = {
      url = "github:ewtodd/Analysis-Utilities";
    };
    # TALYS, driven by talys-xs for the Hauser-Feshbach curve on the
    # cross-section plot. Its store path is compiled into the tooling.
    talys-nix = {
      url = "github:ewtodd/talys-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      utils,
      talys-nix,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        isLaptop = false;
        includePython = true;
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
                      # accelerate's MultiCPUTester spawns worker
                      # processes, which torch.multiprocessing cannot do in
                      # the nix sandbox; the other 238 tests pass. Disable
                      # the file rather than the whole checkPhase.
                      accelerate = python-prev.accelerate.overridePythonAttrs (old: {
                        disabledTestPaths = (old.disabledTestPaths or [ ]) ++ [
                          "tests/test_cpu.py"
                        ];
                      });
                      torch-bin = python-prev.torch-bin.overridePythonAttrs (old: {
                        pythonRelaxDeps = (old.pythonRelaxDeps or [ ]) ++ [
                          "setuptools"
                        ];
                      });
                      # The 12B rung of the VLM ladder is the gemma4_unified
                      # architecture, which transformers gained in 5.10.0;
                      # this pin carries 5.5.4. 5.13.1 is the newest release
                      # whose dependency bounds the pinned tokenizers
                      # (0.22.2), safetensors (0.8.0) and huggingface-hub
                      # (1.16) still satisfy -- 5.14+ wants tokenizers 0.23.
                      transformers = python-prev.transformers.overridePythonAttrs (old: rec {
                        version = "5.13.1";
                        src = final.fetchFromGitHub {
                          owner = "huggingface";
                          repo = "transformers";
                          tag = "v${version}";
                          hash = "sha256-7khrrnATvSl7Wo8yvsZ2Shyzv6saXUkcs8lvF23Fbe4=";
                        };
                      });
                      # transformers and accelerate depend on `torch`, while
                      # the shell below asks for `torch-bin`. Without this,
                      # buildEnv is handed two different torch-2.12.0 store
                      # paths and refuses with a conflicting-subpath error.
                      # Aliasing collapses the set onto one derivation.
                      torch = python-final.torch-bin;
                      torchvision = python-final.torchvision-bin;
                      # bitsandbytes reads its CUDA settings off `torch`
                      # (cudaSupport ? torch.cudaSupport), which the wheel
                      # build above does not carry. `.override` on the stock
                      # package is too late -- the python builder's disabled
                      # check forces the derivation, and with it that
                      # default, before the override lands -- so the file is
                      # called directly with the settings given outright.
                      # CUDA 13 moved the crt/ headers (crt/host_config.h,
                      # crt/host_defines.h) out of cuda_nvcc into cuda_crt,
                      # which the derivation knows nothing about: nvcc's
                      # compiler test wants them via CUDA_HOME, and the
                      # host-side pythonInterface.cpp via the include path,
                      # so cuda_crt goes into both. The
                      # kernels compile for the 8.9 capability set at the
                      # top; doCheck is off because the test suite wants a
                      # GPU the sandbox does not have.
                      bitsandbytes =
                        (python-final.callPackage
                          (nixpkgs + "/pkgs/development/python-modules/bitsandbytes") {
                            cudaSupport = true;
                            cudaPackages = final.cudaPackages;
                            rocmSupport = false;
                          }).overridePythonAttrs (old:
                          let
                            cudaHome = final.symlinkJoin {
                              name = "cuda-native-redist-crt";
                              paths = [ old.env.CUDA_HOME final.cudaPackages.cuda_crt ];
                            };
                          in
                          {
                            doCheck = false;
                            buildInputs = (old.buildInputs or [ ]) ++ [
                              final.cudaPackages.cuda_crt
                            ];
                            env = old.env // {
                              CUDA_HOME = cudaHome;
                              NVCC_PREPEND_FLAGS = "-I${cudaHome}/include -L${cudaHome}/lib";
                            };
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
        tabfm = pkgs.python3Packages.buildPythonPackage {
          pname = "tabfm";
          version = "1.0.0-unstable-2026-08-18";
          pyproject = true;
          src = pkgs.fetchFromGitHub {
            owner = "google-research";
            repo = "tabfm";
            rev = "fbb665569425fd2f490c6576b3af967876fe11ff";
            hash = "sha256-yhzAzSSD3A9OWihA/bU+mkk29uAjbsrIO121HgMVFZw=";
          };
          build-system = [ pkgs.python3Packages.flit-core ];
          postPatch = ''
            substituteInPlace tabfm/src/classifier_and_regressor.py \
              --replace-fail 'jt.typed = jt.jaxtyped(typechecker=typeguard.typechecked)' \
                             'jt.typed = lambda function: function'
          '';
          dependencies = with pkgs.python3Packages; [
            absl-py
            huggingface-hub
            jaxtyping
            numpy
            pandas
            scikit-learn
            scipy
            torch-bin
            typeguard
          ];
          pythonRemoveDeps = [
            "jaxtyping"
            "typeguard"
          ];
          doCheck = false;
        };
        root = if isLaptop then pkgs.root else utils.packages.${system}.rootCuda;
        talys = talys-nix.packages.${system}.default;
        talys-potentials = talys-nix.packages.${system}.talys-atomki-v2-potentials;
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
              clang-tools
              # API docs: scripts/build_docs.sh
              doxygen
              graphviz
            ];
            buildInputs = [
              analysis-utils
              root
              talys
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
                  # VLM event classification (python/vlm.py). transformers 5.5
                  # is the first release carrying gemma4; this nixpkgs pin has
                  # 5.5.4, so the model loads with no overlay.
                  transformers
                  accelerate
                  safetensors
                  pillow
                  # AutoProcessor for gemma4 instantiates Gemma4VideoProcessor
                  # even for a still image, and that hard-requires torchvision.
                  torchvision
                  # 8-bit / 4-bit weights for the 12B rung of the ladder
                  # (config.VLM_LOAD_IN); bf16 12B does not fit a 24 GB card.
                  bitsandbytes
                  peft
                  tabfm
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

              # doxygen-awesome-css ships only data files, so it cannot be
              # found on PATH; scripts/build_docs.sh reads this.
              export DOXYGEN_AWESOME_CSS="${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css"

              # --- dataset selection ---
              export MUSIC_DATASET="${dataset}"
              export MUSIC_DATASET_DIR="$git_root/analysis/${dataset}"
              echo "MUSIC dataset: ${dataset}  ($MUSIC_DATASET_DIR)"

              # Where GENERATED outputs (root_files, plots) land. Defaults to the
              # in-repo dataset dir (on /home); override this one var to redirect
              # processed output to a scratch drive without rebuilding.
              export MUSIC_RESULTS_DIR="''${MUSIC_RESULTS_DIR:-$git_root/analysis/${dataset}}"
              echo "MUSIC results: $MUSIC_RESULTS_DIR"
              export TALYS_BIN="${talys}/bin/talys"
              export TALYS_ATOMKI_V2_DIR="${talys-potentials}/share/talys/atomki-v2"

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

              # for using argo
              anl-opencode() {
                if ! (exec 3<>"/dev/tcp/127.0.0.1/8118") 2>/dev/null; then
                  echo "anl-opencode: nothing on 127.0.0.1:8118 - start scripts/anl-tunnel.sh first." >&2
                  return 1
                fi
                HTTP_PROXY="http://127.0.0.1:8118" HTTPS_PROXY="http://127.0.0.1:8118" \
                  NO_PROXY="localhost,127.0.0.1,::1" opencode "$@"
              }
            '';
          };
      in
      let
        mkPackage =
          dataset:
          let
            extraBuild = if (!isLaptop) then ''GPU_LIB_OUT="$out/lib/libgpuaccel.so"'' else "";
            extraInstall = if (!isLaptop) then "cp tooling/gpu/libgpuaccel.so $out/lib/" else "";
            # The Makefile derives this from `git rev-parse`, which cannot work
            # in the sandbox: the flake source carries no .git and git is not a
            # build input, so it always fell back to "unknown". The flake knows
            # the revision, so pass it in and let the command-line value win.
            gitHash =
              if self ? shortRev then
                self.shortRev
              else if self ? dirtyShortRev then
                self.dirtyShortRev
              else
                "unknown";
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
              make -C tooling -j GIT_HASH="${gitHash}" DATASET_DIR_OUT="$out/analysis/${dataset}" ASSETS_DIR_OUT="$out/assets" TALYS_BIN="${talys}/bin/talys" TALYS_ATOMKI_V2_DIR="${talys-potentials}/share/talys/atomki-v2" ${extraBuild}
            '';

            installPhase = ''
              mkdir -p $out/bin $out/lib $out/assets $out/analysis/${dataset}/config
              cp -r tooling/assets/. $out/assets/
              cp analysis/${dataset}/bin/* $out/bin/
              ${extraInstall}
              cp -r analysis/${dataset}/config/* $out/analysis/${dataset}/config/
            '';
          };
        # The API reference. Dataset-independent: it documents tooling/, which
        # is identical across datasets, so there is one docs output rather than
        # one per dataset.
        docs = pkgs.stdenv.mkDerivation {
          name = "music-docs";
          src = ./.;

          nativeBuildInputs = with pkgs; [
            doxygen
            graphviz
            git
          ];

          dontConfigure = true;

          buildPhase = ''
            runHook preBuild
            export DOXYGEN_AWESOME_CSS="${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css"
            bash scripts/build_docs.sh "$PWD/docs"

            # A warning here means a broken \ref, a malformed doc block, or a
            # parameter that no longer exists. Failing the build is what keeps
            # the published site honest, and gives CI the check for free.
            if [ -s docs/doxygen-warnings.log ]; then
              echo "doxygen emitted warnings:" >&2
              cat docs/doxygen-warnings.log >&2
              exit 1
            fi
            runHook postBuild
          '';

          installPhase = ''
            runHook preInstall
            mkdir -p $out
            cp -r docs/html/. $out/
            runHook postInstall
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
          inherit docs;
        };
      }
    );
}
