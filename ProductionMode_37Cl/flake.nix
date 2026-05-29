{
  description = "ROOT Analysis Development Environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    agenix = {
      url = "github:ryantm/agenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    utils = {
      url = "/home/e-work/Analysis-Utilities";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      agenix,
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
        root = if !isCUDA then pkgs.root else utils.packages.${system}.rootCuda;
        agenixPkg = agenix.packages.${system}.default;
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
      in
      {
        devShells.default = pkgs.mkShell {
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
            agenixPkg
          ]
          ++ pkgs.lib.optionals isCUDA [
            pkgs.cudaPackages.cuda_nvcc
            pkgs.cudaPackages.cuda_cudart
            pkgs.cudaPackages.cuda_cccl
          ];
          shellHook = ''
             echo "Analysis-Utilities version: ${analysis-utils.version}${pkgs.lib.optionalString isCUDA " (CUDA)"}"
             flake_root="$PWD"
            ${pkgs.lib.optionalString isCUDA ''
              export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
              export LD_LIBRARY_PATH="/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
              install -m 644 ${clangdConfigFile} "$flake_root/gpu/.clangd"
            ''}

            git_root="$(git -C "$flake_root" rev-parse --show-toplevel)"


            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

            mkdir -p "$flake_root/.claude"
            (
              cd "$git_root/secrets"
              ${agenixPkg}/bin/agenix -d settings.json.age \
                -i "$HOME/.ssh/id_ed25519"
            ) > "$flake_root/.claude/settings.json"
            chmod 644 "$flake_root/.claude/settings.json"

            find "$flake_root" -type d \
              \( -name .git -o -name .claude -o -name result -o -name node_modules -o -name .direnv -o -name plots -o -name root_files \) -prune \
              -o -type d -not -path "$flake_root" -print \
            | while read -r dir; do
                if [ ! -e "$dir/.claude" ]; then
                  ln -s "$flake_root/.claude" "$dir/.claude"
                fi
              done
            alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm *_cxx.so *_cxx.d *_cxx_ACLiC_dict_rdict.pcm AutoDict_*'
          '';
        };
      }
    );
}
