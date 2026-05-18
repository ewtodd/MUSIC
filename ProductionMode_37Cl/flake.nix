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
        analysis-utils = utils.packages.${system}.cuda;
        root = utils.packages.${system}.rootCuda;
        agenixPkg = agenix.packages.${system}.default;
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
            agenixPkg
          ];
          shellHook = ''
            echo "Analysis-Utilities version: ${analysis-utils.version} (CUDA)"
            export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
            flake_root="$PWD"
            git_root="$(git -C "$flake_root" rev-parse --show-toplevel)"

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
            cd macros
          '';
        };
      }
    );
}
