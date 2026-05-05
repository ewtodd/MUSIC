{
  description = "ROOT Analysis Development Environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    utils = {
      inputs.nixpkgs.follows = "nixpkgs";
      url = "github:ewtodd/Analysis-Utilities";
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
        pkgs = nixpkgs.legacyPackages.${system};
        analysis-utils = utils.packages.${system}.default;
      in
      {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
            clang-tools
          ];
          buildInputs = with pkgs; [
            analysis-utils
            root
            bash
          ];
          shellHook = ''
            echo "Analysis-Utilities version: ${analysis-utils.version}"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };
      }
    );
}
