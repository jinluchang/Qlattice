# Used to set up nix shell.
#
# With this file, you can build with:
#
# $ nix-shell
# $ ./build default-nix

{
  pkgs ? import <nixpkgs> {}
}:

let
  local-pkgs = import ./default.nix;
  qlat-pkgs = local-pkgs."qlat-pkgs";
  env = pkgs.mkShell {
    name = "qlat-build-sh";
    packages = qlat-pkgs;
    inputsFrom = qlat-pkgs;
  };
in env
