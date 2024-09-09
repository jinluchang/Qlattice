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
in pkgs.mkShell {
  name = "qlat-build-sh";
  packages = with pkgs; [
    local-pkgs.qlat-full-cuda-local
  ];
}
