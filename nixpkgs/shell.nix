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
  many-qlat-pkgs = import ./default.nix;
  #
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-cuda-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-std-clang-local;
  #
  qlat-pkgs = let pkgs = many-qlat-pkgs.all-pkgs-cuda-local; in [
    (pkgs.python3.withPackages (ps: with pkgs; [
      qlat_utils
    ]))
  ];
  # qlat-pkgs = let pkgs = many-qlat-pkgs.all-pkgs-local; in [
  #   (pkgs.python3.withPackages (ps: with pkgs; [
  #     qlat_utils
  #   ]))
  # ];
  env = pkgs.mkShell rec {
    name = "qlat-build-sh";
    packages = qlat-pkgs;
    inputsFrom = packages;
  };
in env
