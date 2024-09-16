{
  pkgs ? import <nixpkgs> {}
}:

let
  many-qlat-pkgs = import ./default.nix;
  #
  qlat-pkgs = many-qlat-pkgs.qlat-pkgs-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-cuda-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-std-clang-local;
  #
  # qlat-pkgs = let pkgs = many-qlat-pkgs.all-pkgs-cuda-local; in [
  #   # pkgs.qlat-nixgl
  #   (pkgs.python3.withPackages (ps: with pkgs; [
  #     qlat_utils
  #     qlat
  #     gpt-lehner
  #   ]))
  # ];
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
