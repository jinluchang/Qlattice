{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  #
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-cuda-local;
  # qlat-pkgs = many-qlat-pkgs.qlat-pkgs-std-clang-local;
  #
  # qlat-pkgs = let pkgs = many-qlat-pkgs.all-pkgs-cuda-local; in [
  #   # pkgs.qlat-nixgl
  #   (pkgs.python3.withPackages (ps: with pkgs; [
  #     qlat_utils
  #     qlat
  #     qlat_cps
  #     qlat_grid
  #     gpt-lehner
  #   ]))
  # ];
  #
  # qlat-pkgs = let pkgs = many-qlat-pkgs.all-pkgs-local; in [
  #   (pkgs.python3.withPackages (ps: with pkgs; [
  #     qlat_utils
  #   ]))
  # ];
   #env = pkgs.mkShell rec {
  #  name = "qlat-build-sh";
  #  packages = qlat-pkgs;
  #  inputsFrom = packages;
  #};
  env1 = many-qlat-pkgs-all.qlat-env-local;
  env2 = many-qlat-pkgs-all.qlat-env-cuda-local;
  #
in env2
