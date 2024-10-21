{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  qlat-tests1 = many-qlat-pkgs-all.qlat-tests;
  qlat-tests2 = many-qlat-pkgs-all.qlat-tests-std;
  qlat-tests3 = many-qlat-pkgs-all.qlat-tests-std-clang;
  qlat-tests4 = many-qlat-pkgs-all.qlat-tests-cuda;
  qlat-tests5 = many-qlat-pkgs-all.qlat-tests-std-cuda;
in qlat-tests1
