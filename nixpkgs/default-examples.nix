{
  nixpkgs ? import ./nixpkgs.nix
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  tests1 = many-qlat-pkgs-all.qlat-tests;
  tests2 = many-qlat-pkgs-all.qlat-tests-std;
  tests3 = many-qlat-pkgs-all.qlat-tests-std-clang;
  tests4 = many-qlat-pkgs-all.qlat-tests-cuda;
  tests5 = many-qlat-pkgs-all.qlat-tests-std-cuda;
in tests1
