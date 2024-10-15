{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  qlat-py-tests1 = many-qlat-pkgs-all.qlat-py-tests-local;
  qlat-py-tests2 = many-qlat-pkgs-all.qlat-py-tests-std-local;
  qlat-py-tests3 = many-qlat-pkgs-all.qlat-py-tests-std-clang-local;
  qlat-py-tests4 = many-qlat-pkgs-all.qlat-py-tests-cuda-local;
  qlat-py-tests5 = many-qlat-pkgs-all.qlat-py-tests-std-cuda-local;
in qlat-py-tests1
