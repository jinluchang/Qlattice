{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  env1 = many-qlat-pkgs-all.qlat-env;
  env2 = many-qlat-pkgs-all.qlat-env-std;
  env3 = many-qlat-pkgs-all.qlat-env-std-clang;
in env1
