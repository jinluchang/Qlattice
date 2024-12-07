{
  nixpkgs ? import ./nixpkgs.nix
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  sh1 = many-qlat-pkgs-all.qlat-sh;
  sh2 = many-qlat-pkgs-all.qlat-sh-std;
  sh3 = many-qlat-pkgs-all.qlat-sh-std-clang;
in sh1
