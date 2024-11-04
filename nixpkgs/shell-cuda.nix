{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  sh1 = many-qlat-pkgs-all.qlat-sh-cuda;
  sh2 = many-qlat-pkgs-all.qlat-sh-std-cuda;
in sh1
