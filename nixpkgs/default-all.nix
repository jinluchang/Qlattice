{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in with qlat-pkgs;
  many-qlat-pkgs-all
