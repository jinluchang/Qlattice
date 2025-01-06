{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "many-qlat-pkgs",
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in
  qlat-pkgs."${name}"
