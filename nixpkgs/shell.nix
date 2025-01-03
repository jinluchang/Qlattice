{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-sh",
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in
  qlat-pkgs."${name}"
