{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "many-qlat-pkgs-more-w-cuda",
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in
  qlat-pkgs."${name}"
