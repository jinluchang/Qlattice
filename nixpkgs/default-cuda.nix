{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "many-q-pkgs-more-w-cuda",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
