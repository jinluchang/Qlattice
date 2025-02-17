{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "q-pkgs-more-w-cuda",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
