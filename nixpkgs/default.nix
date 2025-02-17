{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "q-pkgs",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
