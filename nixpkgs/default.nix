{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "many-q-pkgs",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
