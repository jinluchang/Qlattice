{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-jhub-sh",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
