{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-jhub-sh-cuda",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}"
