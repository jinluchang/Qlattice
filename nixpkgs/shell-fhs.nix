{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-jhub-fhs",
}:

let
  q-pkgs = import ./q-pkgs.nix { inherit nixpkgs; };
in
  q-pkgs."${name}".env
