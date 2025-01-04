{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-fhs",
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in
  qlat-pkgs."${name}".env
