{
  nixpkgs ? import ./nixpkgs.nix,
  name ? "qlat-jhub-sh-cuda",
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
in
  qlat-pkgs."${name}"
