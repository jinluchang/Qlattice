{
  nixpkgs ? null,
  version ? null,
  use-gitee ? null,
  name ? "q-pkgs-cuda",
}:

let
  all-q-pkgs = import ./q-pkgs.nix { inherit nixpkgs version use-gitee; };
in
  all-q-pkgs."${name}"
