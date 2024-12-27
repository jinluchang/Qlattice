{
  nixpkgs ? import ./nixpkgs.nix
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  pkgs = many-qlat-pkgs-all.pkgs;
in (pkgs.buildFHSEnv {
  name = "qlat-build-env";
  targetPkgs = pkgs: (with pkgs; [
    qlat-jhub-env
  ]);
  runScript = "bash";
  extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
}).env
