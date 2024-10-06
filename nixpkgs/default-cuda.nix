let
  qlat-pkgs = import ./qlat-pkgs.nix;
in with qlat-pkgs;
  many-qlat-pkgs-core-w-cuda
