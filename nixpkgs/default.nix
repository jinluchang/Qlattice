let
  nixpkgs = <nixpkgs>;
  pkgs = import nixpkgs {
    config = {};
    overlays = [];
  };
  qlat-pkgs = {
    inherit qlat_utils qlat;
  };
  qlat_utils = pkgs.python3Packages.callPackage ./qlat_utils.nix {};
  qlat = pkgs.python3Packages.callPackage ./qlat.nix { inherit qlat_utils; };
in qlat-pkgs
