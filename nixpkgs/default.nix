let
  nixpkgs = <nixpkgs>;
  pkgs = import nixpkgs {
    config = {};
    overlays = [];
  };
  qlat-pkgs = {
    inherit cuba qlat_utils qlat;
  };
  cuba = pkgs.callPackage ./cuba.nix {};
  qlat_utils = pkgs.python3Packages.callPackage ./qlat_utils.nix { inherit cuba; };
  qlat = pkgs.python3Packages.callPackage ./qlat.nix { inherit cuba qlat_utils; };
in qlat-pkgs
