let
  nixpkgs = <nixpkgs>;
  pkgs = import nixpkgs {
    config = {};
    overlays = [];
  };
  qlat-pkgs = {
    inherit cuba qlat_utils qlat grid-lehner;
  };
  cuba = pkgs.callPackage ./cuba.nix {};
  qlat_utils = pkgs.python3Packages.callPackage ./qlat_utils.nix { inherit cuba; };
  qlat = pkgs.python3Packages.callPackage ./qlat.nix { inherit cuba qlat_utils; };
  grid-lehner = pkgs.callPackage ./grid-lehner.nix {};
in qlat-pkgs
