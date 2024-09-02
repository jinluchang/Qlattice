let
  nixpkgs = <nixpkgs>;
  pkgs = import nixpkgs {
    config = {};
    overlays = [];
  };
  qlat-pkgs = {
    inherit cuba qlat_utils qlat c-lime grid-lehner;
  };
  cuba = pkgs.callPackage ./cuba.nix {};
  qlat_utils = pkgs.python3Packages.callPackage ./qlat_utils.nix { inherit cuba; };
  qlat = pkgs.python3Packages.callPackage ./qlat.nix { inherit cuba qlat_utils; };
  c-lime = pkgs.callPackage ./c-lime.nix {};
  grid-lehner = pkgs.callPackage ./grid-lehner.nix { inherit c-lime; };
in qlat-pkgs
