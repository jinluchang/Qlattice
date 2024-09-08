let
  nixpkgs = <nixpkgs>;
  pkgs = import nixpkgs {
    config = {};
    overlays = [];
  };
  qlat-pkgs = {
    inherit cuba qlat_utils qlat;
    inherit c-lime grid-lehner gpt-lehner qlat_grid;
    inherit qmp qio cps qlat_cps;
  };
  #
  eigen = grid-lehner;
  # eigen = pkgs.eigen;
  #
  cuba = pkgs.callPackage ./cuba.nix {};
  qlat_utils = pkgs.python3Packages.callPackage ./qlat_utils.nix { inherit cuba eigen; };
  qlat = pkgs.python3Packages.callPackage ./qlat.nix { inherit cuba eigen qlat_utils; };
  qlat_grid = pkgs.python3Packages.callPackage ./qlat_grid.nix { inherit cuba qlat_utils qlat c-lime grid-lehner; };
  qlat_cps = pkgs.python3Packages.callPackage ./qlat_cps.nix { inherit cuba qlat_utils qlat c-lime qmp qio cps; };
  c-lime = pkgs.callPackage ./c-lime.nix {};
  qmp = pkgs.callPackage ./qmp.nix {};
  qio = pkgs.callPackage ./qio.nix { inherit qmp; };
  cps = pkgs.callPackage ./cps.nix { inherit qmp qio; };
  grid-lehner = pkgs.callPackage ./grid-lehner.nix { inherit c-lime; };
  gpt-lehner = pkgs.python3Packages.callPackage ./gpt-lehner.nix { inherit c-lime grid-lehner; };
in qlat-pkgs
