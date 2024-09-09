let

  nixpkgs = import <nixpkgs>;

  pkgs = nixpkgs {
    config = {
      allowUnfree = true;
    };
    overlays = [
      overlay
    ];
  };

  overlay = final: prev: let
    pkgs = final;
  in rec {
    is-pypi-src = true;
    qlat-eigen = pkgs.grid-lehner;
    grid-lehner-c-lime = pkgs.qio;
    qlat-stdenv = pkgs.stdenv;
    mpi = prev.mpi.override { cudaSupport = true; };
    call-pkg = pkgs.callPackage;
    py-call-pkg = (pkgs.python3Packages.override { stdenv = pkgs.qlat-stdenv; });
    #
    cuba = call-pkg ./cuba.nix { stdenv = pkgs.qlat-stdenv; };
    qlat_utils = py-call-pkg ./qlat_utils.nix { eigen = pkgs.qlat-eigen; };
    qlat = py-call-pkg ./qlat.nix {};
    qlat_grid = py-call-pkg ./qlat_grid.nix {};
    qlat_cps = py-call-pkg ./qlat_cps.nix {};
    c-lime = call-pkg ./c-lime.nix { stdenv = pkgs.qlat-stdenv; };
    qmp = call-pkg ./qmp.nix { stdenv = pkgs.qlat-stdenv; };
    qio = call-pkg ./qio.nix { stdenv = pkgs.qlat-stdenv; };
    cps = call-pkg ./cps.nix { stdenv = pkgs.qlat-stdenv; };
    grid-lehner = call-pkg ./grid-lehner.nix { stdenv = pkgs.qlat-stdenv; c-lime = pkgs.grid-lehner-c-lime; };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix {};
    #
    qlat-std = pkgs.buildEnv {
      name = "qlat-std";
      paths = with pkgs ; [
        cuba
        (pkgs.python3.withPackages (ps: with ps; [
          qlat_utils
          qlat
        ]))
      ];
      extraOutputsToInstall = [ "man" "doc" "info" "dev" "static" ];
    };
    #
    qlat-full = pkgs.buildEnv {
      name = "qlat-full";
      paths = with pkgs ; [
        grid-lehner
        cps
        cuba
        qmp
        qio
        (pkgs.python3.withPackages (ps: with ps; [
          qlat_utils
          qlat
          qlat_cps
          qlat_grid
          gpt-lehner
        ]))
      ];
      extraOutputsToInstall = [ "man" "doc" "info" "dev" "static" ];
    };
  };

  overlay-local = final: prev: let
    pkgs = final;
  in rec {
    is-pypi-src = false;
  };

  overlay-std = final: prev: let
    pkgs = final;
  in rec {
    qlat-eigen = pkgs.eigen;
  };

  overlay-clang = final: prev: let
    pkgs = final;
  in rec {
    qlat-stdenv = prev.clangStdenv;
  };

  mk-qlat-full = overlays: let
    pkgs = nixpkgs {
      config = {
        allowUnfree = true;
        cudaSupport = false;
      };
      overlays = [
        overlay
      ] ++ overlays;
    };
  in pkgs.qlat-full;

  mk-qlat-std = overlays: let
    pkgs = nixpkgs {
      config = {
        allowUnfree = true;
        cudaSupport = false;
      };
      overlays = [
        overlay
        overlay-std
      ] ++ overlays;
    };
  in pkgs.qlat-std;

  qlat-pkgs = {
    qlat-full = mk-qlat-full [];
    qlat-full-local = mk-qlat-full [ overlay-local ];
    qlat-std = mk-qlat-std [];
    qlat-std-local = mk-qlat-std [ overlay-local ];
    qlat-std-clang-local = mk-qlat-std [ overlay-local overlay-clang ];
  };

in qlat-pkgs
