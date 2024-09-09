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
    #
    is-pypi-src = true;
    qlat-cudaSupport = false;
    qlat-eigen = pkgs.grid-lehner;
    qlat-stdenv = pkgs.stdenv;
    mpi = prev.mpi.override { cudaSupport = pkgs.qlat-cudaSupport; };
    grid-lehner-c-lime = pkgs.qio;
    call-pkg = pkgs.callPackage;
    py-call-pkg = pkgs.python3Packages.callPackage;
    #
    cuba = call-pkg ./cuba.nix { stdenv = pkgs.qlat-stdenv; };
    qlat_utils = py-call-pkg ./qlat_utils.nix { stdenv = pkgs.qlat-stdenv; eigen = pkgs.qlat-eigen; };
    qlat = py-call-pkg ./qlat.nix { stdenv = pkgs.qlat-stdenv; };
    qlat_grid = py-call-pkg ./qlat_grid.nix { stdenv = pkgs.qlat-stdenv; };
    qlat_cps = py-call-pkg ./qlat_cps.nix { stdenv = pkgs.qlat-stdenv; };
    c-lime = call-pkg ./c-lime.nix { stdenv = pkgs.qlat-stdenv; };
    qmp = call-pkg ./qmp.nix { stdenv = pkgs.qlat-stdenv; };
    qio = call-pkg ./qio.nix { stdenv = pkgs.qlat-stdenv; };
    cps = call-pkg ./cps.nix { stdenv = pkgs.qlat-stdenv; };
    grid-lehner = call-pkg ./grid-lehner.nix { stdenv = pkgs.qlat-stdenv; c-lime = pkgs.grid-lehner-c-lime; };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix { stdenv = pkgs.qlat-stdenv; };
    #
    qlat-std = pkgs.python3.withPackages (ps: with ps; [
      qlat_utils
      qlat
    ]);
    #
    qlat-full = pkgs.python3.withPackages (ps: with ps; [
      qlat_utils
      qlat
      qlat_cps
      qlat_grid
      gpt-lehner
    ]);
    #
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
    qlat-stdenv = pkgs.clangStdenv;
    openmp = pkgs.llvmPackages.openmp;
  };

  overlay-cuda = final: prev: let
    pkgs = final;
  in rec {
    qlat-cudaSupport = true;
  };

  mk-qlat-full = overlays: let
    pkgs = nixpkgs {
      config = {
        allowUnfree = true;
        # cudaSupport = true;
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
        # cudaSupport = true;
      };
      overlays = [
        overlay
        overlay-std
      ] ++ overlays;
    };
  in pkgs.qlat-std;

  qlat-pkgs = {
    qlat-full = mk-qlat-full [];
    qlat-full-clang = mk-qlat-full [ overlay-clang ];
    qlat-full-cuda = mk-qlat-full [ overlay-cuda ];
    qlat-full-local = mk-qlat-full [ overlay-local ];
    qlat-full-clang-local = mk-qlat-full [ overlay-clang overlay-local ];
    qlat-full-cuda-local = mk-qlat-full [ overlay-cuda overlay-local ];
    qlat-std = mk-qlat-std [];
    qlat-std-clang = mk-qlat-std [ overlay-clang ];
    qlat-std-cuda = mk-qlat-std [ overlay-cuda ];
    qlat-std-local = mk-qlat-std [ overlay-local ];
    qlat-std-clang-local = mk-qlat-std [ overlay-clang overlay-local ];
    qlat-std-cuda-local = mk-qlat-std [ overlay-cuda overlay-local ];
  };

in qlat-pkgs
