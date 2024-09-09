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
    qlat-name = "";
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
    qlat-dep-pkgs = with pkgs; [
      git pkg-config zlib gsl fftw fftwFloat hdf5-cpp openssl gmp mpfr
    ];
    #
    qlat-py = pkgs.python3.withPackages (ps: with pkgs; [
      qlat_utils
      qlat
      qlat_cps
      qlat_grid
      gpt-lehner
    ]);
    qlat-pkgs = with pkgs; [
      mpi cuba qlat-eigen cps qmp qio grid-lehner qlat-py
    ] ++ pkgs.qlat-dep-pkgs;
    #
  };

  overlay-local = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-local";
    is-pypi-src = false;
  };

  overlay-std = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-std";
    qlat-eigen = pkgs.eigen;
    qlat-py = pkgs.python3.withPackages (ps: with pkgs; [
      qlat_utils
      qlat
    ]);
    qlat-pkgs = with pkgs; [
      mpi cuba qlat-eigen qlat-py
    ] ++ pkgs.qlat-dep-pkgs;
  };

  overlay-clang = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-clang";
    qlat-stdenv = pkgs.clangStdenv;
    openmp = pkgs.llvmPackages.openmp;
  };

  overlay-cuda = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-cuda";
    qlat-cudaSupport = true;
  };

  mk-qlat-pkgs = overlays: let
    pkgs = nixpkgs {
      config = {
        allowUnfree = true;
        # cudaSupport = true;
      };
      overlays = [
        overlay
      ] ++ overlays;
    };
  in {
    "qlat-py${pkgs.qlat-name}" = pkgs.qlat-py;
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
  };

  qlat-pkgs = {}
  // mk-qlat-pkgs []
  // mk-qlat-pkgs [ overlay-local ]
  // mk-qlat-pkgs [ overlay-cuda ]
  // mk-qlat-pkgs [ overlay-cuda overlay-local ]
  // mk-qlat-pkgs [ overlay-clang ]
  // mk-qlat-pkgs [ overlay-clang overlay-local ]
  // mk-qlat-pkgs [ overlay-clang overlay-cuda ]
  // mk-qlat-pkgs [ overlay-clang overlay-cuda overlay-local ]
  // mk-qlat-pkgs [ overlay-std ]
  // mk-qlat-pkgs [ overlay-std overlay-local ]
  // mk-qlat-pkgs [ overlay-std overlay-cuda ]
  // mk-qlat-pkgs [ overlay-std overlay-cuda overlay-local ]
  // mk-qlat-pkgs [ overlay-std overlay-clang ]
  // mk-qlat-pkgs [ overlay-std overlay-clang overlay-local ]
  // mk-qlat-pkgs [ overlay-std overlay-clang overlay-cuda ]
  // mk-qlat-pkgs [ overlay-std overlay-clang overlay-cuda overlay-local ]
  ;

in qlat-pkgs
