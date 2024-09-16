let

  nixpkgs = import <nixpkgs>;

  nixgl-src = (nixpkgs {}).fetchFromGitHub {
    owner = "nix-community";
    repo = "nixGL";
    rev = "310f8e49a149e4c9ea52f1adf70cdc768ec53f8a";
    hash = "sha256-lnzZQYG0+EXl/6NkGpyIz+FEOc/DSEG57AP1VsdeNrM=";
  };

  nixgl = import nixgl-src {};

  overlay = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "";
    #
    is-pypi-src = true;
    qlat-cudaSupport = false;
    qlat-ngpu = "2"; # adjust with actual number of GPUs
    qlat-eigen = pkgs.grid-lehner;
    qlat-stdenv = pkgs.stdenv;
    mpi = prev.mpi.override { cudaSupport = pkgs.qlat-cudaSupport; };
    grid-lehner-c-lime = pkgs.qio;
    call-pkg = pkgs.callPackage;
    py-call-pkg = pkgs.python3Packages.callPackage;
    #
    qlat-nixgl = if pkgs.qlat-cudaSupport then nixgl.auto.nixGLDefault else "";
    #
    cuba = call-pkg ./cuba.nix { stdenv = pkgs.qlat-stdenv; };
    qlat_utils = py-call-pkg ./qlat_utils.nix {
      stdenv = pkgs.qlat-stdenv; eigen = pkgs.qlat-eigen;
      cudaSupport = pkgs.qlat-cudaSupport;
    };
    qlat = py-call-pkg ./qlat.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
    };
    qlat_grid = py-call-pkg ./qlat_grid.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
    };
    qlat_cps = py-call-pkg ./qlat_cps.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
    };
    c-lime = call-pkg ./c-lime.nix { stdenv = pkgs.qlat-stdenv; };
    qmp = call-pkg ./qmp.nix { stdenv = pkgs.qlat-stdenv; };
    qio = call-pkg ./qio.nix { stdenv = pkgs.qlat-stdenv; };
    cps = call-pkg ./cps.nix { stdenv = pkgs.qlat-stdenv; };
    grid-lehner = call-pkg ./grid-lehner.nix {
      stdenv = pkgs.qlat-stdenv;
      c-lime = pkgs.grid-lehner-c-lime;
      cudaSupport = pkgs.qlat-cudaSupport;
    };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
    };
    #
    qlat-examples-cpp = py-call-pkg ./qlat-examples-cpp.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-cpp-grid = py-call-pkg ./qlat-examples-cpp-grid.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py = py-call-pkg ./qlat-examples-py.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py-gpt = py-call-pkg ./qlat-examples-py-gpt.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py-cps = py-call-pkg ./qlat-examples-py-cps.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    #
    qlat-dep-pkgs = with pkgs; ([
      git pkg-config zlib gsl fftw fftwFloat hdf5-cpp openssl gmp mpfr
    ]
    ++ (if qlat-cudaSupport then [ qlat-nixgl ] else [])
    );
    #
    qlat-py = pkgs.python3.withPackages (ps: with pkgs; [
      qlat_utils
      qlat
      qlat_cps
      qlat_grid
      gpt-lehner
      #
      qlat-examples-cpp
      qlat-examples-cpp-grid
      qlat-examples-py
      qlat-examples-py-gpt
      qlat-examples-py-cps
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
      #
      qlat-examples-cpp
      qlat-examples-py
    ]);
    qlat-pkgs = with pkgs; [
      mpi cuba qlat-eigen qlat-py
    ] ++ pkgs.qlat-dep-pkgs;
    inherit pkgs;
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
    "all-pkgs${pkgs.qlat-name}" = pkgs.pkgs;
  };

  many-qlat-pkgs-core = {}
  # // mk-qlat-pkgs [ overlay-local ]
  // mk-qlat-pkgs [ overlay-cuda overlay-local ]
  # // mk-qlat-pkgs [ overlay-std overlay-clang overlay-local ]
  ;

  many-qlat-pkgs-all = many-qlat-pkgs-core
  // mk-qlat-pkgs []
  # // mk-qlat-pkgs [ overlay-cuda ]
  # // mk-qlat-pkgs [ overlay-clang ]
  # // mk-qlat-pkgs [ overlay-clang overlay-local ]
  // mk-qlat-pkgs [ overlay-std ]
  // mk-qlat-pkgs [ overlay-std overlay-local ]
  # // mk-qlat-pkgs [ overlay-std overlay-cuda ]
  # // mk-qlat-pkgs [ overlay-std overlay-cuda overlay-local ]
  // mk-qlat-pkgs [ overlay-std overlay-clang ]
  ;

  many-qlat-pkgs = many-qlat-pkgs-core;
  # many-qlat-pkgs = many-qlat-pkgs-all;

in many-qlat-pkgs
