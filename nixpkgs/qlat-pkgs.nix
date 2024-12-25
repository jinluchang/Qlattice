{
  nixpkgs ? import ./nixpkgs.nix,
  ngpu ? "2", # adjust with actual number of GPUs
  nvcc-arch ? "sm_86", # adjust with actual arch of Nvidia GPU
}:

let

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
    is-pypi-src = false;
    qlat-cudaSupport = false;
    qlat-ngpu = ngpu;
    qlat-nvcc-arch = nvcc-arch;
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
      stdenv = pkgs.qlat-stdenv;
      eigen = pkgs.qlat-eigen;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
    };
    qlat = py-call-pkg ./qlat.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
    };
    qlat_grid = py-call-pkg ./qlat_grid.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
    };
    qlat_cps = py-call-pkg ./qlat_cps.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
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
      nvcc-arch = pkgs.qlat-nvcc-arch;
    };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
    };
    #
    qlat-examples-cpp = py-call-pkg ./qlat-examples-cpp.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-cpp-grid = py-call-pkg ./qlat-examples-cpp-grid.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py = py-call-pkg ./qlat-examples-py.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py-gpt = py-call-pkg ./qlat-examples-py-gpt.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-examples-py-cps = py-call-pkg ./qlat-examples-py-cps.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
      ngpu = pkgs.qlat-ngpu;
    };
    qlat-pypi = py-call-pkg ./qlat-pypi.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
    };
    #
    qlat-dep-pkgs = with pkgs; ([
      git pkg-config zlib gsl fftw fftwFloat hdf5-cpp openssl gmp mpfr
    ]
    ++ (if qlat-cudaSupport then [ qlat-nixgl ] else [])
    );
    #
    qlat-dep-pkgs-extra = with pkgs; [
      mpi cuba qlat-eigen cps qmp qio grid-lehner
    ];
    qlat-py-pkgs = with pkgs; [
      qlat_utils
      qlat
      qlat_cps
      qlat_grid
      gpt-lehner
    ];
    qlat-tests-pkgs = with pkgs; [
      qlat-pypi
      qlat-examples-cpp
      qlat-examples-cpp-grid
      qlat-examples-py
      qlat-examples-py-gpt
      qlat-examples-py-cps
    ];
    #
    qlat-py = pkgs.python3.withPackages (ps: [
      ps.build
      ps.wheel
    ] ++ pkgs.qlat-py-pkgs);
    qlat-pkgs = with pkgs; [
      qlat-py
    ] ++ pkgs.qlat-dep-pkgs ++ pkgs.qlat-dep-pkgs-extra;
    qlat-tests = pkgs.buildEnv {
      name = "qlat-tests${pkgs.qlat-name}";
      paths = pkgs.qlat-tests-pkgs;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-env = pkgs.buildEnv {
      name = "qlat-env${pkgs.qlat-name}";
      paths = pkgs.qlat-pkgs;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-sh = pkgs.mkShell rec {
      name = "qlat-sh${pkgs.qlat-name}";
      packages = [ pkgs.qlat-env ];
      inputsFrom = packages;
    };
    #
    qlat-jhub-py = pkgs.python3.withPackages (ps: with ps; [
      ipykernel
      pip
      numpy
      scipy
      sympy
      jax
      jaxlib
      meson
      ninja
      mpi4py
      psutil
      cython
      pybind11
      pythran
      poetry-core
      pkgconfig
      meson-python
      scikit-build
      setuptools-scm
      pyproject-metadata
      build
      wheel
      pyproject-hooks
      pep517
      packaging
      tomli
      flit-core
      virtualenv
      h5py
      pandas
      scikit-learn
      xarray
      matplotlib
      plotly
      seaborn
      jupyter-server-mathjax
      numba
      transformers
      torch
      sphinx
      linkify-it-py
      myst-parser
      pycuda
      pytools
      lz4
      torchvision
      torchaudio
      xformers
      jupyterlab
      jupyterhub
      jupyterhub-systemdspawner
    ] ++ pkgs.qlat-py-pkgs);
    qlat-jhub-env = pkgs.buildEnv {
      name = "qlat-jhub-env${pkgs.qlat-name}";
      paths = with pkgs; [
        qlat-jhub-py
        bashInteractive
        coreutils
        openssh
        linux-pam
        findutils
        gcc
        clang-tools
        git
        gnumake
        zlib
        pkg-config
        mpi
        killall
        wget
        rsync
        automake
        autoconf
        gsl
        fftw
        fftwFloat
        openssl
        gnuplot
        texliveFull
        pipx
        twine
        poppler_utils
        file
        zip
        unzip
      ] ++ pkgs.qlat-dep-pkgs ++ pkgs.qlat-dep-pkgs-extra;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-jhub-sh = pkgs.mkShell rec {
      name = "qlat-jhub-sh${pkgs.qlat-name}";
      packages = [ pkgs.qlat-jhub-env ];
      inputsFrom = packages;
    };
    #
  };

  overlay-pypi = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-pypi";
    is-pypi-src = true;
  };

  overlay-std = final: prev: let
    pkgs = final;
  in rec {
    qlat-name = "${prev.qlat-name}-std";
    qlat-eigen = pkgs.eigen;
    qlat-dep-pkgs-extra = with pkgs; [
      mpi cuba qlat-eigen
    ];
    qlat-py-pkgs = with pkgs; [
      qlat_utils
      qlat
    ];
    qlat-tests-pkgs = with pkgs; [
      qlat-examples-cpp
      qlat-examples-py
    ];
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

  mk-qlat-pkgs-gen = config: overlays: let
    pkgs = nixpkgs {
      config = {
        allowUnfree = true;
      } // config;
      overlays = [
        overlay
      ] ++ overlays;
    };
  in {
    "qlat-py${pkgs.qlat-name}" = pkgs.qlat-py;
    "qlat-env${pkgs.qlat-name}" = pkgs.qlat-env;
    "qlat-sh${pkgs.qlat-name}" = pkgs.qlat-sh;
    "qlat-tests${pkgs.qlat-name}" = pkgs.qlat-tests;
    "qlat-jhub-py${pkgs.qlat-name}" = pkgs.qlat-jhub-py;
    "qlat-jhub-env${pkgs.qlat-name}" = pkgs.qlat-jhub-env;
    "qlat-jhub-sh${pkgs.qlat-name}" = pkgs.qlat-jhub-sh;
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
    "pkgs${pkgs.qlat-name}" = pkgs;
  };

  mk-qlat-pkgs = mk-qlat-pkgs-gen {};
  mk-qlat-pkgs-cuda = mk-qlat-pkgs-gen { cudaSupport = false; }; # cudaSupport does not compile yet.

  many-qlat-pkgs-core = {}
  // mk-qlat-pkgs []
  // mk-qlat-pkgs [ overlay-std overlay-clang ]
  ;

  many-qlat-pkgs-core-w-cuda = many-qlat-pkgs-core
  // mk-qlat-pkgs-cuda [ overlay-cuda ]
  ;

  many-qlat-pkgs-all = many-qlat-pkgs-core-w-cuda
  // mk-qlat-pkgs [ overlay-pypi ]
  // mk-qlat-pkgs-cuda [ overlay-cuda overlay-pypi ]
  // mk-qlat-pkgs [ overlay-clang ]
  // mk-qlat-pkgs [ overlay-clang overlay-pypi ]
  // mk-qlat-pkgs [ overlay-std ]
  // mk-qlat-pkgs [ overlay-std overlay-pypi ]
  // mk-qlat-pkgs-cuda [ overlay-std overlay-cuda ]
  // mk-qlat-pkgs-cuda [ overlay-std overlay-cuda overlay-pypi ]
  // mk-qlat-pkgs [ overlay-std overlay-clang overlay-pypi ]
  ;

in {
  #
  inherit mk-qlat-pkgs-gen;
  inherit mk-qlat-pkgs;
  inherit mk-qlat-pkgs-cuda;
  inherit overlay;
  inherit overlay-pypi;
  inherit overlay-std;
  inherit overlay-clang;
  inherit overlay-cuda;
  #
  inherit many-qlat-pkgs-core;
  inherit many-qlat-pkgs-core-w-cuda;
  inherit many-qlat-pkgs-all;
  #
}
