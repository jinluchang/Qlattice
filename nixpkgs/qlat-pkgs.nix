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
    call-pkg = pkgs.callPackage;
    py-call-pkg = pkgs.python3.pkgs.callPackage;
  in {
    qlat-name = "";
    #
    is-pypi-src = false;
    qlat-cudaSupport = false;
    qlat-ngpu = ngpu;
    qlat-nvcc-arch = nvcc-arch;
    qlat-eigen = pkgs.grid-lehner;
    qlat-stdenv = pkgs.stdenv;
    qlat-cc = [ pkgs.gcc ];
    ucx = prev.ucx.override { enableCuda = pkgs.qlat-cudaSupport; };
    ucx-dev = pkgs.buildEnv {
      name = "ucx-dev";
      paths = [ pkgs.ucx ];
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    mpi = (prev.mpi.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags ++ (
        let
          cudaSupport = pkgs.qlat-cudaSupport;
          lib = pkgs.lib;
          cudaPackages = pkgs.cudaPackages;
          ucx-dev = pkgs.ucx-dev;
        in [
          "--with-ucx=${lib.getDev ucx-dev}"
          (lib.withFeatureAs cudaSupport "cuda-libdir" "${cudaPackages.cuda_cudart.stubs}/lib")
        ]);
    })).override { cudaSupport = pkgs.qlat-cudaSupport; };
    python3 = prev.python3.override {
      packageOverrides = final: prev: {
        mpi4py = prev.mpi4py.overridePythonAttrs (prev: {
          doCheck = true;
          nativeBuildInputs = (prev.nativeBuildInputs or [])
          ++ pkgs.lib.optionals pkgs.qlat-cudaSupport [
            pkgs.qlat-nixgl
            pkgs.which
          ];
          preInstallCheck = pkgs.lib.optionalString pkgs.qlat-cudaSupport ''
            which nixGL
            echo
            echo "run with nixGL"
            echo
            cat $(which nixGL) | grep -v 'exec ' | grep -v '^#!' > nix-gl.sh
            echo
            echo cat nix-gl.sh
            cat nix-gl.sh
            source nix-gl.sh
            echo
            echo $LD_LIBRARY_PATH
          '';
        });
      };
    };
    #
    grid-lehner-c-lime = pkgs.qio;
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
    qlat-docs = py-call-pkg ./qlat-docs.nix {
      stdenv = pkgs.qlat-stdenv;
      cudaSupport = pkgs.qlat-cudaSupport;
      nvcc-arch = pkgs.qlat-nvcc-arch;
      nixgl = pkgs.qlat-nixgl;
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
      qlat_grid
      qlat_cps
      gpt-lehner
      qlat-docs
      qlat-pypi
    ];
    qlat-tests-pkgs = with pkgs; [
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
    qlat-fhs = pkgs.buildFHSEnv {
      name = "qlat-fhs${pkgs.qlat-name}";
      targetPkgs = pkgs: [ pkgs.qlat-env ];
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
      profile=''
        PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
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
        clang-tools
        git
        gnumake
        zlib
        pkg-config
        mpi
        killall
        wget
        which
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
      ] ++ pkgs.qlat-cc ++ pkgs.qlat-dep-pkgs ++ pkgs.qlat-dep-pkgs-extra;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-jhub-sh = pkgs.mkShell rec {
      name = "qlat-jhub-sh${pkgs.qlat-name}";
      packages = [ pkgs.qlat-jhub-env ];
      inputsFrom = packages;
    };
    qlat-jhub-fhs = pkgs.buildFHSEnv {
      name = "qlat-jhub-fhs${pkgs.qlat-name}";
      targetPkgs = pkgs: [ pkgs.qlat-jhub-env ];
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
      profile=''
        PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
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
    qlat-cc = [ pkgs.clang openmp ];
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
    "qlat_utils${pkgs.qlat-name}" = pkgs.qlat_utils;
    "qlat${pkgs.qlat-name}" = pkgs.qlat;
    "qlat_grid${pkgs.qlat-name}" = pkgs.qlat_grid;
    "qlat_cps${pkgs.qlat-name}" = pkgs.qlat_cps;
    "qlat-py${pkgs.qlat-name}" = pkgs.qlat-py;
    "qlat-tests${pkgs.qlat-name}" = pkgs.qlat-tests;
    "qlat-env${pkgs.qlat-name}" = pkgs.qlat-env;
    "qlat-sh${pkgs.qlat-name}" = pkgs.qlat-sh;
    "qlat-fhs${pkgs.qlat-name}" = pkgs.qlat-fhs;
    "qlat-jhub-py${pkgs.qlat-name}" = pkgs.qlat-jhub-py;
    "qlat-jhub-env${pkgs.qlat-name}" = pkgs.qlat-jhub-env;
    "qlat-jhub-sh${pkgs.qlat-name}" = pkgs.qlat-jhub-sh;
    "qlat-jhub-fhs${pkgs.qlat-name}" = pkgs.qlat-jhub-fhs;
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
    "pkgs${pkgs.qlat-name}" = pkgs;
  };

  mk-qlat-pkgs = mk-qlat-pkgs-gen {};
  mk-qlat-pkgs-cuda = mk-qlat-pkgs-gen { cudaSupport = false; }; # cudaSupport does not compile yet.

  many-qlat-pkgs = {}
  // mk-qlat-pkgs []
  ;

  many-qlat-pkgs-core = many-qlat-pkgs
  // mk-qlat-pkgs [ overlay-std overlay-clang ]
  ;

  many-qlat-pkgs-core-w-cuda = many-qlat-pkgs-core
  // mk-qlat-pkgs-cuda [ overlay-cuda ]
  ;

  many-qlat-pkgs-core-w-cuda-pypi = {}
  // mk-qlat-pkgs [ overlay-pypi ]
  // mk-qlat-pkgs-cuda [ overlay-cuda overlay-pypi ]
  // mk-qlat-pkgs [ overlay-std overlay-clang overlay-pypi ]
  ;

  many-qlat-pkgs-all = many-qlat-pkgs-core-w-cuda
  // many-qlat-pkgs-core-w-cuda-pypi
  // mk-qlat-pkgs [ overlay-clang ]
  // mk-qlat-pkgs [ overlay-clang overlay-pypi ]
  // mk-qlat-pkgs [ overlay-std ]
  // mk-qlat-pkgs [ overlay-std overlay-pypi ]
  // mk-qlat-pkgs-cuda [ overlay-std overlay-cuda ]
  // mk-qlat-pkgs-cuda [ overlay-std overlay-cuda overlay-pypi ]
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
  inherit many-qlat-pkgs;
  inherit many-qlat-pkgs-core;
  inherit many-qlat-pkgs-core-w-cuda;
  inherit many-qlat-pkgs-core-w-cuda-pypi;
  inherit many-qlat-pkgs-all;
  #
  inherit (many-qlat-pkgs-all) qlat_utils;
  inherit (many-qlat-pkgs-all) qlat_utils-std;
  inherit (many-qlat-pkgs-all) qlat_utils-std-clang;
  inherit (many-qlat-pkgs-all) qlat_utils-cuda;
  inherit (many-qlat-pkgs-all) qlat_utils-std-cuda;
  inherit (many-qlat-pkgs-all) qlat_utils-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat;
  inherit (many-qlat-pkgs-all) qlat-pypi;
  inherit (many-qlat-pkgs-all) qlat-std;
  inherit (many-qlat-pkgs-all) qlat-std-clang;
  inherit (many-qlat-pkgs-all) qlat-cuda;
  inherit (many-qlat-pkgs-all) qlat-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat_grid;
  inherit (many-qlat-pkgs-all) qlat_grid-pypi;
  inherit (many-qlat-pkgs-all) qlat_grid-cuda;
  inherit (many-qlat-pkgs-all) qlat_grid-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat_cps;
  inherit (many-qlat-pkgs-all) qlat_cps-pypi;
  inherit (many-qlat-pkgs-all) qlat_cps-cuda;
  inherit (many-qlat-pkgs-all) qlat_cps-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-tests;
  inherit (many-qlat-pkgs-all) qlat-tests-pypi;
  inherit (many-qlat-pkgs-all) qlat-tests-std;
  inherit (many-qlat-pkgs-all) qlat-tests-std-clang;
  inherit (many-qlat-pkgs-all) qlat-tests-cuda;
  inherit (many-qlat-pkgs-all) qlat-tests-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-tests-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-env;
  inherit (many-qlat-pkgs-all) qlat-env-pypi;
  inherit (many-qlat-pkgs-all) qlat-env-std;
  inherit (many-qlat-pkgs-all) qlat-env-std-clang;
  inherit (many-qlat-pkgs-all) qlat-env-cuda;
  inherit (many-qlat-pkgs-all) qlat-env-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-env-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-jhub-env;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-pypi;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-std;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-std-clang;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-env-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-sh;
  inherit (many-qlat-pkgs-all) qlat-sh-pypi;
  inherit (many-qlat-pkgs-all) qlat-sh-std;
  inherit (many-qlat-pkgs-all) qlat-sh-std-clang;
  inherit (many-qlat-pkgs-all) qlat-sh-cuda;
  inherit (many-qlat-pkgs-all) qlat-sh-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-sh-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-jhub-sh;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-pypi;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-std;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-std-clang;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-sh-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-fhs;
  inherit (many-qlat-pkgs-all) qlat-fhs-pypi;
  inherit (many-qlat-pkgs-all) qlat-fhs-std;
  inherit (many-qlat-pkgs-all) qlat-fhs-std-clang;
  inherit (many-qlat-pkgs-all) qlat-fhs-cuda;
  inherit (many-qlat-pkgs-all) qlat-fhs-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-fhs-clang; # not working
  #
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-pypi;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-std;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-std-clang;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-std-cuda;
  inherit (many-qlat-pkgs-all) qlat-jhub-fhs-clang; # not working
  #
  inherit (many-qlat-pkgs-all) pkgs;
  inherit (many-qlat-pkgs-all) pkgs-pypi;
  inherit (many-qlat-pkgs-all) pkgs-std;
  inherit (many-qlat-pkgs-all) pkgs-std-clang;
  inherit (many-qlat-pkgs-all) pkgs-cuda;
  inherit (many-qlat-pkgs-all) pkgs-std-cuda;
  inherit (many-qlat-pkgs-all) pkgs-clang;
  #
}
