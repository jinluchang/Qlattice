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

  options-default = {
    qlat-name = "";
    ngpu = "1";
    nvcc-arch = "sm_86";
    use-grid-gpt = true;
    use-cuda = false;
    use-cudasupport = false;
    use-cubaquad = true;
    use-clang = false;
    use-ucx = true;
    use-pypi = false;
  } // {
    ngpu = ngpu;
    nvcc-arch = nvcc-arch;
  };

  mk-overlay = options: final: prev: let
    opts = options-default // options;
    #
    pkgs = prev;
    lib = pkgs.lib;
    #
    call-pkg = final.callPackage;
    py-call-pkg = final.python3.pkgs.callPackage;
    #
    qlat-name = opts.qlat-name
    + lib.optionalString (! opts.use-grid-gpt) "-std"
    + lib.optionalString opts.use-cuda "-cuda"
    + lib.optionalString opts.use-cudasupport (assert opts.use-cuda; "support")
    + lib.optionalString (! opts.use-cubaquad) (assert (! opts.use-grid-gpt); "-cubaquadless")
    + lib.optionalString opts.use-clang "-clang"
    + lib.optionalString (! opts.use-ucx) "-ucxless"
    + lib.optionalString opts.use-pypi "-pypi"
    ;
    #
    qlat-stdenv = if ! opts.use-clang
    then pkgs.stdenv
    else pkgs.clangStdenv;
    #
    openmp = if ! opts.use-clang
    then null
    else pkgs.llvmPackages.openmp;
    #
    qlat-eigen = if opts.use-grid-gpt
    then grid-lehner
    else pkgs.eigen;
    #
    qlat-cc = if ! opts.use-clang
    then [ pkgs.gcc ]
    else [ pkgs.clang openmp ];
    #
    ucx = pkgs.ucx.override {
      enableCuda = opts.use-cuda;
    };
    ucx-dev = pkgs.buildEnv {
      name = "ucx-dev";
      paths = [ ucx ];
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    mpi = (pkgs.mpi.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags ++ (
        let
          cudaSupport = opts.use-cuda;
          cudaPackages = pkgs.cudaPackages;
        in [
          (lib.withFeatureAs opts.use-ucx "with-ucx" "${lib.getDev ucx-dev}")
          (lib.withFeatureAs cudaSupport "cuda-libdir" "${cudaPackages.cuda_cudart.stubs}/lib")
        ]);
    })).override { cudaSupport = opts.use-cuda; };
    python3 = pkgs.python3.override {
      packageOverrides = final: prev: {
        mpi4py = prev.mpi4py.overridePythonAttrs (prev: {
          doCheck = true;
          nativeBuildInputs = (prev.nativeBuildInputs or [])
          ++ lib.optionals opts.use-cuda [
            qlat-nixgl
            pkgs.which
          ];
          preInstallCheck = lib.optionalString opts.use-cuda ''
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
    grid-lehner-c-lime = qio;
    #
    qlat-nixgl = if opts.use-cuda then nixgl.auto.nixGLDefault else null;
    #
    cubaquad = if opts.use-cubaquad
    then call-pkg ./cubaquad.nix { stdenv = qlat-stdenv; }
    else null;
    #
    c-lime = call-pkg ./c-lime.nix { stdenv = qlat-stdenv; };
    qmp = call-pkg ./qmp.nix { stdenv = qlat-stdenv; };
    qio = call-pkg ./qio.nix { stdenv = qlat-stdenv; };
    cps = call-pkg ./cps.nix { stdenv = qlat-stdenv; };
    #
    qlat_utils = py-call-pkg ./qlat_utils.nix {
      stdenv = qlat-stdenv;
      eigen = qlat-eigen;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
    };
    qlat = py-call-pkg ./qlat.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
    };
    qlat_grid = py-call-pkg ./qlat_grid.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
    };
    qlat_cps = py-call-pkg ./qlat_cps.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
    };
    grid-lehner = call-pkg ./grid-lehner.nix {
      stdenv = qlat-stdenv;
      c-lime = grid-lehner-c-lime;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
    };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
    };
    #
    qlat-examples-cpp = py-call-pkg ./qlat-examples-cpp.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      ngpu = opts.ngpu;
    };
    qlat-examples-cpp-grid = py-call-pkg ./qlat-examples-cpp-grid.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      ngpu = opts.ngpu;
    };
    qlat-examples-py = py-call-pkg ./qlat-examples-py.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      ngpu = opts.ngpu;
    };
    qlat-examples-py-gpt = py-call-pkg ./qlat-examples-py-gpt.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      ngpu = opts.ngpu;
    };
    qlat-examples-py-cps = py-call-pkg ./qlat-examples-py-cps.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      ngpu = opts.ngpu;
    };
    qlat_docs = py-call-pkg ./qlat_docs.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
    };
    qlat_pypipkgs = py-call-pkg ./qlat_pypipkgs.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
    };
    #
    qlat-dep-pkgs = with pkgs; ([
      git pkg-config zlib gsl fftw fftwFloat hdf5-cpp openssl gmp mpfr
    ] ++ (if opts.use-cuda then [ qlat-nixgl ] else [])
    );
    #
    qlat-dep-pkgs-extra = with pkgs;
    if opts.use-grid-gpt then [
      mpi cubaquad qlat-eigen cps qmp qio grid-lehner
    ] else [
      mpi cubaquad qlat-eigen
    ];
    qlat-py-pkgs = with pkgs;
    if opts.use-grid-gpt then [
      qlat_utils
      qlat
      qlat_grid
      qlat_cps
      gpt-lehner
      qlat_docs
      qlat_pypipkgs
    ] else [
      qlat_utils
      qlat
    ];
    qlat-tests-pkgs = with pkgs;
    if opts.use-grid-gpt then [
      qlat-examples-cpp
      qlat-examples-cpp-grid
      qlat-examples-py
      qlat-examples-py-gpt
      qlat-examples-py-cps
    ] else [
      qlat-examples-cpp
      qlat-examples-py
    ];
    #
    qlat-py = python3.withPackages (ps: [
      ps.build
      ps.wheel
    ] ++ qlat-py-pkgs);
    qlat-pkgs = with pkgs; [
      qlat-py
    ] ++ qlat-dep-pkgs ++ qlat-dep-pkgs-extra;
    qlat-tests = pkgs.buildEnv {
      name = "qlat-tests${qlat-name}";
      paths = qlat-tests-pkgs;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-env = pkgs.buildEnv {
      name = "qlat-env${qlat-name}";
      paths = qlat-pkgs;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-sh = pkgs.mkShell rec {
      name = "qlat-sh${qlat-name}";
      packages = [ qlat-env ];
      inputsFrom = packages;
    };
    qlat-fhs = pkgs.buildFHSEnv {
      name = "qlat-fhs${qlat-name}";
      targetPkgs = pkgs: [ qlat-env ];
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
      profile=''
        PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
    #
    qlat-jhub-py = python3.withPackages (ps: with ps; [
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
      pytools
      lz4
      torchvision
      torchaudio
      xformers
      jupyterlab
      jupyterhub
      jupyterhub-systemdspawner
    ]
    ++ qlat-py-pkgs
    ++ lib.optionals opts.use-cuda [
      pycuda
    ]);
    qlat-jhub-env = pkgs.buildEnv {
      name = "qlat-jhub-env${qlat-name}";
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
      ] ++ qlat-cc ++ qlat-dep-pkgs ++ qlat-dep-pkgs-extra;
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
    };
    qlat-jhub-sh = pkgs.mkShell rec {
      name = "qlat-jhub-sh${qlat-name}";
      packages = [ qlat-jhub-env ];
      inputsFrom = packages;
    };
    qlat-jhub-fhs = pkgs.buildFHSEnv {
      name = "qlat-jhub-fhs${qlat-name}";
      targetPkgs = pkgs: [ qlat-jhub-env ];
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "out" "doc" ];
      profile=''
        PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
    #
  in {
    inherit qlat-name;
    inherit python3 mpi openmp ucx;
    inherit c-lime qmp qio cps cubaquad grid-lehner gpt-lehner;
    inherit (opts) use-pypi;
    inherit qlat_utils qlat qlat_grid qlat_cps;
    inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt qlat-examples-py-cps;
    inherit qlat_docs qlat_pypipkgs;
    inherit qlat-py qlat-pkgs qlat-tests qlat-env qlat-sh qlat-fhs;
    inherit qlat-jhub-py qlat-jhub-env qlat-jhub-sh qlat-jhub-fhs;
  };

  mk-qlat-pkgs = options: let
    opts = options-default // options;
    pkgs = nixpkgs {
      config = {
        allowUnfree = opts.use-cuda;
        cudaSupport = opts.use-cudasupport;
      };
      overlays = [
        (mk-overlay options)
      ];
    };
  in {
    "qlat_utils${pkgs.qlat-name}" = pkgs.qlat_utils;
    "qlat${pkgs.qlat-name}" = pkgs.qlat;
    "qlat_grid${pkgs.qlat-name}" = pkgs.qlat_grid;
    "qlat_cps${pkgs.qlat-name}" = pkgs.qlat_cps;
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
    "qlat-py${pkgs.qlat-name}" = pkgs.qlat-py;
    "qlat-tests${pkgs.qlat-name}" = pkgs.qlat-tests;
    "qlat-env${pkgs.qlat-name}" = pkgs.qlat-env;
    "qlat-sh${pkgs.qlat-name}" = pkgs.qlat-sh;
    "qlat-fhs${pkgs.qlat-name}" = pkgs.qlat-fhs;
    "qlat-jhub-py${pkgs.qlat-name}" = pkgs.qlat-jhub-py;
    "qlat-jhub-env${pkgs.qlat-name}" = pkgs.qlat-jhub-env;
    "qlat-jhub-sh${pkgs.qlat-name}" = pkgs.qlat-jhub-sh;
    "qlat-jhub-fhs${pkgs.qlat-name}" = pkgs.qlat-jhub-fhs;
    "pkgs${pkgs.qlat-name}" = pkgs;
  };

  many-qlat-pkgs = {}
  // mk-qlat-pkgs {}
  ;
  many-qlat-pkgs-core = many-qlat-pkgs
  // mk-qlat-pkgs { use-grid-gpt = false; use-clang = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cubaquad = false; }
  // mk-qlat-pkgs { use-ucx = false; }
  ;
  many-qlat-pkgs-core-w-cuda = many-qlat-pkgs-core
  // mk-qlat-pkgs { use-cuda = true; }
  // mk-qlat-pkgs { use-cuda = true; use-ucx = false; }
  ;
  many-qlat-pkgs-core-w-cuda-pypi = {}
  // mk-qlat-pkgs { use-pypi = true; }
  // mk-qlat-pkgs { use-ucx = false; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-clang = true; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cubaquad = false; use-pypi = true; }
  // mk-qlat-pkgs { use-cuda = true; use-pypi = true; }
  // mk-qlat-pkgs { use-cuda = true; use-ucx = false; use-pypi = true; }
  ;
  many-qlat-pkgs-all = many-qlat-pkgs-core-w-cuda
  // many-qlat-pkgs-core-w-cuda-pypi
  // mk-qlat-pkgs { use-grid-gpt = false; use-ucx = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; use-ucx = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; use-pypi = true; }
  // mk-qlat-pkgs { use-clang = true; }
  // mk-qlat-pkgs { use-clang = true; use-pypi = true; }
  // mk-qlat-pkgs { use-cuda = true; use-cudasupport = true; }
  // mk-qlat-pkgs { use-cuda = true; use-cudasupport = true; use-pypi = true; }
  ;

in {
  #
  inherit mk-qlat-pkgs;
  inherit mk-overlay;
  #
  inherit many-qlat-pkgs;
  inherit many-qlat-pkgs-core;
  inherit many-qlat-pkgs-core-w-cuda;
  inherit many-qlat-pkgs-core-w-cuda-pypi;
  inherit many-qlat-pkgs-all;
  #
  inherit (many-qlat-pkgs-all) qlat_utils;
  inherit (many-qlat-pkgs-all) qlat_utils-pypi;
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
  inherit (many-qlat-pkgs-all) pkgs-cudasupport;
  #
}
