{
  nixpkgs ? import ./nixpkgs.nix,
  ngpu ? "2", # adjust with actual number of GPUs
  nvcc-arch ? "sm_86", # adjust with actual arch of Nvidia GPU
  cudaCapability ? "8.6",
  cudaForwardCompat ? false,
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
    cudaCapabilities = [ "8.6" ];
    cudaForwardCompat = false;
  } // {
    ngpu = ngpu;
    nvcc-arch = nvcc-arch;
    cudaCapabilities = [ cudaCapability ];
    cudaForwardCompat = cudaForwardCompat;
  };

  mk-qlat-name = options:
  let
    opts = options;
    lib = (nixpkgs {}).lib;
  in opts.qlat-name
  + lib.optionalString (! opts.use-grid-gpt) "-std"
  + lib.optionalString opts.use-cuda "-cuda"
  + lib.optionalString opts.use-cudasupport (assert opts.use-cuda; "support")
  + lib.optionalString (! opts.use-cubaquad) (assert (! opts.use-grid-gpt); "-cubaquadless")
  + lib.optionalString opts.use-clang "-clang"
  + lib.optionalString (! opts.use-ucx) "-ucxless"
  + lib.optionalString opts.use-pypi "-pypi"
  ;

  mk-overlay = options: final: prev: let
    opts = options-default // options;
    #
    pkgs = prev;
    lib = pkgs.lib;
    #
    call-pkg = final.callPackage;
    py-call-pkg = final.python3.pkgs.callPackage;
    #
    qlat-name = mk-qlat-name opts;
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
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
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
    qlat-py = python3.withPackages (ps: qlat-py-pkgs);
    qlat-pkgs = with pkgs; [
      qlat-py
    ] ++ qlat-dep-pkgs ++ qlat-dep-pkgs-extra;
    qlat-tests = pkgs.buildEnv {
      name = "qlat-tests${qlat-name}";
      paths = qlat-tests-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
    };
    qlat-env = pkgs.buildEnv {
      name = "qlat-env${qlat-name}";
      paths = qlat-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
    };
    qlat-sh = pkgs.mkShell rec {
      name = "qlat-sh${qlat-name}";
      packages = [ qlat-env ];
      inputsFrom = packages;
      buildInputs = [ pkgs.pkg-config ];
    };
    qlat-fhs = pkgs.buildFHSEnv {
      name = "qlat-fhs${qlat-name}";
      targetPkgs = pkgs: [
        qlat-env
        pkgs.pkg-config
      ] ++ qlat-cc;
      multiPkgs = pkgs: [
        qlat-env
        pkgs.pkg-config
      ] ++ qlat-cc;
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "static" "man" "doc" "info" ];
      profile=''
        # PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
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
      ipywidgets
      accelerate
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
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
    };
    qlat-jhub-sh = pkgs.mkShell rec {
      name = "qlat-jhub-sh${qlat-name}";
      packages = [ qlat-jhub-env ];
      inputsFrom = packages;
      buildInputs = [ pkgs.pkg-config ];
    };
    qlat-jhub-fhs = pkgs.buildFHSEnv {
      name = "qlat-jhub-fhs${qlat-name}";
      targetPkgs = pkgs: [
        qlat-jhub-env
        pkgs.pkg-config
      ] ++ qlat-cc;
      multiPkgs = pkgs: [
        qlat-jhub-env
        pkgs.pkg-config
      ] ++ qlat-cc;
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "static" "man" "doc" "info" ];
      profile=''
        # PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
    #
    qlat-jhub-tests = {
      inherit qlat-tests qlat-jhub-env;
    };
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
    inherit qlat-jhub-tests;
  };

  mk-qlat-pkgs = options: let
    opts = options-default // options;
    pkgs = nixpkgs {
      config = {
        allowUnfree = opts.use-cuda;
        cudaSupport = opts.use-cudasupport;
        ${if opts.use-cudasupport then "cudaCapabilities" else null} = opts.cudaCapabilities;
        ${if opts.use-cudasupport then "cudaForwardCompat" else null} = opts.cudaForwardCompat;
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
    "qlat-jhub-tests${pkgs.qlat-name}" = pkgs.qlat-jhub-tests;
    "pkgs${pkgs.qlat-name}" = pkgs;
  };

  many-qlat-pkgs = mk-qlat-pkgs {}
  ;
  many-qlat-pkgs-pypi = mk-qlat-pkgs { use-pypi = true; }
  ;
  many-qlat-pkgs-more = many-qlat-pkgs
  // mk-qlat-pkgs { use-ucx = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-clang = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cubaquad = false; }
  ;
  many-qlat-pkgs-more-pypi = many-qlat-pkgs-pypi
  // mk-qlat-pkgs { use-ucx = false; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-clang = true; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cubaquad = false; use-pypi = true; }
  ;
  many-qlat-pkgs-more-w-cuda = many-qlat-pkgs-more
  // mk-qlat-pkgs { use-cuda = true; }
  // mk-qlat-pkgs { use-cuda = true; use-ucx = false; }
  ;
  many-qlat-pkgs-more-w-cuda-pypi = many-qlat-pkgs-more-pypi
  // mk-qlat-pkgs { use-cuda = true; use-pypi = true; }
  // mk-qlat-pkgs { use-cuda = true; use-ucx = false; use-pypi = true; }
  ;
  many-qlat-pkgs-extra = many-qlat-pkgs-more-w-cuda
  // mk-qlat-pkgs { use-grid-gpt = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-ucx = false; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; use-ucx = false; }
  // mk-qlat-pkgs { use-cuda = true; use-cudasupport = true; }
  ;
  many-qlat-pkgs-extra-pypi = many-qlat-pkgs-more-w-cuda-pypi
  // mk-qlat-pkgs { use-grid-gpt = false; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-ucx = false; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; use-pypi = true; }
  // mk-qlat-pkgs { use-grid-gpt = false; use-cuda = true; use-ucx = false; use-pypi = true; }
  // mk-qlat-pkgs { use-cuda = true; use-cudasupport = true; use-pypi = true; }
  ;
  many-qlat-pkgs-all = many-qlat-pkgs-extra
  // many-qlat-pkgs-extra-pypi
  // mk-qlat-pkgs { use-clang = true; }
  // mk-qlat-pkgs { use-clang = true; use-pypi = true; }
  ;

in many-qlat-pkgs-all // {
  #
  inherit mk-qlat-pkgs;
  inherit mk-overlay;
  #
  inherit many-qlat-pkgs;
  inherit many-qlat-pkgs-pypi;
  inherit many-qlat-pkgs-more;
  inherit many-qlat-pkgs-more-pypi;
  inherit many-qlat-pkgs-more-w-cuda;
  inherit many-qlat-pkgs-more-w-cuda-pypi;
  inherit many-qlat-pkgs-extra;
  inherit many-qlat-pkgs-extra-pypi;
  inherit many-qlat-pkgs-all;
}
