{
  # Path or fetchTarball result to nixpkgs. Auto-detected from `version` if null.
  # E.g. (fetchTarball "https://channels.nixos.org/nixos-26.05/nixexprs.tar.xz")
  nixpkgs ? null,
  # NixOS release version string. Used to fetch matching nixpkgs when `nixpkgs` is null.
  # E.g. "26.05"
  version ? null,
  # Number of NVIDIA GPUs as a string. Auto-detected from /dev/nvidia* if null.
  # E.g. "2"
  ngpu ? null,
  # CUDA compute capability string (major.minor). Auto-detected via nvidia-smi if null.
  # E.g. "8.6"
  cudaCapability ? null,
  # Enable CUDA forward compatibility (bool). Defaults to false if null.
  cudaForwardCompat ? null,
  # Use gitee mirrors instead of github (bool). Defaults to false if null.
  use-gitee ? null,
}:

let

  opts = import ./options.nix {
    inherit nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;
  };

in let

  inherit (opts) nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;
  inherit (opts) options-default mk-options mk-qlat-name;
  inherit (opts) import-nixpkgs import-nixgl;
  inherit (opts) nixpkgs-release;
  inherit (opts) version-pypi;
  inherit (opts) options-list qlat-name-list qlat-name-list-file-from-str qlat-name-list-file;

  mk-overlay = options: final-g: prev-g: let
    opts = mk-options options;
    #
    pkgs = prev-g;
    #
    lib = pkgs.lib;
    #
    is-linux = (lib.lists.elem builtins.currentSystem lib.platforms.linux);
    #
    call-pkg = pkgs.callPackage;
    py-call-pkg = python3.pkgs.callPackage;
    #
    qlat-name = mk-qlat-name opts;
    #
    qlat-stdenv = if opts.use-clang
    then pkgs.clangStdenv
    else if opts.use-cuda
    then pkgs.cudaPackages.backendStdenv
    else pkgs.gccStdenv;
    #
    openmp = if qlat-stdenv.cc.isClang
    then pkgs.llvmPackages.openmp
    else null;
    #
    qlat-eigen = if opts.use-grid-gpt
    then grid-lehner
    else pkgs.eigen;
    #
    qlat-cc = {
      cc = qlat-stdenv.cc;
      # cc-c = qlat-stdenv.cc.libc;
      # cc-cc = qlat-stdenv.cc.cc;
      # cc-cc-c = qlat-stdenv.cc.cc.lib;
    }
    //
    { inherit (pkgs) pkg-config; }
    //
    (if opts.use-clang then { inherit openmp; } else {})
    //
    (if opts.use-cuda-software then ({
      cuda_cudart = lib.meta.hiPrio pkgs.cudaPackages.cuda_cudart;
    }
    // {
      inherit (pkgs.cudaPackages)
      cuda_nvcc
      cuda_cccl
      cuda_profiler_api
      libcufft
      libcublas
      libnpp
      cudnn
      nccl
      cudatoolkit
      fabricmanager
      ;
    }) else {}
    );
    #
    nvidia_x11_bin = if opts.ngpu == "0" || ! opts.use-cuda-software
    then null
    else pkgs.linuxPackages.nvidia_x11.bin;
    #
    ucx-mt = (pkgs.ucx.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags ++ [
        "--enable-mt"
      ];
    })).override {
      enableCuda = opts.use-cuda-software;
    };
    ucx-mt-dev = pkgs.buildEnv {
      name = "ucx-mt-dev";
      paths = [ ucx-mt ];
      extraOutputsToInstall = [ "out" "bin" "dev" "lib" "static" "man" "doc" "info" ];
    };
    mpi = if ! opts.use-ucx && ! opts.use-cuda-software
    then pkgs.mpi
    else if ! opts.use-cuda-software
    then pkgs.mpi.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags
      ++ lib.optionals opts.use-ucx [ "--with-ucx=${lib.getDev ucx-mt-dev}" ];
    })
    else (pkgs.mpi.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags
      ++ (let
        cudaPackages = pkgs.cudaPackages;
        cudaSupport = opts.use-cuda-software;
        use-ucx = opts.use-ucx;
      in [
        (lib.withFeatureAs use-ucx "with-ucx" "${lib.getDev ucx-mt-dev}")
        (lib.withFeatureAs cudaSupport "cuda-libdir" "${lib.getLib cudaPackages.cuda_cudart}/lib/stubs")
      ]
      );
      env.NIX_CFLAGS_COMPILE = lib.concatStringsSep " " [
        "-Wno-error=int-conversion"
        "-Wno-error=incompatible-pointer-types"
        "-Wno-error=implicit-function-declaration"
      ];
    })).override {
      cudaSupport = opts.use-cuda-software;
    };
    python3 = pkgs.python3.override {
      packageOverrides = final: prev: rec {
        mpi4py = if ! opts.use-cuda-software
        then prev.mpi4py
        else prev.mpi4py.overridePythonAttrs (py-prev: {
          doCheck = false;
          nativeBuildInputs = (py-prev.nativeBuildInputs or []) ++ [
            qlat-nixgl
            pkgs.which
          ];
          preInstallCheck = ''
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
        jax = if ! opts.use-cudasupport
        then prev.jax
        else prev.jax.overridePythonAttrs (py-prev: {
          doCheck = false;
        });
        keras = if ! opts.use-cudasupport
        then prev.keras
        else prev.keras.overridePythonAttrs (py-prev: {
          doCheck = false;
        });
        accelerate = if ! opts.use-cudasupport
        then prev.accelerate
        else prev.accelerate.overridePythonAttrs (py-prev: {
          doCheck = false;
        });
        gvar = pkgs.python3.pkgs.callPackage ./gvar.nix {};
        vegas = pkgs.python3.pkgs.callPackage ./vegas.nix { gvar = gvar; };
        lsqfit = pkgs.python3.pkgs.callPackage ./lsqfit.nix { gvar = gvar; vegas = vegas; };
        corrfitter = pkgs.python3.pkgs.callPackage ./corrfitter.nix { lsqfit = lsqfit; gvar = gvar; };
        qcd_ml = pkgs.python3.pkgs.callPackage ./qcd_ml.nix {};
        # diffusers = pkgs.python3.pkgs.callPackage
        # (import "${n-pkgs-src}/pkgs/development/python-modules/diffusers/default.nix") {
        #   jax = jax; diffusers = diffusers;
        # };
        # ollama = pkgs.python3.pkgs.callPackage
        # (import "${n-pkgs-src}/pkgs/development/python-modules/ollama/default.nix") {
        # };
      };
    };
    #
    grid-lehner-c-lime = if opts.use-cps then qio else c-lime;
    #
    nixgl = (import-nixgl { pkgs = pkgs; }).auto.nixGLDefault;
    #
    qlat-nixgl = if opts.use-cuda-software then nixgl else null;
    #
    cuba-quad = if opts.use-cubaquad
    then call-pkg ./cuba-quad.nix { stdenv = qlat-stdenv; }
    else null;
    #
    c-lime = call-pkg ./c-lime.nix { stdenv = qlat-stdenv; };
    qmp = call-pkg ./qmp.nix { stdenv = qlat-stdenv; };
    qio = call-pkg ./qio.nix { stdenv = qlat-stdenv; };
    cps = call-pkg ./cps.nix { stdenv = qlat-stdenv; };
    grid-lehner = call-pkg ./grid-lehner.nix {
      stdenv = qlat-stdenv;
      c-lime = grid-lehner-c-lime;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      cpuinfo-sys = opts.cpuinfo-sys;
    };
    #
    qlat-version = if opts.use-pypi != null
    then version-pypi
    else builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
    #
    version-pypi-git = let
      parts = builtins.split "\\." version-pypi;
      major = builtins.elemAt parts 0;
      minor = builtins.elemAt parts 2;
    in "${major}.${lib.fixedWidthNumber 2 (lib.toInt minor)}";
    #
    git-src = builtins.fetchGit {
      url = if use-gitee then "https://gitee.com/jinluchang/Qlattice" else "https://github.com/jinluchang/Qlattice";
      ref = "refs/tags/v${version-pypi-git}";
    };
    #
    pypi-qlat-utils = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat_utils/qlat_utils-${version-pypi}.tar.gz";
    pypi-qlat = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat/qlat-${version-pypi}.tar.gz";
    pypi-qlat-grid = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat_grid/qlat_grid-${version-pypi}.tar.gz";
    pypi-qlat-cps = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat_cps/qlat_cps-${version-pypi}.tar.gz";
    #
    qlat-src = if opts.use-pypi != null
    then {
      qlat-utils = pypi-qlat-utils;
      qlat = pypi-qlat;
      qlat-grid = pypi-qlat-grid;
      qlat-cps = pypi-qlat-cps;
      examples-cpp = "${git-src}/examples-cpp";
      examples-cpp-grid = "${git-src}/examples-cpp-grid";
      examples-py = "${git-src}/examples-py";
      examples-py-gpt = "${git-src}/examples-py-gpt";
      examples-py-cps = "${git-src}/examples-py-cps";
      docs = "${git-src}/docs";
      qcore = "${git-src}/qcore";
    }
    else {
      qlat-utils = ../qlat-utils;
      qlat = ../qlat;
      qlat-grid = ../qlat-grid;
      qlat-cps = ../qlat-cps;
      examples-cpp = ../examples-cpp;
      examples-cpp-grid = ../examples-cpp-grid;
      examples-py = ../examples-py;
      examples-py-gpt = ../examples-py-gpt;
      examples-py-cps = ../examples-py-cps;
      docs = ../docs;
      qcore = ../qcore;
    };
    #
    qlat_utils = py-call-pkg ./qlat_utils.nix {
      stdenv = qlat-stdenv;
      eigen = qlat-eigen;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      nixgl = qlat-nixgl;
      version = qlat-version;
    };
    qlat = py-call-pkg ./qlat.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      version = qlat-version;
    };
    qlat_grid = py-call-pkg ./qlat_grid.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      version = qlat-version;
    };
    qlat_cps = py-call-pkg ./qlat_cps.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      version = qlat-version;
    };
    gpt-lehner = py-call-pkg ./gpt-lehner.nix {
      stdenv = qlat-stdenv;
      cudaSupport = opts.use-cuda;
    };
    #
    qlat-examples-cpp = py-call-pkg ./qlat-examples-cpp.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda-software;
      nvcc-arch = opts.nvcc-arch;
      ngpu = opts.ngpu;
      version = qlat-version;
    };
    qlat-examples-cpp-grid = py-call-pkg ./qlat-examples-cpp-grid.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda-software;
      nvcc-arch = opts.nvcc-arch;
      ngpu = opts.ngpu;
      version = qlat-version;
    };
    qlat-examples-py = py-call-pkg ./qlat-examples-py.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda-software;
      nvcc-arch = opts.nvcc-arch;
      ngpu = opts.ngpu;
      version = qlat-version;
    };
    qlat-examples-py-gpt = py-call-pkg ./qlat-examples-py-gpt.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda-software;
      nvcc-arch = opts.nvcc-arch;
      ngpu = opts.ngpu;
      version = qlat-version;
    };
    qlat-examples-py-cps = py-call-pkg ./qlat-examples-py-cps.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda-software;
      nvcc-arch = opts.nvcc-arch;
      ngpu = opts.ngpu;
      version = qlat-version;
    };
    qlat_docs = py-call-pkg ./qlat_docs.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      version = qlat-version;
    };
    qlat_pypipkgs = py-call-pkg ./qlat_pypipkgs.nix {
      stdenv = qlat-stdenv;
      qlat-src = qlat-src;
      cudaSupport = opts.use-cuda;
      nvcc-arch = opts.nvcc-arch;
      version = qlat-version;
    };
    #
    qlat-dep-pkgs = {
      inherit (pkgs) git zlib gsl fftw fftwFloat hdf5-cpp openssl gmp mpfr;
    } // (if opts.use-cuda-software then { inherit qlat-nixgl; } else {});
    #
    qlat-dep-pkgs-extra = if ! opts.use-grid-gpt && ! opts.use-cps then {
      inherit mpi cuba-quad qlat-eigen;
    } else if opts.use-grid-gpt && ! opts.use-cps then {
      inherit mpi cuba-quad qlat-eigen grid-lehner c-lime;
    } else if ! opts.use-grid-gpt && opts.use-cps then {
      inherit mpi cuba-quad qlat-eigen cps qmp qio;
    } else {
      inherit mpi cuba-quad qlat-eigen cps qmp qio grid-lehner;
    };
    qlat-py-pkgs = if ! opts.use-grid-gpt && ! opts.use-cps then {
      inherit qlat_utils qlat;
    } else if opts.use-grid-gpt && ! opts.use-cps then {
      inherit qlat_utils qlat qlat_grid gpt-lehner;
    } else if ! opts.use-grid-gpt && opts.use-cps then {
      inherit qlat_utils qlat qlat_cps;
    } else {
      inherit qlat_utils qlat qlat_grid qlat_cps gpt-lehner qlat_docs qlat_pypipkgs;
    };
    qlat-tests-pkgs = if ! opts.use-grid-gpt && ! opts.use-cps then {
      inherit qlat-examples-cpp qlat-examples-py;
    } else if opts.use-grid-gpt && ! opts.use-cps then {
      inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt;
    } else if ! opts.use-grid-gpt && opts.use-cps then {
      inherit qlat-examples-cpp qlat-examples-py qlat-examples-py-cps;
    } else {
      inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt qlat-examples-py-cps;
    };
    #
    qlat-py = python3.withPackages (ps: builtins.attrValues qlat-py-pkgs);
    qlat-pkgs = {
      qlat-py = lib.meta.hiPrio qlat-py;
    }
    // {
      inherit
      qlat-nixgl
      qlat-stdenv
      mpi
      nvidia_x11_bin
      ;
    }
    // qlat-dep-pkgs
    // qlat-dep-pkgs-extra;
    qlat-tests = pkgs.buildEnv {
      name = "qlat-tests${qlat-name}";
      paths = builtins.attrValues qlat-tests-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "lib" "static" "man" "doc" "info" ];
    };
    qlat-env = pkgs.buildEnv {
      name = "qlat-env${qlat-name}";
      paths = builtins.attrValues qlat-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "lib" "static" "man" "doc" "info" ];
      # ignoreCollisions = true;
    };
    qlat-sh = pkgs.mkShell rec {
      name = "qlat-sh${qlat-name}";
      packages = [ qlat-env ];
      inputsFrom = packages;
      buildInputs = builtins.attrValues qlat-cc;
      shellHook = ''
      '';
    };
    qlat-fhs = pkgs.buildFHSEnv {
      name = "qlat-fhs${qlat-name}";
      targetPkgs = pkgs: [
        qlat-env
      ] ++ builtins.attrValues qlat-cc;
      multiPkgs = pkgs: [
        qlat-env
      ] ++ builtins.attrValues qlat-cc;
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "lib" "static" "man" "doc" "info" ];
      profile=''
        # PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
    qlat-jhub-py = python3.withPackages (ps: builtins.attrValues ({
      inherit (ps)
      ipykernel
      pip
      numpy
      scipy
      sympy
      jax
      jaxlib
      gvar
      vegas
      qcd_ml
      lsqfit
      corrfitter
      meson
      ninja
      mpi4py
      openai
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
      sentencepiece
      sphinx
      linkify-it-py
      myst-parser
      pytest
      pytools
      lz4
      soundfile
      pypdf
      reportlab
      pdfminer-six
      awkward
      uproot
      iminuit
      hepunits
      networkx
      joblib
      dask
      tables
      zarr
      tqdm
      pyarrow
      bokeh
      ipython
      mypy
      black
      memory-profiler
      tabulate
      pyyaml
      click
      tiktoken
      ;
    }
    // (if is-linux then {
      inherit (ps)
      xarray
      numba
      torch
      flax
      keras
      matplotlib
      ipywidgets
      scikit-learn
      accelerate
      jupyter-server-mathjax
      jupyterlab
      jupyterhub
      jupyterhub-systemdspawner
      xformers
      rich
      transformers
      plotly
      seaborn
      pandas
      torchvision
      torchaudio
      tensorboard
      optuna
      datasets
      tokenizers
      safetensors
      ;
    } else {}
    )
    // (if opts.use-cuda-software then {
      inherit (ps)
      pycuda
      bitsandbytes
      diffusers
      ;
    } else {}
    )
    // qlat-py-pkgs
    ));
    qlat-lib-set = ({
      qlat-jhub-py = lib.meta.hiPrio qlat-jhub-py;
      clang-tools = lib.meta.lowPrio pkgs.clang-tools;
    }
    // {
      inherit
      qlat-nixgl
      mpi
      nvidia_x11_bin
      ;
      inherit (pkgs)
      bashInteractive
      bash-completion
      stdenvNoCC
      coreutils
      openssh
      findutils
      ruff
      m4
      git
      gnumake
      zlib
      hwloc
      killall
      wget
      which
      gnugrep
      rsync
      automake
      autoconf
      gsl
      fftw
      fftwFloat
      openssl
      gnuplot
      twine
      file
      zip
      unzip
      bzip2
      curl
      expat
      icu
      libsodium
      libssh
      nspr
      nss
      zstd
      pdftk
      qpdf
      poppler-utils
      ghostscript
      ;
      pipx = pkgs.pipx.overridePythonAttrs (py-prev: { doCheck = false; });
    } // (if is-linux
    then {
      util-linux = lib.meta.hiPrio pkgs.util-linux;
      inherit (pkgs)
      linux-pam
      fuse3
      pipewire
      cups
      dbus
      acl
      alsa-lib
      at-spi2-atk
      at-spi2-core
      attr
      gdk-pixbuf
      glib
      atk
      gtk3
      pango
      cairo
      libappindicator-gtk3
      libdrm
      libGL
      libglvnd
      libnotify
      libpulseaudio
      libunwind
      libusb1
      libuuid
      libxkbcommon
      mesa
      vulkan-loader
      texliveFull
      fontconfig
      freetype
      ;
    }
    else {}
    )
    // qlat-cc
    // qlat-dep-pkgs
    // qlat-dep-pkgs-extra
    );
    qlat-jhub-env = pkgs.buildEnv {
      name = "qlat-jhub-env${qlat-name}";
      paths = builtins.attrValues qlat-lib-set;
      extraOutputsToInstall = [ "out" "bin" "dev" "lib" "static" "man" "doc" "info" ];
      # ignoreCollisions = true;
      postBuild = ''
        cat >$out/bin/setenv-qlat.sh <<EOF
          # Set up environment variables for qlat
          export PATH="$out/bin\''${PATH:+:\$PATH}"
          export LD_LIBRARY_PATH="$out/lib:$out/lib64\''${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
          export LIBRARY_PATH="$out/lib:$out/lib64\''${LIBRARY_PATH:+:\$LIBRARY_PATH}"
          export CPATH="$out/include\''${CPATH:+:\$CPATH}"
          export PKG_CONFIG_PATH="$out/lib/pkgconfig\''${PKG_CONFIG_PATH:+:\$PKG_CONFIG_PATH}"
          #
          source $out/bin/cuda-mpi-qlat.sh :
          #
          "\$@"
        EOF
      '';
    };
    qlat-jhub-sh = pkgs.mkShell rec {
      name = "qlat-jhub-sh${qlat-name}";
      packages = [ qlat-jhub-env ];
      inputsFrom = packages;
      buildInputs = builtins.attrValues qlat-cc;
      shellHook = ''
      '';
    };
    qlat-jhub-fhs = pkgs.buildFHSEnv {
      name = "qlat-jhub-fhs${qlat-name}";
      targetPkgs = pkgs: [
        qlat-jhub-env
      ] ++ builtins.attrValues qlat-cc;
      multiPkgs = pkgs: [
        qlat-jhub-env
      ] ++ builtins.attrValues qlat-cc;
      runScript = "bash";
      extraOutputsToInstall = [ "bin" "dev" "lib" "static" "man" "doc" "info" ];
      profile=''
        # PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig"
      '';
    };
    #
    q-pkgs = {
      inherit qlat-py qlat-env qlat-sh qlat-fhs;
      inherit qlat-jhub-py qlat-jhub-env qlat-jhub-sh qlat-jhub-fhs;
      inherit qlat-tests;
    };
    #
  in {
    inherit qlat-name;
    inherit use-gitee;
    inherit qlat-nixgl qlat-stdenv;
    inherit python3 mpi openmp ucx-mt ucx-mt-dev;
    inherit c-lime qmp qio cps cuba-quad grid-lehner gpt-lehner;
    inherit qlat-src qlat-version;
    inherit qlat_utils qlat qlat_grid qlat_cps;
    inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt qlat-examples-py-cps;
    inherit qlat_docs qlat_pypipkgs;
    inherit qlat-py qlat-tests qlat-env qlat-sh qlat-fhs;
    inherit qlat-jhub-py qlat-jhub-env qlat-jhub-sh qlat-jhub-fhs;
    inherit qlat-pkgs q-pkgs;
    inherit qlat-lib-set;
    qlat-options = opts;
    qlat-use-pypi = opts.use-pypi;
  };

  mk-q-pkgs = options: let
    opts = mk-options options;
    qlat-name = mk-qlat-name opts;
    pkgs = import-nixpkgs {
      config = {
        allowUnfree = opts.use-cuda-software;
        cudaSupport = opts.use-cudasupport;
        ${if opts.use-cudasupport then "cudaCapabilities" else null} = opts.cudaCapabilities;
        ${if opts.use-cudasupport then "cudaForwardCompat" else null} = opts.cudaForwardCompat;
      };
      overlays = [
        (mk-overlay opts)
      ];
    };
  in {
    "qlat-pkgs${qlat-name}" = pkgs.qlat-pkgs;
    "q-pkgs${qlat-name}" = pkgs.q-pkgs;
    "pkgs${qlat-name}" = pkgs;
    #
    "qlat-env${qlat-name}" = pkgs.qlat-env;
    "qlat-tests${qlat-name}" = pkgs.qlat-tests;
  };

  q-pkgs-list = builtins.map mk-q-pkgs options-list;

  all-q-pkgs = builtins.foldl' (s: v: s // v) {} q-pkgs-list;

  all-qlat-env = builtins.listToAttrs (builtins.map (v: {
    name = "qlat-env${v}";
    value = all-q-pkgs."qlat-env${v}";
  }) qlat-name-list);

  all-qlat-tests = builtins.listToAttrs (builtins.map (v: {
    name = "qlat-tests${v}";
    value = all-q-pkgs."qlat-tests${v}";
  }) qlat-name-list);

  q-pkgs = all-q-pkgs // {
    inherit mk-q-pkgs;
    inherit mk-overlay;
    inherit options-list qlat-name-list qlat-name-list-file-from-str;
    inherit qlat-name-list-file;
    inherit all-qlat-env all-qlat-tests;
    inherit options-default;
    inherit mk-options mk-qlat-name;
    inherit nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;
  };

in builtins.deepSeq [
  nixpkgs-release
  options-default
  qlat-name-list
  qlat-name-list-file-from-str
] q-pkgs
