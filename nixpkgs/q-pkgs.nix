{
  nixpkgs ? null, # nixpkgs. E.g. (fetchTarball "https://channels.nixos.org/nixos-25.05/nixexprs.tar.xz")
  version ? null, # version of the nixpkgs. E.g. "25.05"
  ngpu ? null, # adjust with desired number of GPUs. E.g. "2"
  cudaCapability ? null, # adjust with desired cudaCapability. E.g. "8.6"
  cudaForwardCompat ? null, # adjust with desired cudaForwardCompat. E.g. false
  use-gitee ? null, # true or false (default false)
}:

let

  version-wd = if version == null then "25.05" else version;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
  nixpkgs-default = import ./nixpkgs.nix {
    version = version-wd;
    use-gitee = use-gitee-wd;
  };
  nixpkgs-wd = if nixpkgs == null then nixpkgs-default else nixpkgs;
  import-nixpkgs-wd = import nixpkgs-wd;

  version-pypi = "0.83";

  o-pkgs = import-nixpkgs-wd {
    config.allowUnfree = true;
  };

  lib = o-pkgs.lib;

  nixpkgs-release = lib.trivial.release;

  runCommandLocal = o-pkgs.runCommandLocal;

  runCommand = o-pkgs.runCommand;

  # n-pkgs-src = builtins.fetchTarball "https://channels.nixos.org/nixos-unstable-small/nixexprs.tar.xz";

  # n-pkgs = import n-pkgs-src {
  #   config.allowUnfree = true;
  # };

  nixgl-src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/nixGL" else "https://github.com/jinluchang/nixGL";
    ref = "main";
    rev = "f72374ff8c0dccef46b7b02e7cde09f91007aac6";
  };

  nixgl = (import nixgl-src {}).auto.nixGLDefault;

  is-linux = lib.lists.elem builtins.currentSystem lib.platforms.linux;

  cpuinfo-sys = builtins.readFile (runCommandLocal
  "impure-cpuinfo-file"
  {
    time = builtins.currentTime;
  }
  ''
    cat /proc/cpuinfo >$out 2>/dev/null || echo >$out
    echo "cpuinfo="
    echo "$(grep '^flags' $out 2>/dev/null | head -n 1)"
  ''
  );

  ngpu-sys = builtins.head (builtins.match
  "(.*)\n"
  (builtins.readFile (runCommandLocal
  "impure-ngpu-file"
  {
    time = builtins.currentTime;
  }
  ''
    mkdir tmp
    cd tmp
    ls /dev/nvidia{?,??} 2>/dev/null | wc -l >$out 2>/dev/null || echo "0" >$out
    echo "ngpu=$(cat $out)"
  ''
  )));

  nvidia_x11_bin = if ngpu-sys == "0"
  then null
  else o-pkgs.linuxPackages.nvidia_x11.bin;

  cudaCapability-sys = if ngpu-sys == "0"
  then null
  else builtins.head (builtins.match
  "(.*)\n"
  (builtins.readFile (runCommandLocal
  "impure-cuda-capability-file"
  {
    time = builtins.currentTime;
  }
  ''
    ${nixgl}/bin/nixGL ${nvidia_x11_bin}/bin/nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
      | head -n 1 >$out 2>/dev/null \
      || echo >$out
    echo "cudaCapability=$(cat $out)"
  ''
  )));

  cudaCapabilities-sys = if cudaCapability-sys == null
  then []
  else [ cudaCapability-sys ];

  get-nvcc-arch-from-cudaCapability = cudaCapability:
  "sm_" + builtins.replaceStrings [ "." ] [ "" ] cudaCapability;

  nvcc-arch-sys = if cudaCapability-sys == null
  then null
  else get-nvcc-arch-from-cudaCapability cudaCapability-sys;

  options-default = {
    qlat-name = "";
    ngpu = ngpu-sys;
    nvcc-arch = nvcc-arch-sys;
    cudaCapabilities = cudaCapabilities-sys;
    cudaForwardCompat = false;
    use-cuda-software = false;
    use-grid-gpt = true;
    use-cps = true;
    use-cuda = false;
    use-cudasupport = false;
    use-cubaquad = true;
    use-clang = false;
    use-ucx = true;
    use-pypi = null;
  } // {
    ${if ngpu != null then "ngpu" else null} = ngpu;
    ${if cudaCapability != null then "nvcc-arch" else null} = get-nvcc-arch-from-cudaCapability cudaCapability;
    ${if cudaCapability != null then "cudaCapabilities" else null} = [ cudaCapability ];
    ${if cudaForwardCompat != null then "cudaForwardCompat" else null} = cudaForwardCompat;
  };

  mk-options = options:
  # Order the precedence:
  # options from arguments of this function.
  # options from arguments of this entire `q-pkgs.nix` file.
  # options-default from system parameters.
  let
    opts-0 = options-default // options;
    corrections-1 = opts: opts // {
      ${if opts.use-grid-gpt then "use-cubaquad" else null} = true;
    };
    corrections-2 = opts: opts // {
      ${if opts.use-cudasupport then "use-cuda" else null} = true;
    };
    corrections-3 = opts: opts // {
      ${if opts.use-cuda then "use-cuda-software" else null} = true;
      ${if opts.use-cuda then "use-clang" else null} = false;
    };
    corrections-4 = opts: opts // {
      ${if ! opts.use-cuda-software then "ngpu" else null} = 0;
      ${if ! opts.use-cuda-software then "nvcc-arch" else null} = null;
      ${if ! opts.use-cuda-software then "cudaCapabilities" else null} = [];
    };
    corrections-5 = opts: opts // {
      ${if opts.ngpu == "0" then "use-cuda-software" else null} = false;
      ${if opts.ngpu == "0" then "use-cuda" else null} = false;
      ${if opts.ngpu == "0" then "use-cudasupport" else null} = false;
    };
    opts = corrections-5 (corrections-4 (corrections-3 (corrections-2 (corrections-1 opts-0))));
    opts-1 =
      assert (opts.use-cudasupport -> opts.use-cuda);
      assert (opts.use-cuda -> opts.use-cuda-software);
      assert (opts.use-cuda-software -> opts.nvcc-arch != null);
      assert (opts.nvcc-arch != null -> opts.cudaCapabilities != []);
      assert (opts.use-grid-gpt -> opts.use-cubaquad);
      opts;
  in opts-1;

  mk-qlat-name = options:
  let
    opts = options;
  in opts.qlat-name
  + lib.optionalString (! opts.use-grid-gpt && ! opts.use-cps) "-std"
  + lib.optionalString (opts.use-grid-gpt && ! opts.use-cps) "-cpsless"
  + lib.optionalString (! opts.use-grid-gpt && opts.use-cps) "-gridless"
  + lib.optionalString opts.use-cuda-software "-cu"
  + lib.optionalString opts.use-cuda "da"
  + lib.optionalString opts.use-cudasupport "support"
  + lib.optionalString (! opts.use-cubaquad) "-cubaquadless"
  + lib.optionalString opts.use-clang "-clang"
  + lib.optionalString (! opts.use-ucx) "-ucxless"
  + lib.optionalString (opts.use-pypi != null) "-pypi"
  ;

  mk-overlay = options: final: prev: let
    opts = mk-options options;
    #
    pkgs = prev;
    #
    call-pkg = prev.callPackage;
    py-call-pkg = python3.pkgs.callPackage;
    #
    use-gitee = use-gitee-wd;
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
    qlat-cc = (if ! opts.use-clang
    then { inherit (pkgs) gcc; }
    else { inherit (pkgs) clang; inherit openmp; }
    )
    //
    { inherit (pkgs) pkg-config; }
    //
    (if opts.use-cuda-software then {
      inherit (pkgs.cudaPackages)
      cuda_nvcc
      cuda_cccl
      cuda_cudart
      cuda_profiler_api
      libcufft
      libcublas
      libnpp
      cudnn
      cudatoolkit
      ;
    } else {}
    );
    #
    # ollama = n-pkgs.ollama;
    ollama = pkgs.ollama;
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
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
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
      in lib.optionals opts.use-ucx [ "--with-ucx=${lib.getDev ucx-mt-dev}" ]
      ++ [
        (lib.withFeatureAs cudaSupport "cuda" (lib.getOutput "include" cudaPackages.cuda_cudart))
        (lib.withFeatureAs cudaSupport "cuda-libdir" "${lib.getLib cudaPackages.cuda_cudart}/lib")
      ]);
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
          doCheck = true;
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
          nativeBuildInputs = (py-prev.nativeBuildInputs or []) ++ [
            prev.jaxlib
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
        accelerate = if ! opts.use-cudasupport
        then prev.accelerate
        else prev.accelerate.overridePythonAttrs (py-prev: {
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
        gvar = pkgs.python3.pkgs.callPackage ./gvar.nix {};
        vegas = pkgs.python3.pkgs.callPackage ./vegas.nix { gvar = gvar; };
        lsqfit = pkgs.python3.pkgs.callPackage ./lsqfit.nix { gvar = gvar; vegas = vegas; };
        corrfitter = pkgs.python3.pkgs.callPackage ./corrfitter.nix { lsqfit = lsqfit; gvar = gvar; };
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
      cpuinfo = cpuinfo-sys;
    };
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
      inherit
      qlat-py
      qlat-nixgl
      qlat-stdenv
      mpi
      nvidia_x11_bin
      ;
    } // qlat-dep-pkgs // qlat-dep-pkgs-extra;
    qlat-tests = pkgs.buildEnv {
      name = "qlat-tests${qlat-name}";
      paths = builtins.attrValues qlat-tests-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
    };
    qlat-env = pkgs.buildEnv {
      name = "qlat-env${qlat-name}";
      paths = builtins.attrValues qlat-pkgs;
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
      ignoreCollisions = true;
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
      extraOutputsToInstall = [ "bin" "dev" "static" "man" "doc" "info" ];
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
      lsqfit
      corrfitter
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
      sentencepiece
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
      jupyterlab
      ollama
      ;
    }
    // (if is-linux then {
      inherit (ps)
      jupyterhub
      jupyterhub-systemdspawner
      xformers
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
    qlat-jhub-env = pkgs.buildEnv {
      name = "qlat-jhub-env${qlat-name}";
      paths = builtins.attrValues ({
        inherit
        qlat-jhub-py
        qlat-nixgl
        qlat-stdenv
        mpi
        nvidia_x11_bin
        ollama
        ;
        inherit (pkgs)
        bashInteractive
        bash-completion
        stdenvNoCC
        coreutils
        openssh
        findutils
        clang-tools
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
        texliveFull
        pipx
        twine
        poppler_utils
        file
        zip
        unzip
        ;
      }
      // (if is-linux then {
        inherit (pkgs)
        linux-pam
        ;
      } else {}
      )
      // qlat-cc
      // qlat-dep-pkgs
      // qlat-dep-pkgs-extra
      );
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
      ignoreCollisions = true;
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
      extraOutputsToInstall = [ "bin" "dev" "static" "man" "doc" "info" ];
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
    inherit qlat_utils qlat qlat_grid qlat_cps;
    inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt qlat-examples-py-cps;
    inherit qlat_docs qlat_pypipkgs;
    inherit qlat-py qlat-tests qlat-env qlat-sh qlat-fhs;
    inherit qlat-jhub-py qlat-jhub-env qlat-jhub-sh qlat-jhub-fhs;
    inherit qlat-pkgs q-pkgs;
    qlat-options = opts;
    qlat-use-pypi = opts.use-pypi;
  };

  mk-q-pkgs = options: let
    opts = mk-options options;
    qlat-name = mk-qlat-name opts;
    pkgs = import-nixpkgs-wd {
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
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
    "q-pkgs${pkgs.qlat-name}" = pkgs.q-pkgs;
    "pkgs${pkgs.qlat-name}" = pkgs;
    #
    "qlat-env${pkgs.qlat-name}" = pkgs.qlat-env;
    "qlat-tests${pkgs.qlat-name}" = pkgs.qlat-tests;
  };

  options-list = [
    {}
    { use-cps = false; use-grid-gpt = false; }
    { use-cps = false; }
    { use-grid-gpt = false; }
    { use-cuda-software = true; }
    { use-cuda = true; }
    { use-cudasupport = true; }
    { use-ucx = false; }
    { use-clang = true; }
    { use-pypi = version-pypi; }
    #
    { use-clang = true; use-ucx = false; }
    { use-cuda = true; use-ucx = false; }
    { use-grid-gpt = false; use-cubaquad = false; }
    { use-grid-gpt = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-ucx = false; }
    { use-cps = false; use-grid-gpt = false; use-clang = true; use-ucx = false; }
    #
    { use-cps = false; use-ucx = false; }
    { use-cps = false; use-clang = true; use-ucx = false; }
    { use-cps = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-cuda-software = true; }
    { use-cps = false; use-grid-gpt = false; use-cuda = true; }
    { use-cps = false; use-grid-gpt = false; use-cudasupport = true; }
  ];

  qlat-name-list = lib.lists.unique (builtins.map mk-qlat-name (builtins.map mk-options options-list));

  qlat-name-list-file-from-str = builtins.toFile "qlat-name-list"
  (builtins.foldl' (s: v: s + "${v}\n") "" qlat-name-list);

  qlat-name-list-file = runCommand
  "qlat-name-list"
  {}
  ''
    cp -v ${qlat-name-list-file-from-str} $out
  '';

  q-pkgs-list = builtins.map mk-q-pkgs options-list;

  all-q-pkgs = builtins.foldl' (s: v: s // v) {} q-pkgs-list;

  all-qlat-env = builtins.foldl' (s: v: s // { "qlat-env${v}" = all-q-pkgs."qlat-env${v}"; }) {} qlat-name-list;

  all-qlat-tests = builtins.foldl' (s: v: s // { "qlat-tests${v}" = all-q-pkgs."qlat-tests${v}"; }) {} qlat-name-list;

in all-q-pkgs // {
  inherit mk-q-pkgs;
  inherit mk-overlay;
  inherit options-list qlat-name-list qlat-name-list-file-from-str;
  inherit qlat-name-list-file;
  inherit all-qlat-env all-qlat-tests;
}
