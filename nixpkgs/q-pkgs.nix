{
  nixpkgs ? import ./nixpkgs.nix,
  ngpu ? null, # adjust with desired number of GPUs. E.g. "2"
  cudaCapability ? null, # adjust with desired cudaCapability. E.g. "8.6"
  cudaForwardCompat ? null, # adjust with desired cudaForwardCompat. E.g. false
}:

let

  o-pkgs = nixpkgs {
    config.allowUnfree = true;
  };

  lib = o-pkgs.lib;

  runCommandLocal = o-pkgs.runCommandLocal;

  runCommand = o-pkgs.runCommand;

  version-pypi = "0.76";

  nixgl-src = builtins.fetchTarball "https://github.com/jinluchang/nixGL/archive/0666d975fb80a7e8bedb4caafff8f113c5b14072.tar.gz";

  nixgl = (import nixgl-src {}).auto.nixGLDefault;

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

  cudaCapability-sys = if ngpu-sys == "0"
  then null
  else let
    nvidia_x11 = o-pkgs.linuxPackages.nvidia_x11;
  in builtins.head (builtins.match
  "(.*)\n"
  (builtins.readFile (runCommandLocal
  "impure-cuda-capability-file"
  {
    time = builtins.currentTime;
  }
  ''
    ${nixgl}/bin/nixGL ${nvidia_x11.bin}/bin/nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
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
      cudatoolkit
      ;
    } else {}
    );
    #
    ucx = pkgs.ucx.override {
      enableCuda = opts.use-cuda-software;
    };
    ucx-dev = pkgs.buildEnv {
      name = "ucx-dev";
      paths = [ ucx ];
      extraOutputsToInstall = [ "out" "bin" "dev" "static" "man" "doc" "info" ];
    };
    mpi = (pkgs.mpi.overrideAttrs (final: prev: {
      configureFlags = prev.configureFlags ++
      (let
        cudaPackages = pkgs.cudaPackages;
      in lib.optionals opts.use-ucx [ "--with-ucx=${lib.getDev ucx-dev}" ]
      ++ lib.optionals opts.use-cuda-software [ "--with-cuda-libdir=${cudaPackages.cuda_cudart.stubs}/lib" ]
      );
    })).override { cudaSupport = opts.use-cuda-software; };
    python3 = pkgs.python3.override {
      packageOverrides = final: prev: rec {
        mpi4py = prev.mpi4py.overridePythonAttrs (py-prev: {
          doCheck = true;
          nativeBuildInputs = (py-prev.nativeBuildInputs or [])
          ++ lib.optionals opts.use-cuda-software [
            qlat-nixgl
            pkgs.which
          ];
          preInstallCheck = lib.optionalString opts.use-cuda-software ''
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
        jax = prev.jax.overridePythonAttrs (py-prev: {
          doCheck = ! opts.use-cudasupport;
          nativeBuildInputs = (py-prev.nativeBuildInputs or [])
          ++ lib.optionals opts.use-cudasupport [
            prev.jaxlib
          ];
        });
        accelerate = prev.accelerate.overridePythonAttrs (py-prev: {
          doCheck = ! opts.use-cudasupport;
        });
        gvar = pkgs.python3.pkgs.callPackage ./gvar.nix {};
        vegas = pkgs.python3.pkgs.callPackage ./vegas.nix { gvar = gvar; };
        lsqfit = pkgs.python3.pkgs.callPackage ./lsqfit.nix { gvar = gvar; vegas = vegas; };
      };
    };
    #
    grid-lehner-c-lime = if opts.use-cps then qio else c-lime;
    #
    qlat-nixgl = if opts.use-cuda-software then nixgl else null;
    #
    cubaquad = if opts.use-cubaquad
    then call-pkg ./cubaquad.nix { stdenv = qlat-stdenv; }
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
      inherit mpi cubaquad qlat-eigen;
    } else if opts.use-grid-gpt && ! opts.use-cps then {
      inherit mpi cubaquad qlat-eigen grid-lehner c-lime;
    } else if ! opts.use-grid-gpt && opts.use-cps then {
      inherit mpi cubaquad qlat-eigen cps qmp qio;
    } else {
      inherit mpi cubaquad qlat-eigen cps qmp qio grid-lehner;
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
      inherit qlat-py;
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
      ollama
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
      ;
    }
    //
    qlat-py-pkgs
    //
    (if opts.use-cuda-software then {
      inherit (ps) pycuda;
    } else {}
    )
    ));
    qlat-jhub-env = pkgs.buildEnv {
      name = "qlat-jhub-env${qlat-name}";
      paths = builtins.attrValues ({
        inherit qlat-jhub-py mpi;
        inherit (pkgs)
        bashInteractive
        bash-completion
        coreutils
        openssh
        linux-pam
        findutils
        clang-tools
        git
        gnumake
        zlib
        hwloc
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
        ollama
        ;
      }
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
    inherit python3 mpi openmp ucx;
    inherit c-lime qmp qio cps cubaquad grid-lehner gpt-lehner;
    inherit (opts) use-pypi;
    inherit qlat_utils qlat qlat_grid qlat_cps;
    inherit qlat-examples-cpp qlat-examples-cpp-grid qlat-examples-py qlat-examples-py-gpt qlat-examples-py-cps;
    inherit qlat_docs qlat_pypipkgs;
    inherit qlat-py qlat-tests qlat-env qlat-sh qlat-fhs;
    inherit qlat-jhub-py qlat-jhub-env qlat-jhub-sh qlat-jhub-fhs;
    inherit qlat-pkgs q-pkgs;
  };

  mk-q-pkgs = options: let
    opts = mk-options options;
    pkgs = nixpkgs {
      config = {
        allowUnfree = opts.use-cuda-software;
        cudaSupport = opts.use-cudasupport;
        ${if opts.use-cudasupport then "cudaCapabilities" else null} = opts.cudaCapabilities;
        ${if opts.use-cudasupport then "cudaForwardCompat" else null} = opts.cudaForwardCompat;
      };
      overlays = [
        (mk-overlay options)
      ];
    };
  in {
    "qlat-pkgs${pkgs.qlat-name}" = pkgs.qlat-pkgs;
    "q-pkgs${pkgs.qlat-name}" = pkgs.q-pkgs;
    "pkgs${pkgs.qlat-name}" = pkgs;
    #
    "qlat-py${pkgs.qlat-name}" = pkgs.qlat-py;
    "qlat-env${pkgs.qlat-name}" = pkgs.qlat-env;
    "qlat-sh${pkgs.qlat-name}" = pkgs.qlat-sh;
    "qlat-fhs${pkgs.qlat-name}" = pkgs.qlat-fhs;
    #
    "qlat-jhub-py${pkgs.qlat-name}" = pkgs.qlat-jhub-py;
    "qlat-jhub-env${pkgs.qlat-name}" = pkgs.qlat-jhub-env;
    "qlat-jhub-sh${pkgs.qlat-name}" = pkgs.qlat-jhub-sh;
    "qlat-jhub-fhs${pkgs.qlat-name}" = pkgs.qlat-jhub-fhs;
    #
    "qlat-tests${pkgs.qlat-name}" = pkgs.qlat-tests;
  };

  options-list = [
    {}
    { use-grid-gpt = false; use-cps = false; }
    { use-grid-gpt = false; }
    { use-cps = false; }
    { use-cuda-software = true; }
    { use-cuda = true; }
    { use-cudasupport = true; }
    { use-grid-gpt = false; use-cubaquad = false; }
    { use-grid-gpt = false; use-clang = true; }
    { use-ucx = false; }
    { use-pypi = version-pypi; }
    #
    { use-grid-gpt = false; use-cps = false; use-ucx = false; }
    { use-grid-gpt = false; use-cps = false; use-clang = true; use-ucx = false; }
    #
    # { use-clang = true; } # not working yet
    # { use-clang = true; use-ucx = false; } # not working yet
    # { use-cps = false; use-clang = true; use-ucx = false; } # not working yet
    # { use-cps = false; use-clang = true; } # not working yet
    #
    { use-grid-gpt = false; use-cps = false; use-cuda-software = true; }
    { use-grid-gpt = false; use-cps = false; use-cuda = true; }
    { use-grid-gpt = false; use-cps = false; use-cudasupport = true; }
  ];

  qlat-name-list = lib.lists.unique (map mk-qlat-name (map mk-options options-list));

  qlat-name-list-file-from-str = builtins.toFile "qlat-name-list"
  (builtins.foldl' (s: v: s + "\"${v}\"\n") "" qlat-name-list);

  qlat-name-list-file = runCommand
  "qlat-name-list"
  {}
  ''
    cp -v ${qlat-name-list-file-from-str} $out
  '';

  q-pkgs = (builtins.foldl' (s: v: s // v) {}
  (builtins.map mk-q-pkgs options-list))
  // {
    inherit mk-q-pkgs;
    inherit mk-overlay;
    inherit options-list qlat-name-list qlat-name-list-file-from-str;
    inherit qlat-name-list-file;
  };

in q-pkgs
