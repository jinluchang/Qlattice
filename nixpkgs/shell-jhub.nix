{
  nixpkgs ? import <nixpkgs>
}:

let
  qlat-pkgs = import ./qlat-pkgs.nix { inherit nixpkgs; };
  many-qlat-pkgs-all = qlat-pkgs.many-qlat-pkgs-all;
  pkgs = many-qlat-pkgs-all.pkgs;
  python3-with-packages =
    pkgs.python3.withPackages (ps: with ps; [
      pkgs.qlat_utils
      pkgs.qlat
      pkgs.qlat_cps
      pkgs.qlat_grid
      pkgs.gpt-lehner
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
      pycuda
      pytools
      lz4
      torchvision
      torchaudio
      xformers
      jupyterlab
      jupyterhub
    ]);
    env = (pkgs.buildFHSEnv {
      name = "jhub-env";
      targetPkgs = pkgs: (with pkgs; [
        python3-with-packages
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
      ]);
      runScript = "bash";
      extraOutputsToInstall = [ "man" "doc" "info" "dev" "static" ];
    }).env;
in env
