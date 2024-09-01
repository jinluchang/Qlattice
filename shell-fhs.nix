{
  pkgs ? import <nixpkgs> {},
}:

(pkgs.buildFHSEnv {
  name = "qlat-build-env";
  targetPkgs = pkgs: (with pkgs; [
    gcc
    mpi
    pkg-config
    automake
    autoconf
    rsync
    zlib
    gsl
    fftw
    fftwFloat
    openssl
    (python3.withPackages (ps: with ps; [
      meson
      ninja
      meson-python
      numpy
      scipy
      sympy
      jax
      jaxlib
      mpi4py
      psutil
      cython
      pkgconfig
      h5py
      matplotlib
      build
      wheel
    ]))
  ]);
  runScript = "bash";
  extraOutputsToInstall = [ "man" "doc" "info" "dev" "static" ];
}).env
