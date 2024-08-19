# Used to set up nix shell.
#
# With this file, you can build with:
#
# $ nix-shell
# $ ./build default-nix
#

let
  pkgs = import <nixpkgs> {};
in pkgs.mkShell {
  packages = with pkgs; [
    gcc
    mpi
    pkg-config
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
    ]))
  ];
}
