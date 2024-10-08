# Used to set up nix shell.
#
# With this file, you can build with:
#
# $ nix-shell
# $ ./build default-nix
#
# May need to have
# services.envfs.enable = true;
# in /etc/nixos/configuration.nix

{
  pkgs ? import <nixpkgs> {},
}:

pkgs.mkShell {
  name = "qlat-build-sh";
  packages = with pkgs; [
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
  ];
}
