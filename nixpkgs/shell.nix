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
  pkgs ? import <nixpkgs> {}
}:

let
  local-pkgs = import ./default.nix;
in pkgs.mkShell {
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
    eigen
    local-pkgs.cuba
    local-pkgs.grid-lehner
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
      local-pkgs.qlat_utils
      local-pkgs.qlat
      local-pkgs.gpt-lehner
    ]))
  ];
}
