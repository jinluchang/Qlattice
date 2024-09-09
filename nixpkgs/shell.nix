# Used to set up nix shell.
#
# With this file, you can build with:
#
# $ nix-shell
# $ ./build default-nix

{
  pkgs ? import <nixpkgs> {}
}:

let
  local-pkgs = import ./default.nix;
  env = pkgs.mkShell {
    name = "qlat-build-sh";
    packages = with pkgs; [
      pkg-config
      local-pkgs.qlat-full-cuda-local
    ];
    inputsFrom = with pkgs; [
      fftw
      fftwFloat
      hdf5-cpp
      gsl
      cuba
      zlib
      local-pkgs.qlat-eigen
      local-pkgs.qlat-full-cuda-local
    ];
  };
in env
