# hello.nix
{ stdenv
, fetchurl
, fetchFromGitHub
, zlib
, eigen
, openssl
, mpi
, gsl
, hdf5
, gmp
, mpfr
, fftw
, fftwFloat
, git
, autoconf
, automake
}:

stdenv.mkDerivation rec {

  pname = "Grid-lehner";
  version = "9f89486df5e65c873308df23240a3b826c257d76";

  src = fetchFromGitHub {
    owner = "lehner";
    repo = "Grid";
    rev = version;
    hash = "sha256-2qdp/HcLEB/lbQ5ipm/3FGW5gofkLfcJj6d4HYtP0fs=";
  };

  enableParallelBuilding = true;

  buildInputs = [
    mpi
    zlib
    eigen
    fftw
    fftwFloat
    gsl
    openssl
    hdf5
    gmp
    mpfr
    autoconf
    automake
  ];

  preConfigure = let
    eigen-file-name = "eigen-3.3.7.tar.bz2";
    eigen-src = fetchurl {
      url = "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2";
      hash = "sha256-aFrfFL2OnAFbeAl8HcIvLwE0N1bxlqzcdqZ44a41LhE=";
    };
  in ''
    echo "-- deploying Eigen source..."
    cp -pv '${eigen-src}' '${eigen-file-name}'
    bash ./scripts/update_eigen.sh '${eigen-file-name}'
    echo '-- generating Make.inc files...'
    bash ./scripts/filelist
    echo '-- generating configure script...'
    autoreconf -fvi
    echo "set configure options"
    configureFlagsArray+=(
      "--enable-simd=AVX2"
      "--enable-alloc-align=4k"
      "--enable-comms=mpi-auto"
      "--enable-gparity=no"
    )
  '';

}
