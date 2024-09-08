# hello.nix
{ stdenv
, fetchFromGitHub
, mpi
, qmp
, qio
, gmp
, mpfr
, zlib
, fftw
, fftwFloat
, gsl
, which
, git
, c-lime
}:

stdenv.mkDerivation rec {

  pname = "cps";
  version = "d3c8dd5e8a3ea6a315fd2fba963bf32585ed6331";

  src = fetchFromGitHub {
    owner = "RBC-UKQCD";
    repo = "CPS_public";
    rev = version;
    hash = "sha256-8+11yToZRcexNuGQHXf5qcl6cs1/Ru9xQGt1MoKzna4=";
  };

  enableParallelBuilding = true;

  nativeBuildInputs = [
    mpi
    git
    which
  ];

  propagatedBuildInputs = [
    mpi
    c-lime
    qmp
    qio
    mpfr
    zlib
    fftw
    fftwFloat
    gsl
    gmp
  ];

  preConfigure = ''
    export CFLAGS="-fPIC -w -Wno-psabi"
    export CXXFLAGS="-fPIC -w -Wno-psabi"
    export CC="mpicc"
    export CXX="mpic++"
    cd cps_pp
    pwd
    ls -l
  '';

  configureFlags = [
    "--enable-gmp"
    "--enable-mpfr"
    "--enable-qmp=${qmp}"
    "--enable-qio=${qio}"
  ];

  postBuild = ''
    pwd
    ls -l
    mkdir -p "$out/include"
    mkdir -p "$out/lib"
    cp -rpv include/* "$out/include/"
    cp -rpv *.h "$out/include/"
    cp -rpv libcps.a "$out/lib/"
  '';

  dontInstall = true;

}
