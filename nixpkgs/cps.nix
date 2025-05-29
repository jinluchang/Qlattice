{ stdenv
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
, lib
}:

stdenv.mkDerivation rec {

  env.NIX_CFLAGS_COMPILE = lib.concatStringsSep " " [
    "-Wno-error=int-conversion"
    "-Wno-error=incompatible-pointer-types"
    "-Wno-error=implicit-function-declaration"
    ];

  pname = "cps";
  version = "d3c8dd5e8a3ea6a315fd2fba963bf32585ed6331";

  src = builtins.fetchGit {
    url = "https://github.com/RBC-UKQCD/CPS_public";
    ref = "master";
    rev = version;
  };

  enableParallelBuilding = true;

  nativeBuildInputs = [
    mpi
    git
    which
  ];

  propagatedBuildInputs = [
    mpi
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
