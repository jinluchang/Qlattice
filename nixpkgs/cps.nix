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
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  env.NIX_CFLAGS_COMPILE = lib.concatStringsSep " " [
    "-Wno-error=int-conversion"
    "-Wno-error=incompatible-pointer-types"
    "-Wno-error=implicit-function-declaration"
    ];

  pname = "cps";
  version = "d3c8dd5e8a3ea6a315fd2fba963bf32585ed6331";

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/cps" else "https://github.com/RBC-UKQCD/CPS_public";
    ref = "master";
    rev = version;
  };

  # Darwin portability fixes for the pinned upstream CPS (shared with the
  # scripts/cps.sh flow):
  #   - only use Linux <sys/sysinfo.h>; portable printMem() fallback
  #   - rename cps_pp/VERSION -> VERSION.txt (it shadows libc++ <version> on
  #     case-insensitive APFS when cps_pp is on the include path)
  #   - portable 128-byte aligned allocation (no <malloc.h>/memalign on macOS)
  patches = [
    ./patches/cps-darwin-sysinfo.patch
    ./patches/cps-darwin-version-shadow.patch
    ./patches/cps-darwin-memalign.patch
  ];

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
