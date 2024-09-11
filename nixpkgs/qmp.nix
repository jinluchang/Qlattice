{ stdenv
, lib
, fetchzip
, mpi
, openmp ? null
}:

stdenv.mkDerivation rec {

  pname = "qmp";
  version = "2.5.4";

  src = fetchzip {
    url = "http://usqcd-software.github.io/downloads/qmp/qmp-${version}.tar.gz";
    hash = "sha256-bUeOtdH2g36SdGK13WVBJ/AiZnp+p8oF31fIShqGKF8=";
  };

  enableParallelBuilding = true;

  propagatedBuildInputs = [
    mpi
  ]
  ++ lib.optional stdenv.cc.isClang openmp;

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
    export CC="mpicc"
    export CXX="mpic++"
  '';

  configureFlags = [
    "--with-qmp-comms-type=mpi"
  ];

}
