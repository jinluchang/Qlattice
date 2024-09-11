{ stdenv
, fetchzip
, mpi
, qmp
}:

stdenv.mkDerivation rec {

  pname = "qio";
  version = "3.0.0";

  src = fetchzip {
    url = "http://usqcd-software.github.io/downloads/qio/qio-${version}.tar.gz";
    hash = "sha256-IDAna0ffc2SE+2+uCOurGXF9+wqikkK8HggVztZX/xM=";
  };

  enableParallelBuilding = true;

  propagatedBuildInputs = [
    mpi
    qmp
  ];

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
    export CC="mpicc"
    export CXX="mpic++"
  '';

  configureFlags = [
    "--with-qmp"
    "--build=none"
    "--enable-largefile"
  ];

}
