{ stdenv
, fetchzip
, mpi
, qmp
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  pname = "qio";
  version = "3.0.0";

  src = fetchzip {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/qio-${version}.tar.gz" else "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/qio-${version}.tar.gz";
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
