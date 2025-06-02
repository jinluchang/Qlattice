{ stdenv
, lib
, fetchzip
, mpi
, openmp ? null
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  pname = "qmp";
  version = "2.5.4";

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/qmp" else "https://github.com/usqcd-software/qmp";
    ref = "refs/tags/qmp2-5-4";
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
