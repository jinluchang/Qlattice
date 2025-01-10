# hello.nix
{ stdenv
, fetchzip
}:

stdenv.mkDerivation rec {

  pname = "Cuba";
  version = "4.2.2";

  src = fetchzip {
    url = "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/Cuba-${version}.tar.gz";
    hash = "sha256-4mPW/mKWN1VPBCkwayif9m1Qd5cXm9Y6PLtUrD/8knM=";
  };

  enableParallelBuilding = false;

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
  '';

}
