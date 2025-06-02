# hello.nix
{ stdenv
, fetchzip
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  pname = "Cuba";
  version = "4.2.2";

  src = fetchzip {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/Cuba-${version}.tar.gz" else "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/Cuba-${version}.tar.gz";
    hash = "sha256-4mPW/mKWN1VPBCkwayif9m1Qd5cXm9Y6PLtUrD/8knM=";
  };

  enableParallelBuilding = false;

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
  '';

}
