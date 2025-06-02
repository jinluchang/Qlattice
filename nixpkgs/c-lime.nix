{ stdenv
, lib
, fetchzip
, openmp ? null
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  pname = "c-lime";
  version = "1.3.2";

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/c-lime" else "https://github.com/usqcd-software/c-lime";
    ref = "refs/tags/c-lime1-3-2";
  };

  enableParallelBuilding = true;

  propagatedBuildInputs = [
  ]
  ++ lib.optional stdenv.cc.isClang openmp;

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
  '';

}
