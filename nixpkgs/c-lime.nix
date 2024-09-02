# hello.nix
{ stdenv
, fetchzip
}:

stdenv.mkDerivation rec {

  pname = "c-lime";
  version = "1.3.2";

  src = fetchzip {
    url = "http://usqcd-software.github.io/downloads/c-lime/lime-${version}.tar.gz";
    hash = "sha256-eOU5fSYRuE0eMLFstqVhfmwzrob0kfd524J+vy0b7ss=";
  };

  enableParallelBuilding = true;

  preConfigure = ''
    export CFLAGS="-fPIC -O2"
    export CXXFLAGS="-fPIC -O2"
  '';

}
