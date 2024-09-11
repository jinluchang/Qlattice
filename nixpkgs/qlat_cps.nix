{ fetchPypi
, stdenv
, buildPythonPackage
, qlat
, cps
, git
, is-pypi-src ? true
, qlat-name ? ""
}:

buildPythonPackage rec {

  pname = "qlat_cps${qlat-name}";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    pname = "qlat_cps";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-X0hCGuFUgsvZ9AKYr7JhhxgM5hCp3zrbHYGpz3zVVj0=";
  };

  version-local = "${../VERSION}-current";
  src-local = ../qlat-cps;

  enableParallelBuilding = true;

  inherit stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
  ];

  propagatedBuildInputs = [
    cps
  ];

  dependencies = [
    qlat
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
