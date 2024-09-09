{ fetchPypi
, stdenv
, buildPythonPackage
, qlat
, grid-lehner
, git
, is-pypi-src ? true
}:

buildPythonPackage rec {

  pname = "qlat_grid";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    inherit pname;
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-ZnsO4Trkihq9fP8Y3viJj14IyFQgXlx99WcwORV2rMY=";
  };

  version-local = "${../VERSION}-current";
  src-local = ../qlat-grid;

  enableParallelBuilding = true;

  inherit stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
    grid-lehner
  ];

  propagatedBuildInputs = [
    grid-lehner
  ];

  dependencies = [
    qlat
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
