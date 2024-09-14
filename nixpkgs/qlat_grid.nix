{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
, grid-lehner
, git
, is-pypi-src ? true
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat_grid${qlat-name}";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    pname = "qlat_grid";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-ZnsO4Trkihq9fP8Y3viJj14IyFQgXlx99WcwORV2rMY=";
  };

  version-local = builtins.readFile ../VERSION + "current";
  src-local = ../qlat-grid;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
    grid-lehner
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    grid-lehner
  ];

  dependencies = [
    qlat
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = ''
    export OMPI_CXX=c++
    export OMPI_CC=cc
    #
    export
  '';

}
