{ fetchPypi
, stdenv
, config
, lib
, buildPythonPackage
, cython
, meson-python
, pkg-config
, numpy
, psutil
, zlib
, eigen
, git
, openmp ? null
, is-pypi-src ? true
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat_utils${qlat-name}";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    pname = "qlat_utils";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-cfvie6jiE1VWyNPjKgQK26P6hpCgYr008WIIIpEipPA=";
  };

  version-local = builtins.readFile ../VERSION + "current";
  src-local = ../qlat-utils;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    meson-python
    pkg-config
    cython
    numpy
  ];

  nativeBuildInputs = [
    git
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    zlib
    eigen
  ]
  ++ lib.optional stdenv.cc.isClang openmp;

  dependencies = [
    meson-python
    pkg-config
    cython
    numpy
    psutil
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
