{ lib
, buildPythonPackage
, fetchPypi
# build-system
, setuptools
, wheel
, numpy
, cython
, oldest-supported-numpy
# dependencies
, scipy
# tests
#
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in buildPythonPackage rec {

  env.NIX_CFLAGS_COMPILE = lib.concatStringsSep " " [
    "-Wno-error=int-conversion"
    "-Wno-error=incompatible-pointer-types"
    "-Wno-error=implicit-function-declaration"
    ];

  pname = "gvar";
  version = "13.1.6";
  pyproject = true;

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/${pname}" else "https://github.com/gplepage/${pname}";
    ref = "refs/tags/v${version}";
  };

  build-system = [
    setuptools
    wheel
    numpy
    cython
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    scipy
  ];

  meta = {
    description = "Utilities for manipulating correlated Gaussian random variables.";
    homepage = "https://github.com/gplepage/gvar";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
