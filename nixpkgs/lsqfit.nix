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
, gvar
, vegas
# tests
, hypothesis
#
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in buildPythonPackage rec {

  pname = "lsqfit";
  version = "13.2.2";
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
    gvar
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    gvar
    vegas
  ];

  meta = {
    description = "Utilities for nonlinear least-squares fits.";
    homepage = "https://github.com/gplepage/lsqfit";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
