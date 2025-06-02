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
# tests
, hypothesis
#
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in buildPythonPackage rec {

  pname = "vegas";
  version = "6.3";
  pyproject = true;

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/${pname}" else "https://github.com/gplepage/${pname}";
    ref = "refs/tags/v${version}";
  };

  build-system = [
    setuptools
    numpy
    cython
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    gvar
  ];

  meta = {
    description = "Tools for adaptive multidimensional Monte Carlo integration.";
    homepage = "https://github.com/gplepage/vegas";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
