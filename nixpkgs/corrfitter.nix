{ lib
, buildPythonPackage
, fetchPypi
# build-system
, setuptools
, setuptools-scm
, lsqfit
, gvar
, numpy
# dependencies
# tests
, hypothesis
#
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in buildPythonPackage rec {

  pname = "corrfitter";
  version = "8.2";
  pyproject = true;

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/${pname}" else "https://github.com/gplepage/${pname}";
    ref = "refs/tags/v${version}";
  };

  build-system = [
    setuptools
    setuptools-scm
    lsqfit
    gvar
    numpy
  ];

  dependencies = [
    lsqfit
    gvar
    numpy
  ];

  meta = {
    description = "Utilities for fitting correlators in lattice QCD.";
    homepage = "https://github.com/gplepage/corrfitter";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
