{ lib
, buildPythonPackage
, fetchFromGitHub
, meson-python
, ninja
, numpy
, torch
, pytest
, use-gitee ? null
}:

let
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in buildPythonPackage rec {

  pname = "qcd_ml";
  version = "1.0.1";
  pyproject = true;

  src = fetchFromGitHub {
    owner = "daknuett";
    repo = "qcd_ml";
    rev = "v${version}";
    hash = "sha256-gofVim9oC13o7Q1Va/S+5MBv51AfkZBSRSrvc3S7Ou4=";
  };

  build-system = [
    meson-python
    ninja
  ];

  dependencies = [
    numpy
    torch
  ];

  test-dependencies = [
    pytest
  ];

  meta = {
    description = "Machine learning tools for QCD";
    homepage = "https://github.com/daknuett/qcd_ml";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
