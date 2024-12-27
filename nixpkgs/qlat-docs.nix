{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
, qlat_cps
, qlat_grid
, gpt-lehner
, sphinx
, linkify-it-py
, myst-parser
, git
, which
, rsync
, bash
, time
, mpi
, mpiCheckPhaseHook
, openssh
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? ""
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat-docs${qlat-name}";
  version = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  pyproject = false;

  src = ../docs;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    qlat
    qlat_cps
    qlat_grid
    gpt-lehner
    sphinx
    linkify-it-py
    myst-parser
  ];

  nativeBuildInputs = [
    git
    time
    mpi
    mpiCheckPhaseHook
    openssh
    which
    rsync
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ++ lib.optionals cudaSupport [ nixgl ]
  ;

  propagatedBuildInputs = [
  ];

  dependencies = [
  ];

  preConfigure = let
    examples-cpp = ../examples-cpp;
    examples-cpp-grid = ../examples-cpp-grid;
    examples-py = ../examples-py;
    examples-py-cps = ../examples-py-cps;
    examples-py-gpt = ../examples-py-gpt;
    gpu_extra = ''
      which nixGL
      echo
      echo "run with nixGL"
      cat $(which nixGL) | grep -v 'exec ' | grep -v '^#!' > nix-gl.sh
      echo
      echo cat nix-gl.sh
      cat nix-gl.sh
      source nix-gl.sh
      echo
      echo $LD_LIBRARY_PATH
      echo
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    pwd
    echo
    rsync -a "${examples-cpp}/" ../examples-cpp/
    rsync -a "${examples-cpp-grid}/" ../examples-cpp-grid/
    rsync -a "${examples-py}/" ../examples-py/
    rsync -a "${examples-py-cps}/" ../examples-py-cps/
    rsync -a "${examples-py-gpt}/" ../examples-py-gpt/
    #
    ls -ld *
    ls -l ..
    ls -l ../examples-cpp/
    ls -l ../examples-py/
    #
    make html
    #
    mkdir -p "$out/share/doc/qlat"
    rsync -a --delete ./build/ "$out/share/doc/qlat/"
  '';

  dontBuild = true;
  dontInstall = true;

}
