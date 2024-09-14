{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
, qlat_grid
, grid-lehner
, git
, bash
, time
, mpi
, mpiCheckPhaseHook
, openssh
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat-examples-cpp-grid${qlat-name}";
  version = builtins.readFile ../VERSION + "current";

  pyproject = false;

  src = ../examples-cpp;

  enableParallelBuilding = true;

  inherit stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
    grid-lehner
    time
    mpi
    mpiCheckPhaseHook
    openssh
  ];

  propagatedBuildInputs = [
    grid-lehner
  ];

  dependencies = [
    qlat
    qlat_grid
  ];

  preConfigure = ''
    export
    echo
    ls -l
    echo
    pwd
    echo
    #
    export LD_LIBRARY_PATH="$(python3 -m qlat qlat-config --LD_LIBRARY_PATH)"
    export mpi_options="--oversubscribe --bind-to none $mpi_options"
    export SHELL=${bash}/bin/bash
    #
    echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    echo mpi_options=$mpi_options
    echo SHELL=$SHELL
    echo NIX_BUILD_CORES=$NIX_BUILD_CORES
    echo NIX_BUILD_TOP=$NIX_BUILD_TOP
    echo
	#
    export num_proc=$((NIX_BUILD_CORES / 4 + 1))
    echo num_proc=$num_proc
    #
    patchShebangs --build */run.sh
    echo
    #
    make run-grid -j$num_proc SHELL=$SHELL
    #
    echo
    for i in * ; do
      if [ -f "$i"/build/log.check.txt.new ] ; then
        if diff "$i"/build/log.check.txt "$i"/build/log.check.txt.new ; then
          echo "$i" passed
        else
          echo
          cat "$i"/build/log.full
          echo
          echo "$i" failed
          echo
          false
        fi
      fi
    done
    echo
    #
    pwd
  '';

  dontBuild = true;
  dontInstall = true;

}
