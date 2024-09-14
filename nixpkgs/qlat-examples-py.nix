{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
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

  pname = "qlat-examples-py${qlat-name}";
  version = builtins.readFile ../VERSION + "current";

  pyproject = false;

  src = ../examples-py;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
    time
    mpi
    mpiCheckPhaseHook
    openssh
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
  ];

  dependencies = [
    qlat
  ];

  preConfigure = ''
    export
    echo
    ls -l
    echo
    pwd
    echo
    #
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
    patchShebangs --build */run.py
    echo
    #
    make update-sources SHELL=$SHELL
    echo
    make run -j$num_proc SHELL=$SHELL
    echo
    #
    for i in *.p ; do
      if [ -d "$i" ] ; then
        if diff "$i"/log.check.txt "$i"/log.check.txt.new ; then
          echo "$i" passed
        else
          echo
          cat "$i"/log.full.txt
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
