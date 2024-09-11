{ fetchPypi
, stdenv
, buildPythonPackage
, qlat
, qlat_grid
, gpt-lehner
, git
, bash
, time
, mpi
, mpiCheckPhaseHook
, openssh
}:

buildPythonPackage rec {

  pname = "qlat-examples-py-gpt";
  version = "${../VERSION}-current";

  pyproject = false;

  src = ../examples-py;

  enableParallelBuilding = false;

  inherit stdenv;

  build-system = [
    qlat
    qlat_grid
    gpt-lehner
  ];

  nativeBuildInputs = [
    git
    time
    mpi
    mpiCheckPhaseHook
    openssh
  ];

  propagatedBuildInputs = [
  ];

  dependencies = [
    qlat
    qlat_grid
    gpt-lehner
  ];

  preConfigure = ''
    export
    echo
    ls -l
    echo
    pwd
    echo
    #
    export mpi_options="--oversubscribe $mpi_options"
    export SHELL=${bash}/bin/bash
    export num_proc=1
    #
    echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    echo mpi_options=$mpi_options
    echo SHELL=$SHELL
    echo NIX_BUILD_CORES=$NIX_BUILD_CORES
    echo NIX_BUILD_TOP=$NIX_BUILD_TOP
    echo num_proc=$num_proc
    echo
    #
    make update-sources
    echo
    make run-gpt -j1 SHELL=$SHELL
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
