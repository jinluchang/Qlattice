{ fetchPypi
, stdenv
, buildPythonPackage
, qlat
, git
, bash
, time
, mpi
, mpiCheckPhaseHook
, openssh
, qlat-name ? ""
}:

buildPythonPackage rec {

  pname = "qlat-examples-cpp${qlat-name}";
  version = "${../VERSION}-current";

  pyproject = false;

  src = ../examples-cpp;

  enableParallelBuilding = true;

  inherit stdenv;

  build-system = [
    qlat
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
    patchShebangs --build */run.sh
    echo
    #
    make run -j1 SHELL=$SHELL
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
