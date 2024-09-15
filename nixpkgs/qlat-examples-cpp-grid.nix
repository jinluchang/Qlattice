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
, NVCC_ARCH ? "sm_86"
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat-examples-cpp-grid${qlat-name}";
  version = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  pyproject = false;

  src = ../examples-cpp;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

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
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    qlat_grid
  ];

  dependencies = [
    qlat
    qlat_grid
  ];

  preConfigure = let
    gpu_extra = ''
      pwd
      cp -pv "${../qcore/bin/NVCC.py}" "$PWD/NVCC.py"
      patchShebangs --build "$PWD/NVCC.py"
      #
      export NVCC_OPTIONS="-w -std=c++14 -arch=${NVCC_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing" # -D__DEBUG_VECUTILS__
      export QLAT_CXX="$PWD/NVCC.py -ccbin c++ $NVCC_OPTIONS"
      export QLAT_MPICXX="$PWD/NVCC.py -ccbin mpic++ $NVCC_OPTIONS"
      export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
      export QLAT_LDFLAGS="--NVCC-link" # --shared
      #
      export MPICXX="$QLAT_MPICXX"
      export CXX="$MPICXX"
      export CXXFLAGS="$QLAT_CXXFLAGS"
      export LDFLAGS="$QLAT_LDFLAGS"
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in ''
    export OMPI_CXX=c++
    export OMPI_CC=cc
    #
    CXX_ARR=($(grid-config --cxx))
    export CXX="''${CXX_ARR[0]}"
    export CXXFLAGS="''${CXX_ARR[@]:1} $CXXFLAGS"
    export LDFLAGS="''${CXX_ARR[@]:1} $LDFLAGS"
    #
    export
  '' + extra + ''
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
