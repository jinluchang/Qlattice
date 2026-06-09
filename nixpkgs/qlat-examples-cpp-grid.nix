{ stdenv
, lib
, config
, buildPythonPackage
, qlat
, qlat_grid
, grid-lehner
, git
, which
, rsync
, bash
, time
, mpi
, mpiCheckPhaseHook
, openssh
, qlat-src
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, ngpu ? "1"
, version ? "current"
}:

let

  src = qlat-src.examples-cpp-grid;

in (buildPythonPackage.override { stdenv = stdenv; }) rec {

  inherit version src;

  pname = "qlat-examples-cpp-grid${qlat-name}";

  pyproject = false;

  enableParallelBuilding = true;

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
    which
    rsync
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
      #
      source ${qlat}/bin/cuda-mpi-qlat.sh echo
      #
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
      #
      echo "CXX=$CXX"
      echo "CXXFLAGS=$CXXFLAGS"
      echo "LDFLAGS=$LDFLAGS"
      #
      echo "MPICXX=$MPICXX"
      echo "OMPI_CXX=$OMPI_CXX"
      #
      export NGPU=${ngpu}
    '';
    cpu_extra = ''
    '';
    extra = (if cudaSupport then gpu_extra else cpu_extra) + (
      if cudaSupport && lib.hasInfix "cuda" qlat-name then ''
        export num_proc=$((NIX_BUILD_CORES / 16 + 1))
      '' else ''
        export num_proc=$((NIX_BUILD_CORES / 4 + 1))
      ''
    );
  in extra + ''
    #
    # CXX_ARR=($(grid-config --cxx))
    # export CXX="''${CXX_ARR[0]}"
    # export CXXFLAGS="''${CXX_ARR[@]:1} $CXXFLAGS"
    # export LDFLAGS="''${CXX_ARR[@]:1} $LDFLAGS"
    #
    export OMP_NUM_THREADS=2
    echo
    ls -l
    echo
    pwd
    echo
    #
    export SHELL=${bash}/bin/bash
    #
    echo
    export
    echo
    #
    echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    echo SHELL=$SHELL
    echo NIX_BUILD_CORES=$NIX_BUILD_CORES
    echo NIX_BUILD_TOP=$NIX_BUILD_TOP
    echo
    #
    echo num_proc=$num_proc
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
        fi
      fi
    done
    echo
    failed=false
    for i in * ; do
      if [ -f "$i"/build/log.check.txt.new ] ; then
        if diff "$i"/build/log.check.txt "$i"/build/log.check.txt.new >/dev/null 2>&1 ; then
          echo "$i" passed
        else
          echo "$i" failed
          failed=true
        fi
      fi
    done
    if $failed ; then
      false
    fi
    echo
    #
    pwd
    #
    # rm -rfv ./*/build/results
    rm -rf ./*/build
    #
    mkdir -p "$out/share/qlat/examples-cpp-grid"
    rsync -a --delete . "$out/share/qlat/examples-cpp-grid"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
