{ stdenv
, lib
, config
, buildPythonPackage
, qlat
, qlat_cps
, qlat_grid
, gpt-lehner
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

  src = qlat-src.examples-py-cps;

in (buildPythonPackage.override { stdenv = stdenv; }) rec {

  inherit version src;

  pname = "qlat-examples-py-cps${qlat-name}";

  pyproject = false;

  enableParallelBuilding = true;

  build-system = [
    qlat
    qlat_cps
    qlat_grid
    gpt-lehner
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
  ;

  propagatedBuildInputs = [
  ];

  dependencies = [
    qlat
    qlat_grid
    qlat_cps
    gpt-lehner
  ];

  preConfigure = let
    gpu_extra = ''
      pwd
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
    extra = ''
      source ${qlat}/bin/cuda-mpi-qlat.sh echo
    '' + (if cudaSupport then gpu_extra else cpu_extra) + (
      if cudaSupport && lib.hasInfix "cuda" qlat-name then ''
        export num_proc=$((NIX_BUILD_CORES / 16 + 1))
      '' else ''
        export num_proc=$((NIX_BUILD_CORES / 4 + 1))
      ''
    );
  in extra + ''
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
    make update-sources SHELL=$SHELL
    echo
    make run-cps -j$num_proc SHELL=$SHELL
    #
    echo
    for i in *.p sparse-from-cps ; do
      if [ -d "$i" ] ; then
        if diff "$i"/log.check.txt "$i"/log.check.txt.new ; then
          echo "$i" passed
        else
          echo
          cat "$i"/log.full.txt
          echo
          echo "$i" failed
          echo
        fi
      fi
    done
    echo
    failed=false
    for i in *.p sparse-from-cps ; do
      if [ -d "$i" ] ; then
        if diff "$i"/log.check.txt "$i"/log.check.txt.new >/dev/null 2>&1 ; then
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
    # rm -rfv ./*.p/results
    # rm -rfv ./sparse-from-cps/results*
    rm -rf ./sparse-from-cps/results*
    rm -rf ./*.p
    #
    mkdir -p "$out/share/qlat/examples-py-cps"
    rsync -a --delete . "$out/share/qlat/examples-py-cps"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
