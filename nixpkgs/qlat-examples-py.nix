{ stdenv
, lib
, config
, buildPythonPackage
, qlat
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

  src = qlat-src.examples-py;

in (buildPythonPackage.override { stdenv = stdenv; }) rec {

  inherit version src;

  pname = "qlat-examples-py${qlat-name}";

  pyproject = false;

  enableParallelBuilding = true;

  build-system = [
    qlat
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
      export mpi_options="$mpi_options bash bind-gpu-qlat.sh"
      #
      export q_num_mp_processes=0
    '';
    cpu_extra = ''
      if [ "$(uname)" == "Darwin" ]; then
        export q_num_mp_processes=0
      fi
    '';
    extra = (if cudaSupport then gpu_extra else cpu_extra) + (
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
    export mpi_options="--oversubscribe --bind-to none $mpi_options"
    export SHELL=${bash}/bin/bash
    #
    echo
    export
    echo
    #
    echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    echo mpi_options=$mpi_options
    echo SHELL=$SHELL
    echo NIX_BUILD_CORES=$NIX_BUILD_CORES
    echo NIX_BUILD_TOP=$NIX_BUILD_TOP
    echo
    #
    echo num_proc=$num_proc
    #
    make update-sources SHELL=$SHELL
    echo
    make run -j$num_proc SHELL=$SHELL
    #
    echo
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
        fi
      fi
    done
    echo
    failed=false
    for i in *.p ; do
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
    rm -rf ./*.p
    #
    mkdir -p "$out/share/qlat/examples-py"
    rsync -a --delete . "$out/share/qlat/examples-py"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
