{ fetchPypi
, stdenv
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
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? null
, ngpu ? "1"
, use-pypi ? null
}:

let

  orig-stdenv = stdenv;

  version-pypi = use-pypi;
  qlat-src-pypi = builtins.fetchGit {
    url = "https://github.com/jinluchang/Qlattice";
    ref = "refs/tags/v${version-pypi}";
  };

  version = if use-pypi != null then version-pypi else builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  src = if use-pypi != null then "${qlat-src-pypi}/examples-py-cps" else ../examples-py-cps;

in buildPythonPackage rec {

  inherit version src;

  pname = "qlat-examples-py-cps${qlat-name}";

  pyproject = false;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

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
  ++ lib.optionals (nixgl != null) [ nixgl ]
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
      cp -pv "${../qcore/bin/NVCC.py}" "$PWD/NVCC.py"
      patchShebangs --build "$PWD/NVCC.py"
      #
      export NGPU=${ngpu}
      export mpi_options="$mpi_options bind-gpu-qlat.sh"
      #
      GXX=""
      GXX+=" -Xcudafe '--diag_suppress=20014'"
      GXX+=" -Xcudafe '--diag_suppress=20236'"
      GXX+=" -Xcudafe '--diag_suppress=20012'"
      GXX+=" -Xcudafe '--diag_suppress=20011'"
      GXX+=" -Xcudafe '--diag_suppress=177'"
      GXX+=" -Xcudafe '--diag_suppress=550'"
      # GXX="-w"
      #
      export NVCC_OPTIONS="-std=c++17 -arch=${nvcc-arch} --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing $GXX" # -D__DEBUG_VECUTILS__
      export QLAT_CXX="$PWD/NVCC.py -ccbin c++ $NVCC_OPTIONS"
      export QLAT_MPICXX="$PWD/NVCC.py -ccbin mpic++ $NVCC_OPTIONS"
      export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
      export QLAT_LDFLAGS="--NVCC-link" # --shared
      #
      export OMPI_CXX=c++
      export OMPI_CC=cc
      #
      export MPICXX="$QLAT_MPICXX"
      export CXX="$QLAT_MPICXX"
      export CXXFLAGS="$QLAT_CXXFLAGS"
      export LDFLAGS="$QLAT_LDFLAGS"
      #
      export q_num_mp_processes=0
      export num_proc=$((NIX_BUILD_CORES / 16 + 1))
    '';
    cpu_extra = ''
      if [ "$(uname)" == "Darwin" ]; then
        export q_num_mp_processes=0
      fi
      export num_proc=$((NIX_BUILD_CORES / 4 + 1))
    '';
    nixgl_extra = if nixgl == null then "" else ''
      which nixGL
      echo
      echo "run with nixGL"
      echo
      nixGL qlat-utils-config
      echo
      cat $(which nixGL) | grep -v 'exec ' | grep -v '^#!' > nix-gl.sh
      echo
      echo cat nix-gl.sh
      cat nix-gl.sh
      source nix-gl.sh
      echo
      echo $LD_LIBRARY_PATH
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + nixgl_extra + ''
    export OMP_NUM_THREADS=2
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
    rm -rfv ./*.p/results
    rm -rfv ./sparse-from-cps/results*
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
