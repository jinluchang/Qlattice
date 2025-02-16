{ fetchPypi
, stdenv
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
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? ""
, ngpu ? "1"
}:

let
  orig-stdenv = stdenv;
in buildPythonPackage rec {

  pname = "qlat-examples-cpp${qlat-name}";
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
    qlat
  ];

  preConfigure = let
    gpu_extra = ''
      pwd
      cp -pv "${../qcore/bin/NVCC.py}" "$PWD/NVCC.py"
      patchShebangs --build "$PWD/NVCC.py"
      #
      cp -pv "${../qcore/bin/bind-gpu-qlat.sh}" "$PWD/bind-gpu-qlat.sh"
      patchShebangs --build "$PWD/bind-gpu-qlat.sh"
      export NGPU=${ngpu}
      export mpi_options="$mpi_options $PWD/bind-gpu-qlat.sh"
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
      if [ -d /run/opengl-driver/lib ] ; then
        export LD_LIBRARY_PATH=/run/opengl-driver/lib"''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      fi
      echo $LD_LIBRARY_PATH
      echo
      export q_num_mp_processes=0
      export num_proc=$((NIX_BUILD_CORES / 16 + 1))
    '';
    cpu_extra = ''
      export num_proc=$((NIX_BUILD_CORES / 4 + 1))
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
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
    echo num_proc=$num_proc
    #
    make run -j$num_proc SHELL=$SHELL
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
    rm -rfv ./*/build/results
    #
    pwd
    mkdir -p "$out/share/qlat/examples-cpp"
    rsync -a --delete . "$out/share/qlat/examples-cpp"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
