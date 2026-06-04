{
  nixpkgs ? null,
  version ? null,
  testName ? "auto-contract-01",
  cudaSupport ? false,
  cudaCapability ? null,
  cudaForwardCompat ? null,
  nvcc-arch ? "sm_86",
  ngpu ? null,
  use-gitee ? null,
  use-pypi ? null,
}:

let
  q-pkgs-module = import ./q-pkgs.nix {
    inherit nixpkgs version cudaCapability cudaForwardCompat use-gitee;
    ngpu = if ngpu != null then ngpu else if cudaSupport then "1" else null;
  };

  use-gitee-wd = if use-gitee == null then false else use-gitee;

  opts = {
    use-cuda = cudaSupport;
    use-cudasupport = false;
    use-cuda-software = cudaSupport;
    use-pypi = use-pypi;
  } // (if cudaSupport then {
    use-cuda-software = true;
    use-clang = false;
  } else {});

  q-pkgs-result = q-pkgs-module.mk-q-pkgs opts;

  pkgs-key = builtins.head (builtins.filter (k: builtins.substring 0 4 k == "pkgs") (builtins.attrNames q-pkgs-result));
  pkgs = q-pkgs-result.${pkgs-key};
  qlat = pkgs.qlat;

  version-pypi = use-pypi;
  qlat-src-pypi = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/Qlattice" else "https://github.com/jinluchang/Qlattice";
    ref = "refs/tags/v${version-pypi}";
  };

  version-val = if use-pypi != null then version-pypi else builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  src = if use-pypi != null then "${qlat-src-pypi}/examples-py" else ../examples-py;

  lib = pkgs.lib;

  mpi = pkgs.mpi;

  cudaSupport-val = cudaSupport;

in (pkgs.python3.pkgs.buildPythonPackage.override { stdenv = pkgs.qlat-stdenv; }) rec {

  inherit src;

  version = version-val;

  pname = "qlat-example-${testName}";

  pyproject = false;
  enableParallelBuilding = true;

  build-system = [ qlat ];

  nativeBuildInputs = [
    pkgs.git
    pkgs.time
    mpi
    pkgs.mpiCheckPhaseHook
    pkgs.openssh
    pkgs.which
    pkgs.rsync
  ] ++ lib.optionals cudaSupport-val (with pkgs.cudaPackages; [ cuda_nvcc ]);

  dependencies = [ qlat ];

  preConfigure = let
    gpu_extra = ''
      pwd
      source ${qlat}/bin/cuda-mpi-qlat.sh echo
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
      echo "CXX=$CXX"
      echo "MPICXX=$MPICXX"
      export NGPU=${if ngpu != null then ngpu else "1"}
      export mpi_options="$mpi_options bind-gpu-qlat.sh"
      export q_num_mp_processes=0
      export num_proc=$((NIX_BUILD_CORES / 16 + 1))
    '';
    cpu_extra = ''
      if [ "$(uname)" == "Darwin" ]; then
        export q_num_mp_processes=0
      fi
      export num_proc=$((NIX_BUILD_CORES / 4 + 1))
    '';
    extra = if cudaSupport-val then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
    export mpi_options="--oversubscribe --bind-to none $mpi_options"
    export SHELL=${pkgs.bash}/bin/bash
    echo num_proc=$num_proc
    make clean SHELL=$SHELL
    make -B ${testName}.log SHELL=$SHELL
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
    mkdir -p "$out/share/qlat/examples-py"
    rsync -a --delete . "$out/share/qlat/examples-py"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version-val} >"$out"/share/version/${pname}
  '';
}
