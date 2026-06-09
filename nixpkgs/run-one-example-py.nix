{
  test-name ? null,
  qlat-name ? null,
  nixpkgs ? null,
  version ? null,
  ngpu ? null,
  cudaCapability ? null,
  cudaForwardCompat ? null,
  use-gitee ? null,
}:

let

  test-name-ini = test-name;

  qlat-name-ini = qlat-name;

  q-pkgs = import ./q-pkgs.nix {
    inherit nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;
  };

in let

  test-name = if test-name-ini == null then "auto-contract-01" else test-name-ini;

  qlat-name = if qlat-name-ini == null then "" else qlat-name-ini;

  inherit (q-pkgs) nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;

  pkgs = q-pkgs."pkgs${qlat-name}";
  qlat = pkgs.qlat;
  opts = pkgs.qlat-options;
  cudaSupport = opts.use-cuda-software;

  src = ../examples-py;

  lib = pkgs.lib;
  mpi = pkgs.mpi;

in (pkgs.python3.pkgs.buildPythonPackage.override { stdenv = pkgs.qlat-stdenv; }) rec {

  inherit src;

  version = qlat-version;

  pname = "qlat-example-py-${test-name}";

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
  ] ++ lib.optionals cudaSupport (with pkgs.cudaPackages; [ cuda_nvcc ]);

  dependencies = [ qlat ];

  preConfigure = let
    gpu_extra = ''
      pwd
      source ${qlat}/bin/cuda-mpi-qlat.sh echo
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
      echo "CXX=$CXX"
      echo "MPICXX=$MPICXX"
      export NGPU=${ngpu}
      export mpi_options="$mpi_options bash bind-gpu-qlat.sh"
      export q_num_mp_processes=0
    '';
    cpu_extra = ''
      if [ "$(uname)" == "Darwin" ]; then
        export q_num_mp_processes=0
      fi
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
    export mpi_options="--oversubscribe --bind-to none $mpi_options"
    export SHELL=${pkgs.bash}/bin/bash
    make -B ${test-name}.log SHELL=$SHELL
    mkdir -p "$out/share/qlat/examples-py"
    cp ${test-name}.py "$out/share/qlat/examples-py/"
    cp ${test-name}.log.json "$out/share/qlat/examples-py/" 2>/dev/null || true
    cp ${test-name}.log "$out/share/qlat/examples-py/"
    cp -r ${test-name}.py.p "$out/share/qlat/examples-py/" 2>/dev/null || true
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
