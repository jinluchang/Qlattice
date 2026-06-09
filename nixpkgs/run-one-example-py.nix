{
  # Name of the Python example/test to run. Must match a .py file in examples-py/.
  # Defaults to "auto-contract-01" if null.
  # E.g. "utils", "auto-contract-01", "hmc-pions", "gf-utils"
  test-name ? null,
  # Suffix selecting a pkgs variant (constructed by mk-qlat-name in options.nix).
  # Defaults to "" (the default pkgs set with grid-gpt and cps enabled).
  # Available options:
  #   ""                            — default (grid-gpt + cps)
  #   "-std"                        — no grid-gpt, no cps
  #   "-cpsless"                    — grid-gpt, no cps
  #   "-gridless"                   — no grid-gpt, cps
  #   "-clang"                      — use clang compiler
  #   "-ucxless"                    — disable UCX
  #   "-pypi"                       — build from PyPI sources
  #   "-cuda"                       — CUDA with GPU codegen (use-cuda)
  #   "-cudasupport"                — full CUDA support (use-cudasupport)
  #   Combined variants (from options-list in options.nix):
  #   "-clang-ucxless", "-cuda-ucxless", "-gridless-cubaquadless",
  #   "-gridless-clang", "-std-ucxless", "-std-clang-ucxless",
  #   "-cpsless-ucxless", "-cpsless-clang-ucxless", "-cpsless-clang",
  #   "-std-cu", "-std-cuda", "-std-cudasupport"
  qlat-name ? null,
  # Path or fetchTarball result for nixpkgs. Auto-detected from `version` if null.
  # E.g. (builtins.fetchTarball "https://channels.nixos.org/nixos-26.05/nixexprs.tar.xz")
  nixpkgs ? null,
  # NixOS release version string. Used to fetch matching nixpkgs channel when `nixpkgs` is null.
  # E.g. "26.05", "25.11". Defaults to "" (current system nixpkgs) if null.
  version ? null,
  # Number of NVIDIA GPUs as a string. Auto-detected from /dev/nvidia* if null.
  # E.g. "0", "1", "2"
  ngpu ? null,
  # CUDA compute capability string (major.minor). Auto-detected via nvidia-smi if null.
  # E.g. "8.6", "8.0"
  cudaCapability ? null,
  # Enable CUDA forward compatibility (bool). Defaults to false if null.
  cudaForwardCompat ? null,
  # Use gitee mirrors instead of github (bool). Defaults to false if null.
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

  lib = pkgs.lib;
  mpi = pkgs.mpi;

in (pkgs.python3.pkgs.buildPythonPackage.override { stdenv = pkgs.qlat-stdenv; }) rec {

  pname = "qlat-example-py-${test-name}";
  version = pkgs.qlat-version;

  src = pkgs.qlat-src.examples-py;

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
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
    export mpi_options="--oversubscribe --bind-to none $mpi_options"
    export SHELL=${pkgs.bash}/bin/bash
    make -B ${test-name}.log SHELL=$SHELL
    echo
    if [ ! -d "${test-name}.py.p" ] ; then
      echo "${test-name}.py.p" directory not found
      false
    fi
    if diff "${test-name}.py.p"/log.check.txt "${test-name}.py.p"/log.check.txt.new ; then
      echo "${test-name}.py.p" passed
    else
      echo
      cat "${test-name}.py.p"/log.full.txt
      echo
      echo "${test-name}.py.p" failed
      false
    fi
    echo
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
