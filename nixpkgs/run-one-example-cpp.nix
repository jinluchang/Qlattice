{
  # Name of the C++ example/test to run. Must match a subdirectory in examples-cpp/.
  # Defaults to "simple-1" if null.
  # E.g. "simple-1", "hmc", "flowed-hmc", "benchmark", "field-rng-tests"
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

  test-name = if test-name-ini == null then "simple-1" else test-name-ini;

  qlat-name = if qlat-name-ini == null then "" else qlat-name-ini;

  inherit (q-pkgs) nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;

  pkgs = q-pkgs."pkgs${qlat-name}";
  qlat = pkgs.qlat;
  opts = pkgs.qlat-options;
  cudaSupport = opts.use-cuda-software;

  lib = pkgs.lib;
  mpi = pkgs.mpi;

in (pkgs.python3.pkgs.buildPythonPackage.override { stdenv = pkgs.qlat-stdenv; }) rec {

  pname = "qlat-example-cpp-${test-name}";
  version = pkgs.qlat-version;

  src = pkgs.qlat-src.examples-cpp;

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
    pkgs.meson
    pkgs.ninja
    pkgs.pkg-config
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
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
    export SHELL=${pkgs.bash}/bin/bash
    export LD_LIBRARY_PATH="$(python3 -m qlat qlat-config --LD_LIBRARY_PATH)"
    make -f ../Makefile -C ${test-name} run-proj SHELL=$SHELL
    echo
    if [ ! -f "${test-name}/build/log.check.txt.new" ] ; then
      echo "${test-name}/build/log.check.txt.new" not found
      false
    fi
    if diff "${test-name}/build/log.check.txt" "${test-name}/build/log.check.txt.new" ; then
      echo "${test-name}" passed
    else
      echo
      cat "${test-name}/build/log.full"
      echo
      echo "${test-name}" failed
      false
    fi
    echo
    mkdir -p "$out/share/qlat/examples-cpp"
    cp -r ${test-name} "$out/share/qlat/examples-cpp/" 2>/dev/null || true
    cp ${test-name}/log "$out/share/qlat/examples-cpp/${test-name}.log" 2>/dev/null || true
    cp ${test-name}/log.full "$out/share/qlat/examples-cpp/${test-name}.log.full" 2>/dev/null || true
  '';

  configurePhase = "runHook preConfigure";
  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
