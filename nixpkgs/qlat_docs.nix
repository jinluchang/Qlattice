{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
, qlat_cps
, qlat_grid
, gpt-lehner
, sphinx
, linkify-it-py
, myst-parser
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
, use-pypi ? null
, use-gitee ? null
}:

let

  use-gitee-wd = if use-gitee == null then false else use-gitee;

  version-pypi = use-pypi;
  qlat-src-pypi = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/Qlattice" else "https://github.com/jinluchang/Qlattice";
    ref = "refs/tags/v${version-pypi}";
  };

  version = if use-pypi != null then version-pypi else builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  src = if use-pypi != null then "${qlat-src-pypi}/docs" else ../docs;

in buildPythonPackage.override { stdenv = if cudaSupport then cudaPackages.backendStdenv else stdenv; } rec {

  inherit version src;

  pname = "qlat-docs${qlat-name}";

  pyproject = false;

  enableParallelBuilding = true;

  build-system = [
    qlat
    qlat_cps
    qlat_grid
    gpt-lehner
    sphinx
    linkify-it-py
    myst-parser
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
  ];

  preConfigure = let
    examples-cpp = ../examples-cpp;
    examples-cpp-grid = ../examples-cpp-grid;
    examples-py = ../examples-py;
    examples-py-cps = ../examples-py-cps;
    examples-py-gpt = ../examples-py-gpt;
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
      export mpi_options="$mpi_options bind-gpu-qlat.sh"
      #
      export q_num_mp_processes=0
      export num_proc=$((NIX_BUILD_CORES / 16 + 1))
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    pwd
    echo
    rsync -a "${examples-cpp}/" ../examples-cpp/
    rsync -a "${examples-cpp-grid}/" ../examples-cpp-grid/
    rsync -a "${examples-py}/" ../examples-py/
    rsync -a "${examples-py-cps}/" ../examples-py-cps/
    rsync -a "${examples-py-gpt}/" ../examples-py-gpt/
    #
    ls -ld *
    ls -l ..
    ls -l ../examples-cpp/
    ls -l ../examples-py/
    #
    make html
    #
    mkdir -p "$out/share/doc/qlat"
    rsync -a --delete ./build/ "$out/share/doc/qlat/"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
