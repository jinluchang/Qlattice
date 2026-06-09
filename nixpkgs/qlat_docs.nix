{ stdenv
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
, qlat-src
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, version ? "current"
}:

let

  src = "${qlat-src}/docs/";

in buildPythonPackage.override { stdenv = stdenv; } rec {

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
    examples-cpp = "${qlat-src}/examples-cpp/";
    examples-cpp-grid = "${qlat-src}/examples-cpp-grid/";
    examples-py = "${qlat-src}/examples-py/";
    examples-py-cps = "${qlat-src}/examples-py-cps/";
    examples-py-gpt = "${qlat-src}/examples-py-gpt/";
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
      export mpi_options="$mpi_options bash bind-gpu-qlat.sh"
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
