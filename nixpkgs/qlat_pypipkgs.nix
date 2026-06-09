{ stdenv
, lib
, config
, buildPythonPackage
, meson-python
, qlat
, qlat_utils
, qlat_cps
, qlat_grid
, build
, wheel
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

  srcs = [
    (builtins.path { name = "qlat-utils"; path = qlat-src.qlat-utils; })
    (builtins.path { name = "qlat"; path = qlat-src.qlat; })
    (builtins.path { name = "qlat-cps"; path = qlat-src.qlat-cps; })
    (builtins.path { name = "qlat-grid"; path = qlat-src.qlat-grid; })
  ];

in (buildPythonPackage.override { stdenv = stdenv; }) rec {

  inherit version srcs;

  pname = "qlat-pypi${qlat-name}";

  pyproject = false;

  sourceRoot = ".";

  enableParallelBuilding = true;

  build-system = [
    qlat
    qlat_cps
    qlat_grid
    build
    wheel
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
    qlat_utils
    qlat
    qlat_cps
    qlat_grid
    meson-python
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
      export mpi_options="$mpi_options bash bind-gpu-qlat.sh"
      #
      export num_proc=$((NIX_BUILD_CORES / 16 + 1))
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in extra + ''
    export OMP_NUM_THREADS=2
    #
    export
    echo
    ls -l
    echo
    pwd
    #
    git init
    git config user.email "ljin.luchang@gmail.com"
    git config user.name "Luchang Jin"
    git add .
    git commit -m "everything"
    mkdir qlat-pypi
    python3 -m build -ns -o qlat-pypi ./qlat-utils
    python3 -m build -ns -o qlat-pypi ./qlat
    python3 -m build -ns -o qlat-pypi ./qlat-cps
    python3 -m build -ns -o qlat-pypi ./qlat-grid
    #
    mkdir -p "$out/share/qlat-pypi"
    rsync -a --delete qlat-pypi "$out/share/"
  '';

  dontBuild = true;
  dontInstall = true;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
