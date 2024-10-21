{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
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
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, NVCC_ARCH ? "sm_86"
, nixgl ? ""
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat-pypi-pkgs${qlat-name}";
  version = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";

  pyproject = false;

  srcs = [
    ../qlat-utils
    ../qlat
    ../qlat-cps
    ../qlat-grid
  ];
  sourceRoot = ".";

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

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
  ++ lib.optionals cudaSupport [ nixgl ]
  ;

  propagatedBuildInputs = [
  ];

  dependencies = [
    qlat_utils
    qlat
    qlat_cps
    qlat_grid
  ];

  preConfigure = ''
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
    mkdir qlat-pypi-pkgs
    python3 -m build -ns -o qlat-pypi-pkgs ./qlat-utils
    python3 -m build -ns -o qlat-pypi-pkgs ./qlat
    python3 -m build -ns -o qlat-pypi-pkgs ./qlat-cps
    python3 -m build -ns -o qlat-pypi-pkgs ./qlat-grid
    #
    mkdir -p "$out/share/qlat-pypi-pkgs"
    rsync -a --delete qlat-pypi-pkgs "$out/share/qlat-pypi-pkgs"
  '';

  dontBuild = true;
  dontInstall = true;

}
