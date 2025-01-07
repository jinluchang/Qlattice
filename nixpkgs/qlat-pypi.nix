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
, nvcc-arch ? "sm_86"
, nixgl ? ""
}:

let
  orig-stdenv = stdenv;
in buildPythonPackage rec {

  pname = "qlat-pypi${qlat-name}";
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

  preConfigure = let
    gpu_extra = ''
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
      echo
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
