{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, qlat
, grid-lehner
, git
, is-pypi-src ? true
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, NVCC_ARCH ? "sm_86"
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "qlat_grid${qlat-name}";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    pname = "qlat_grid";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-ZnsO4Trkihq9fP8Y3viJj14IyFQgXlx99WcwORV2rMY=";
  };

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat-grid;

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    qlat
  ];

  nativeBuildInputs = [
    git
    grid-lehner
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    qlat
    grid-lehner
  ];

  dependencies = [
    qlat
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = let
    gpu_extra = ''
      pwd
      cp -pv "${../qcore/bin/NVCC.py}" "$PWD/NVCC.py"
      patchShebangs --build "$PWD/NVCC.py"
      #
      export NVCC_OPTIONS="-w -std=c++14 -arch=${NVCC_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing" # -D__DEBUG_VECUTILS__
      export QLAT_CXX="$PWD/NVCC.py -ccbin c++ $NVCC_OPTIONS"
      export QLAT_MPICXX="$PWD/NVCC.py -ccbin mpic++ $NVCC_OPTIONS"
      export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
      export QLAT_LDFLAGS="--NVCC-link" # --shared
      #
      export MPICXX="$QLAT_MPICXX"
      export CXX="$MPICXX"
      export CXXFLAGS="$QLAT_CXXFLAGS"
      export LDFLAGS="$QLAT_LDFLAGS"
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in ''
    export OMPI_CXX=c++
    export OMPI_CC=cc
    #
    CXX_ARR=($(grid-config --cxx))
    export CXX="''${CXX_ARR[0]}"
    export CXXFLAGS="''${CXX_ARR[@]:1} $CXXFLAGS"
    export LDFLAGS="''${CXX_ARR[@]:1} $LDFLAGS"
    #
    export
  '' + extra;

}
