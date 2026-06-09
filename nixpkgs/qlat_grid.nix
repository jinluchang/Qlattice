{ stdenv
, lib
, config
, buildPythonPackage
, meson-python
, qlat
, grid-lehner
, git
, mpi
, which
, qlat-src
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, version ? "current"
}:

let

  pname = "qlat_grid${qlat-name}";

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;

  pyproject = true;

  src = qlat-src.qlat-grid;

  enableParallelBuilding = true;

  build-system = [
    qlat
    meson-python
  ];

  nativeBuildInputs = [
    git
    grid-lehner
    which
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    qlat
    grid-lehner
    mpi
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
      echo
    '';
    cpu_extra = ''
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in ''
    # CXX_ARR=($(grid-config --cxx))
    # export CXX="''${CXX_ARR[0]}"
    # export CXXFLAGS="''${CXX_ARR[@]:1} $CXXFLAGS"
    # export LDFLAGS="''${CXX_ARR[@]:1} $LDFLAGS"
  '' + extra + ''
    echo
    export
    echo
  '';

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
