{ stdenv
, lib
, config
, buildPythonPackage
, meson-python
, qlat
, cps
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

  pname = "qlat_cps";

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;
  # The variant suffix must not live in pname: pythonMetadataCheckPhase
  # resolves the installed dist-info by pname, and only the project name
  # matches it.
  name = "qlat_cps${qlat-name}-${version}";

  pyproject = true;

  src = qlat-src.qlat-cps;

  enableParallelBuilding = true;

  build-system = [
    qlat
    meson-python
  ];

  nativeBuildInputs = [
    git
    which
  ]
  ++ lib.optionals cudaSupport [ cudaPackages.cuda_nvcc ]
  ;

  propagatedBuildInputs = [
    qlat
    cps
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
    extra = ''
      source ${qlat}/bin/cuda-mpi-qlat.sh echo
    '' + (if cudaSupport then gpu_extra else cpu_extra);
  in ''
  '' + extra + ''
    echo
    export
    echo
  '';

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}${qlat-name}
  '';

}
