{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, meson-python
, qlat
, cps
, git
, mpi
, which
, use-pypi ? null
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
}:

let

  version-pypi = use-pypi;
  src-pypi = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat_cps/qlat_cps-${version-pypi}.tar.gz";

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat-cps;

  pname = "qlat_cps${qlat-name}";
  version = if use-pypi != null then version-pypi else version-local;

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;

  pyproject = true;

  src = if use-pypi != null then src-pypi else src-local;

  enableParallelBuilding = true;

  build-system = [
    qlat
    meson-python
  ];

  nativeBuildInputs = [
    git
    which
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
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
