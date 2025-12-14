{ fetchPypi
, stdenv
, lib
, config
, buildPythonPackage
, mpi4py
, sympy
, scipy
, jax
, jaxlib
, qlat_utils
, mpi
, git
, which
, fftw
, fftwFloat
, gsl
, cuba-quad
, use-pypi ? null
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? ""
}:

let

  version-pypi = use-pypi;
  src-pypi = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat/qlat-${version-pypi}.tar.gz";

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat;

  pname = "qlat${qlat-name}";
  version = if use-pypi != null then version-pypi else version-local;

in buildPythonPackage {

  pname = pname;
  version = version;

  pyproject = true;

  src = if use-pypi != null then src-pypi else src-local;

  enableParallelBuilding = true;

  stdenv = stdenv;

  build-system = [
    qlat_utils
  ];

  nativeBuildInputs = [
    git
    mpi
    which
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ++ lib.optionals cudaSupport [ nixgl ]
  ;

  propagatedBuildInputs = [
    mpi
    fftw
    fftwFloat
    gsl
    cuba-quad
    qlat_utils
  ];

  dependencies = [
    qlat_utils
    mpi4py
    sympy
    scipy
    jax
    jaxlib
  ];

  # requiredSystemFeatures = [ "require-cuda" ];

  postPatch = ''
      sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = let
    gpu_extra = ''
      pwd
      cp -pv "${../qcore/bin/NVCC.py}" "$PWD/NVCC.py"
      patchShebangs --build "$PWD/NVCC.py"
      #
      GXX=""
      GXX+=" -Xcudafe '--diag_suppress=20014'"
      GXX+=" -Xcudafe '--diag_suppress=20236'"
      GXX+=" -Xcudafe '--diag_suppress=20012'"
      GXX+=" -Xcudafe '--diag_suppress=20011'"
      GXX+=" -Xcudafe '--diag_suppress=177'"
      GXX+=" -Xcudafe '--diag_suppress=550'"
      # GXX="-w"
      #
      export NVCC_OPTIONS="-std=c++17 -arch=${nvcc-arch} --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing $GXX" # -D__DEBUG_VECUTILS__
      export QLAT_CXX="$PWD/NVCC.py -ccbin c++ $NVCC_OPTIONS"
      export QLAT_MPICXX="$PWD/NVCC.py -ccbin ${mpi.dev}/bin/mpic++ $NVCC_OPTIONS"
      export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
      export QLAT_LDFLAGS="--NVCC-link" # --shared
      #
      export OMPI_CXX=c++
      export OMPI_CC=cc
      #
      export MPICXX="$QLAT_MPICXX"
      export CXX="$QLAT_MPICXX"
      export CXXFLAGS="$QLAT_CXXFLAGS"
      export LDFLAGS="$QLAT_LDFLAGS"
      #
      which nvcc
      echo
      nvcc --version
      echo
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
  in ''
  '' + extra + ''
    # export
  '';

  preFixup = ''
    # echo
    # echo ldd $out
    # ldd $out/lib/python3*/site-packages/cqlat.cpython-*.so
    # echo
    # echo readelf -d $out
    # readelf -d $out/lib/python3*/site-packages/cqlat.cpython-*.so
    # echo
  '';

  postFixup = ''
    # echo
    # echo ldd $out
    # ldd $out/lib/python3*/site-packages/cqlat.cpython-*.so
    # echo
    # echo readelf -d $out
    # readelf -d $out/lib/python3*/site-packages/cqlat.cpython-*.so
    # echo
    patchShebangs --build "$out"/bin/bind-gpu-qlat.sh
    #
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
