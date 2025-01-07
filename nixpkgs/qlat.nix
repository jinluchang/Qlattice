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
, cuba
, is-pypi-src ? true
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? ""
}:

let

  orig-stdenv = stdenv;

  version-pypi = "0.76";
  src-pypi = fetchPypi {
    pname = "qlat";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-5ZTZOGN6M8O/UdRVqO/cVhecuQmHBexTLtbPrSOb+qE=";
  };

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat;

in

  buildPythonPackage rec {

    pname = "qlat${qlat-name}";
    version = if is-pypi-src then version-pypi else version-local;

    pyproject = true;

    src = if is-pypi-src then src-pypi else src-local;

    enableParallelBuilding = true;

    stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

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
      cuba
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
        export QLAT_MPICXX="$PWD/NVCC.py -ccbin mpic++ $NVCC_OPTIONS"
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
    '';

  }
