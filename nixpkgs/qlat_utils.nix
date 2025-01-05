{ stdenv
, fetchPypi
, config
, lib
, buildPythonPackage
, cython
, meson-python
, pkg-config
, numpy
, psutil
, zlib
, eigen
, git
, which
, autoAddDriverRunpath
, openmp ? null
, is-pypi-src ? true
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
}:

let

  orig-stdenv = stdenv;

  version-pypi = "0.75";
  src-pypi = fetchPypi {
    pname = "qlat_utils";
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-Q1/xAOMpaPipVt+ZqKLdkJmWrDQpui0z7bPWgRKWQkg=";
  };

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat-utils;

in

  buildPythonPackage rec {

    pname = "qlat_utils${qlat-name}";
    version = if is-pypi-src then version-pypi else version-local;

    pyproject = true;

    src = if is-pypi-src then src-pypi else src-local;

    enableParallelBuilding = true;

    stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

    build-system = [
      meson-python
      pkg-config
      cython
      numpy
    ];

    nativeBuildInputs = [
      git
      which
    ]
    ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
    ;

    propagatedBuildInputs = [
      zlib
      eigen
    ]
    ++ lib.optional stdenv.cc.isClang openmp
    ++ lib.optionals cudaSupport (with cudaPackages; [
      cuda_cccl
      cuda_cudart
      cuda_profiler_api
      libcufft
    ])
    ++ lib.optionals cudaSupport [ autoAddDriverRunpath ]
    ;

    dependencies = [
      meson-python
      pkg-config
      cython
      numpy
      psutil
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
        GXX=""
        GXX+=" -Xcudafe '--diag_suppress=20014'"
        GXX="-w"
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
        export CXX="$QLAT_CXX"
        export CXXFLAGS="$QLAT_CXXFLAGS"
        export LDFLAGS="$QLAT_LDFLAGS"
        #
        echo $LD_LIBRARY_PATH
        echo
      '';
      cpu_extra = ''
      '';
      extra = if cudaSupport then gpu_extra else cpu_extra;
    in extra + ''
      # export
    '';

    preFixup = ''
      # echo
      # echo ldd $out
      # ldd $out/lib/python3*/site-packages/qlat_utils/timer.cpython-*.so
      # echo
      # echo readelf -d $out
      # readelf -d $out/lib/python3*/site-packages/qlat_utils/timer.cpython-*.so
      # echo
    '';

    postFixup = ''
      # echo
      # echo ldd $out
      # ldd $out/lib/python3*/site-packages/qlat_utils/timer.cpython-*.so
      # echo
      # echo readelf -d $out
      # readelf -d $out/lib/python3*/site-packages/qlat_utils/timer.cpython-*.so
      # echo
    '';

  }
