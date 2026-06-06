{ stdenv
, fetchPypi
, config
, lib
, buildPythonPackage
, cython
, meson-python
, meson
, ninja
, pkg-config
, numpy
, psutil
, zlib
, eigen
, git
, which
, bash
, autoAddDriverRunpath
, openmp ? null
, use-pypi ? null
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? null
}:

let

  version-pypi = use-pypi;

  src-pypi = builtins.fetchTarball "https://files.pythonhosted.org/packages/source/q/qlat_utils/qlat_utils-${version-pypi}.tar.gz";

  version-local = builtins.replaceStrings [ "\n" ] [ "" ] (builtins.readFile ../VERSION) + "-current";
  src-local = ../qlat-utils;

  pname = "qlat_utils${qlat-name}";
  version = if use-pypi != null then version-pypi else version-local;

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;

  pyproject = true;

  src = if use-pypi != null then src-pypi else src-local;

  enableParallelBuilding = true;

  build-system = [
    meson-python
    pkg-config
    cython
    numpy
  ];

  nativeBuildInputs = [
    git
    which
    bash
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ++ lib.optionals cudaSupport [ nixgl ]
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
    libcublas
    libcufft
  ])
  ++ lib.optionals cudaSupport [ autoAddDriverRunpath ]
  ;

  dependencies = [
    bash
    meson
    ninja
    pkg-config
    cython
    numpy
    psutil
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = let
    gl_extra = if nixgl == null
    then ''
      mkdir -pv "$out/bin"
      echo >$out/bin/nixgl-qlat.sh
      echo '"$@"' >>$out/bin/nixgl-qlat.sh
    ''
    else ''
      mkdir -pv "$out/bin"
      ls "${nixgl}/bin/nixGL"
      echo "# Source nixGL env var" >$out/bin/cuda-qlat.sh
      echo >>$out/bin/nixgl-qlat.sh
      cat "${nixgl}/bin/nixGL" | grep -v 'exec ' | grep -v '^#!' >>$out/bin/nixgl-qlat.sh
      echo >>$out/bin/nixgl-qlat.sh
      echo '"$@"' >>$out/bin/nixgl-qlat.sh
    '';
    gpu_extra = ''
      pwd
      mkdir -pv "$out/bin"
      #
      cp -pv "${../qcore/bin/NVCC.py}" "$out/bin/NVCC.py"
      patchShebangs --build "$out/bin/NVCC.py"
      #
      echo "# Source nixGL and NVCC env var" >$out/bin/cuda-qlat.sh
      echo >>$out/bin/cuda-qlat.sh
      #
      cat >>"$out/bin/cuda-qlat.sh" <<EOF
        #
        source $out/bin/nixgl-qlat.sh :
        #
        # GXX="-w"
        #
        GXX=""
        GXX+=" -Xcudafe '--diag_suppress=20014'"
        GXX+=" -Xcudafe '--diag_suppress=20236'"
        GXX+=" -Xcudafe '--diag_suppress=20012'"
        GXX+=" -Xcudafe '--diag_suppress=20011'"
        GXX+=" -Xcudafe '--diag_suppress=1160'"
        GXX+=" -Xcudafe '--diag_suppress=177'"
        GXX+=" -Xcudafe '--diag_suppress=550'"
        #
        export NVCC_OPTIONS="-std=c++17 -arch=${nvcc-arch} --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing \$GXX" # -D__DEBUG_VECUTILS__
        export QLAT_CXX="$out/bin/NVCC.py -ccbin c++ \$NVCC_OPTIONS"
        export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
        export QLAT_LDFLAGS="--NVCC-link" # --shared
        #
        export CXX="\$QLAT_CXX"
        export CXXFLAGS="\$QLAT_CXXFLAGS"
        export LDFLAGS="\$QLAT_LDFLAGS"
        #
      EOF
      #
      echo >>$out/bin/cuda-qlat.sh
      echo '"$@"' >>$out/bin/cuda-qlat.sh
      #
      source $out/bin/cuda-qlat.sh echo
      #
      echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
      #
      echo "CXX=$CXX"
      echo "CXXFLAGS=$CXXFLAGS"
      echo "LDFLAGS=$LDFLAGS"
      #
      echo
    '';
    cpu_extra = ''
      mkdir -pv "$out/bin"
      cat >"$out/bin/cuda-qlat.sh" <<EOF
        source $out/bin/nixgl-qlat.sh :
      EOF
      echo '"$@"' >>$out/bin/cuda-qlat.sh
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in gl_extra + extra + ''
    echo
    export
    echo
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
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
