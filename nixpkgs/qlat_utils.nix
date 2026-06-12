{ stdenv
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
, jax
, jaxlib
, zlib
, eigen
, git
, which
, bash
, autoAddDriverRunpath
, openmp ? null
, qlat-src
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaSupportInLibs ? false
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, nixgl ? null
, version ? "current"
}:

let

  pname = "qlat_utils${qlat-name}";

  is-linux = lib.lists.elem builtins.currentSystem lib.platforms.linux;

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;

  pyproject = true;

  src = qlat-src.qlat-utils;

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
    jax
    jaxlib
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = let
    gl_extra = if nixgl == null
    then ''
      mkdir -pv "$out/bin"
      echo "#!/usr/bin/env bash" >$out/bin/nixgl-qlat.sh
      echo '"$@"' >>$out/bin/nixgl-qlat.sh
    ''
    else ''
      mkdir -pv "$out/bin"
      ls "${nixgl}/bin/nixGL"
      echo "#!/usr/bin/env bash" >$out/bin/cuda-qlat.sh
      echo "# Source nixGL env var" >>$out/bin/cuda-qlat.sh
      echo >>$out/bin/nixgl-qlat.sh
      cat "${nixgl}/bin/nixGL" | grep -v 'exec ' | grep -v '^#!' >>$out/bin/nixgl-qlat.sh
      echo >>$out/bin/nixgl-qlat.sh
      echo '"$@"' >>$out/bin/nixgl-qlat.sh
    '';
    gpu_extra = ''
      pwd
      mkdir -pv "$out/bin"
      #
      cp -pv "${qlat-src.qcore}/bin/NVCC.py" "$out/bin/NVCC.py"
      patchShebangs --build "$out/bin/NVCC.py"
      #
      echo "#!/usr/bin/env bash" >$out/bin/cuda-qlat.sh
      echo "# Source nixGL and NVCC env var" >>$out/bin/cuda-qlat.sh
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
        export JAX_PLATFORMS="${if cudaSupportInLibs then "cuda" else "cpu"}"
        #
        ${if cudaSupport || !is-linux then "export q_num_mp_processes=0" else ""}
      EOF
      #
      echo >>$out/bin/cuda-qlat.sh
      echo '"$@"' >>$out/bin/cuda-qlat.sh
    '';
    cpu_extra = ''
      mkdir -pv "$out/bin"
      echo "#!/usr/bin/env bash" >$out/bin/cuda-qlat.sh
      cat >>"$out/bin/cuda-qlat.sh" <<EOF
        source $out/bin/nixgl-qlat.sh :
        export JAX_PLATFORMS="${if cudaSupportInLibs then "cuda" else "cpu"}"
        ${if cudaSupport || !is-linux then "export q_num_mp_processes=0" else ""}
      EOF
      #
      echo >>$out/bin/cuda-qlat.sh
      echo '"$@"' >>$out/bin/cuda-qlat.sh
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in gl_extra + extra + ''
    source $out/bin/cuda-qlat.sh echo
    #
    echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
    echo "CXX=$CXX"
    echo "CXXFLAGS=$CXXFLAGS"
    echo "LDFLAGS=$LDFLAGS"
    #
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
