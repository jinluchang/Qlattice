{ stdenv
, lib
, config
, buildPythonPackage
, meson-python
, mpi4py
, sympy
, scipy
, qlat_utils
, mpi
, git
, which
, fftw
, fftwFloat
, gsl
, cuba-quad
, qlat-src
, qlat-name ? ""
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, version ? "current"
}:

let

  pname = "qlat${qlat-name}";

in buildPythonPackage.override { stdenv = stdenv; } {

  pname = pname;
  version = version;

  pyproject = true;

  src = qlat-src.qlat;

  enableParallelBuilding = true;

  build-system = [
    qlat_utils
    meson-python
  ];

  nativeBuildInputs = [
    git
    mpi
    which
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
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
  ];

  # requiredSystemFeatures = [ "require-cuda" ];

  postPatch = ''
      sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

  preConfigure = let
    gpu_extra = ''
      pwd
      mkdir -pv "$out/bin"
      #
      echo "#!/usr/bin/env bash" >$out/bin/cuda-mpi-qlat.sh
      echo "# Source nixGL and NVCC + MPI env var" >>$out/bin/cuda-mpi-qlat.sh
      echo >>$out/bin/cuda-mpi-qlat.sh
      #
      cat >>"$out/bin/cuda-mpi-qlat.sh" <<EOF
        #
        source ${qlat_utils}/bin/cuda-qlat.sh :
        #
        export QLAT_MPICXX="${qlat_utils}/bin/NVCC.py -ccbin ${mpi.dev}/bin/mpic++ \$NVCC_OPTIONS"
        export CXX="\$QLAT_MPICXX"
        #
        export MPICXX="\$QLAT_MPICXX"
        export OMPI_CXX=c++
        #
      EOF
      #
      echo >>$out/bin/cuda-mpi-qlat.sh
      echo '"$@"' >>$out/bin/cuda-mpi-qlat.sh
      #
      source $out/bin/cuda-mpi-qlat.sh echo
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
      mkdir -pv "$out/bin"
      echo "#!/usr/bin/env bash" >$out/bin/cuda-mpi-qlat.sh
      cat >>"$out/bin/cuda-mpi-qlat.sh" <<EOF
        source ${qlat_utils}/bin/cuda-qlat.sh :
      EOF
      echo '"$@"' >>$out/bin/cuda-mpi-qlat.sh
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in ''
  '' + extra + ''
    echo
    export
    echo
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
    # patchShebangs --build "$out"/bin/bind-gpu-qlat.sh
    #
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
