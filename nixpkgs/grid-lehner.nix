{ stdenv
, config
, lib
, fetchurl
, fetchFromGitHub
, mpi
, c-lime
, zlib
, openssl
, gsl
, hdf5-cpp
, gmp
, mpfr
, fftw
, fftwFloat
, git
, autoconf
, automake
, which
, autoAddDriverRunpath
, openmp ? null
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
, nvcc-arch ? "sm_86"
, cpuinfo ? "flags : avx2"
}:

let

  grid-stdenv = if cudaSupport then cudaPackages.backendStdenv else stdenv;

  cpuinfo-has-avx2 = (builtins.match ".*flags.*:.*avx2.*" cpuinfo) != null;

in grid-stdenv.mkDerivation rec {

  pname = "Grid-lehner";
  version = "e138258129faed5cb12fc0e5bbf829f69c6a2a3d";

  src = fetchFromGitHub {
    owner = "lehner";
    repo = "Grid";
    rev = version;
    hash = "sha256-m1wpTWgoRFeZJbidy7ckD1iobFnaNOUjEnhYuJkxe4s=";
  };

  enableParallelBuilding = true;

  nativeBuildInputs = [
    mpi
    git
    autoconf
    automake
    which
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    mpi
    c-lime
    zlib
    fftw
    fftwFloat
    gsl
    openssl
    hdf5-cpp
    gmp
    mpfr
  ]
  ++ lib.optional grid-stdenv.cc.isClang openmp
  ++ lib.optionals cudaSupport (with cudaPackages; [
    cuda_cccl
    cuda_cudart
    cuda_profiler_api
    libcublas
  ])
  ++ lib.optionals cudaSupport [ autoAddDriverRunpath ]
  ;

  preConfigure = let
    eigen-file-name = "eigen-3.4.0.tar.bz2";
    eigen-src = fetchurl {
      url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2";
      hash = "sha256-tMGYRg66byjTSJTjpXEJmIGFFRBNbnTlzDMc4x5G5iY=";
    };
    cpu_cxx = "c++";
    gpu_cxx = "nvcc";
    cpu_cflags = "-fPIC -w -Wno-psabi";
    gpu_cflags = "-Xcompiler -fPIC -ccbin mpic++ -arch=${nvcc-arch} -w";
    cpu_ldflags = "";
    gpu_ldflags = "-Xcompiler -fopenmp -ccbin mpic++";
    cxx = if cudaSupport then gpu_cxx else cpu_cxx;
    cflags = if cudaSupport then gpu_cflags else cpu_cflags;
    ldflags = if cudaSupport then gpu_ldflags else cpu_ldflags;
    cpu_extra = ''
    '';
    gpu_extra = ''
      which nvcc
      nvcc --version
    '';
    extra = if cudaSupport then gpu_extra else cpu_extra;
  in ''
    echo "-- deploying Eigen source..."
    cp -pv '${eigen-src}' '${eigen-file-name}'
    bash ./scripts/update_eigen.sh '${eigen-file-name}'
    rm '${eigen-file-name}'
    # patch Eigen/unsupported/Eigen/CXX11/Tensor scripts/eigen-3.3.5.Tensor.patch
    #
    echo '-- generating Make.inc files...'
    bash ./scripts/filelist
    echo '-- generating configure script...'
    autoreconf -fvi
    #
    echo '-- set FLAGS ...'
    export CXX=${cxx}
    export CFLAGS="${cflags}"
    export CXXFLAGS="${cflags}"
    export LDFLAGS="${ldflags}"
    echo CXX="$CXX"
    echo CFLAGS="$CFLAGS"
    echo CXXFLAGS="$CXXFLAGS"
    echo LDFLAGS="$LDFLAGS"
    #
    export OMPI_CXX=c++
    export OMPI_CC=cc
    #
    which mpic++
    mpic++ --version
    #
    which c++
    c++ --version
    #
    echo
    echo 'grid-stdenv=${grid-stdenv.cc}'
    echo
  '' + extra;

  configureFlags = let
    cpu_avx2_flags = [
      "--enable-simd=AVX2"
      "--enable-alloc-align=4k"
      "--enable-comms=mpi-auto"
      "--enable-gparity=no"
    ];
    cpu_gen16_flags = [
      "--enable-simd=GEN"
      "--enable-gen-simd-width=16"
      "--enable-alloc-align=4k"
      "--enable-comms=mpi-auto"
      "--enable-gparity=no"
    ];
    cpu_flags = if cpuinfo-has-avx2 then cpu_avx2_flags else cpu_gen16_flags;
    gpu_flags = [
      "--enable-simd=GPU"
      "--enable-gen-simd-width=32"
      "--enable-alloc-align=4k"
      "--enable-comms=mpi-auto"
      "--enable-gparity=no"
      "--disable-fermion-reps"
      "--enable-unified=no"
      "--enable-accelerator=cuda"
      "--enable-accelerator-cshift"
    ];
    flags = if cudaSupport then gpu_flags else cpu_flags;
  in flags;

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
