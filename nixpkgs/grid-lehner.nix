{ stdenv
, config
, lib
, fetchurl
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
, use-gitee ? null
}:

let
  cpuinfo-has-avx2 = (builtins.match ".*flags.*:.*avx2.*" cpuinfo) != null;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in stdenv.mkDerivation rec {

  pname = "Grid-lehner";
  version = "965cee4372db2a61f0f13546fd521670a614495c";

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/grid" else "https://github.com/lehner/Grid";
    ref = "feature/gpt";
    rev = version;
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
  ++ lib.optional stdenv.cc.isClang openmp
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
    gpu_cflags = "-Xcompiler -fPIC -ccbin mpic++ -arch=${nvcc-arch} -w -cudart shared";
    cpu_ldflags = "";
    gpu_ldflags = "-Xcompiler -fopenmp -ccbin mpic++ -cudart shared -lcublas";
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
    echo 'stdenv=${stdenv.cc}'
    echo
  '' + extra;

  configureFlags = let
    cpu_avx2_flags = [
      "--enable-comms=mpi-auto"
      "--enable-unified=yes"
      # "--enable-shm=shmopen"
      # "--enable-shm-fast-path=shmopen"
      "--enable-accelerator=none"
      "--enable-simd=AVX2"
      "--enable-alloc-align=4k"
      "--disable-fermion-reps"
      "--disable-gparity"
    ];
    cpu_gen16_flags = [
      "--enable-comms=mpi-auto"
      "--enable-unified=yes"
      "--enable-accelerator=none"
      "--enable-simd=GEN"
      "--enable-gen-simd-width=16"
      "--enable-alloc-align=4k"
      "--disable-fermion-reps"
      "--disable-gparity"
    ];
    cpu_flags = if cpuinfo-has-avx2 then cpu_avx2_flags else cpu_gen16_flags;
    gpu_flags = [
      "--enable-comms=mpi-auto"
      "--enable-unified=no"
      "--enable-shm=nvlink"
      "--enable-accelerator=cuda"
      "--enable-simd=GPU"
      "--enable-alloc-align=4k"
      "--disable-fermion-reps"
      "--disable-gparity"
    ];
    flags = if cudaSupport then gpu_flags else cpu_flags;
  in flags;

  preBuild = ''
    cd Grid
    pwd
  '';

  postInstall= ''
    cd ..
    make -j "$NIX_BUILD_CORES" -C benchmarks Benchmark_dwf_fp32
    install -D -m755 grid-config "$out"/bin/grid-config
    install -D -m755 benchmarks/Benchmark_dwf_fp32 "$out"/bin/Benchmark_dwf_fp32
  '';

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
