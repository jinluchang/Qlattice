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
, cpuinfo-sys ? "flags : avx2"
, use-gitee ? null
}:

let
  pname = "Grid-lehner";
  cpuinfo-has-avx2 = (builtins.match "flags *: .*avx2.*" cpuinfo-sys) != null;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
  #
  # version = "f0573d04c76dd67a0af2ef1ce18ddf6b227567e2"; # 2025/04/17
  # version = "faed6d3e9d8bfa45326a2720b8916228b941b69b"; # 2025/12/12
  # version = "afd4423bef303fd176bcdb9f3ea89255cba06d53"; # 2026/03/11
  # version = "e13370f81f553b06b62dd7c2aaf20677b5d16ede"; # 2026/04/18
  # version = "f62689d905a5fe3445b93eadcf9e3b33aaf1b298"; # 2026/05/22
  version = "25035275741cedf372a7c6ab2fe2d5a8b8a0aceb"; # 2026/05/31
  #
  src-remote = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/Grid-clehner" else "https://github.com/jinluchang/Grid-clehner";
    ref = "feature/gpt";
    rev = version;
  };
  src-local = builtins.fetchGit ../distfiles/Grid-ljin;
in stdenv.mkDerivation {

  inherit pname version;

  src = src-remote;

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
    libcufft
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
    gpu_cflags = "-Xcompiler -fPIC -ccbin ${mpi.dev}/bin/mpic++ -arch=${nvcc-arch} -w -cudart shared";
    cpu_ldflags = "";
    gpu_ldflags = "-Xcompiler -fopenmp -ccbin ${mpi.dev}/bin/mpic++ -cudart shared -lcublas -L${lib.getLib cudaPackages.cuda_cudart}/lib/stubs";
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
      "--enable-gen-simd-width=32"
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
