if [ "$(uname)" == "Darwin" ]; then
    echo "Setting for Mac OS X"
    export q_num_mp_processes=0
    if which brew >/dev/null 2>&1 ; then
        echo "Setting for brew in Mac OS X with prefix: $(brew --prefix)"
        export PATH="$(brew --prefix)/opt/openssl@3/bin":"$PATH"
        export PATH="$(brew --prefix)/opt/findutils/libexec/gnubin":"$PATH"
        export PATH="$(brew --prefix)/opt/coreutils/libexec/gnubin":"$PATH"
        export PATH="$(brew --prefix)/opt/llvm/bin":"$PATH"
        export LD_RUN_PATH="$(brew --prefix)/opt/llvm/lib/c++":"$LD_RUN_PATH"
        export LIBRARY_PATH="$(brew --prefix)/Cellar/fftw/3.3.10_1/lib":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/openssl@3/lib":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib/c++":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/zlib/lib":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/lib":"$LIBRARY_PATH"
        export CPATH="$(brew --prefix)/opt/openssl@3/include":"$CPATH"
        export CPATH="$(brew --prefix)/Cellar/fftw/3.3.10_1/include":"$CPATH"
        export CPATH="$(brew --prefix)/opt/llvm/include":"$CPATH"
        export CPATH="$(brew --prefix)/opt/zlib/include":"$CPATH"
        export CPATH="$(brew --prefix)/opt/include":"$CPATH"
        export PKG_CONFIG_PATH="$(brew --prefix)/opt/zlib/lib/pkgconfig":"$PKG_CONFIG_PATH"
    fi
    if [ -z ${USE_COMPILER+x} ] ; then
        export USE_COMPILER=clang
    fi
    export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
elif [ ! "$(uname)" == "Linux" ]; then
    echo "Setting for $(uname) as if it is a Linux"
fi

if [ -n "$NVCC_ARCH" ] ; then
    echo "Setting use GPU."
    export q_num_mp_processes=0
fi

export LC_ALL="C"

if [ -z "${num_proc}" ] ; then
    export num_proc=2
fi

if [ -z "${num_test}" ] ; then
    export num_test=1
fi

export NINJA_NUM_JOBS="$num_proc"

if [ -z "${USE_COMPILER+x}" ] ; then
    export USE_COMPILER=gcc
fi

if [ -z "${CC+x}" ] ; then
    export CC=CC.sh
fi

if [ -z "${CXX+x}" ] ; then
    export CXX=CXX.sh
fi

if [ -z "${MPICC+x}" ] ; then
    export MPICC=MPICC.sh
fi

if [ -z "${MPICXX+x}" ] ; then
    export MPICXX=MPICXX.sh
fi

if [ -z "${OMP_NUM_THREADS+x}" ] ; then
    export OMP_NUM_THREADS=2
fi

if [ -z "${CUBACORES+x}" ] ; then
    export CUBACORES=0
fi

if [ -z "${PYTHONHOME+x}" ] ; then
    if [ "/usr/bin/python3" = "$(which python3)" ] ; then
        export PYTHONHOME=/usr
    fi
fi

export PATH="$setenv_prefix/bin":"$PATH"
