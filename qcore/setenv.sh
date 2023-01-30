if [ "$(uname)" == "Darwin" ]; then
    echo "Setting for Mac OS X"
    export q_num_mp_processes=0
    if which brew >/dev/null 2>&1 ; then
        echo "Setting for brew in Mac OS X with prefix: $(brew --prefix)"
        export PATH="$(brew --prefix)/opt/openssl@3/bin":"$PATH"
        export PATH="$(brew --prefix)/opt/llvm/bin":"$PATH"
        export PATH="$(brew --prefix)/opt/findutils/libexec/gnubin":"$PATH"
        export LD_RUN_PATH="$(brew --prefix)/opt/llvm/lib/c++":"$LD_RUN_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/openssl@3/lib":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib/c++":"$LIBRARY_PATH"
        export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib":"$LIBRARY_PATH"
        export CPATH="$(brew --prefix)/opt/openssl@3/include":"$CPATH"
        export CPATH="$(brew --prefix)/opt/llvm/include":"$CPATH"
    fi
    if [ -z ${USE_COMPILER+x} ] ; then
        export USE_COMPILER=clang
    fi
    export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
elif [ ! "$(uname)" == "Linux" ]; then
    echo "Setting for $(uname) as if it is a Linux"
fi

if [ -z "${num_proc+x}" ] ; then
    export num_proc=2
fi

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

if [ -z "${NINJA_NUM_JOBS+x}" ] ; then
    export NINJA_NUM_JOBS="$num_proc"
fi

if [ -z "${OMP_NUM_THREADS+x}" ] ; then
    export OMP_NUM_THREADS=2
fi

export PATH="$setenv_prefix/bin":"$PATH"
