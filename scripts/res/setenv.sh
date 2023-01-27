echo
echo "prefix=$prefix"
echo
echo "num_proc=$num_proc"

export QLAT_PREFIX="$prefix"

if [ "$(uname)" == "Darwin" ]; then
    echo "Setting for Mac OS X"
    export q_num_mp_processes=0
    export PATH="/usr/local/opt/openssl@3/bin":"$PATH"
    export PATH="/usr/local/opt/llvm/bin":"$PATH"
    export PATH="/usr/local/opt/findutils/libexec/gnubin":"$PATH"
    export LD_RUN_PATH="/usr/local/opt/llvm/lib/c++":"$LD_RUN_PATH"
    export LIBRARY_PATH="/usr/local/opt/openssl@3/lib":"$LIBRARY_PATH"
    export LIBRARY_PATH="/usr/local/opt/llvm/lib/c++":"$LIBRARY_PATH"
    export LIBRARY_PATH="/usr/local/opt/llvm/lib":"$LIBRARY_PATH"
    export CPATH="/usr/local/opt/openssl@3/include":"$CPATH"
    export CPATH="/usr/local/opt/llvm/include":"$CPATH"
    if [ -z ${USE_COMPILER+x} ] ; then
        export USE_COMPILER=clang
    fi
    export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
elif [ "$(uname)" == "Linux" ]; then
    echo "Setting for Linux"
else
    echo "Setting for $(uname) as if it is a Linux"
fi

if [ -z ${USE_COMPILER+x} ] ; then
    export USE_COMPILER=gcc
fi

if [ -z ${CC+x} ] ; then
    export CC=CC.sh
fi

if [ -z ${CXX+x} ] ; then
    export CXX=CXX.sh
fi

if [ -z ${CFLAGS+x} ] ; then
    export CFLAGS=
fi

if [ -z ${CXXFLAGS+x} ] ; then
    export CXXFLAGS=
fi

if [ -z ${LDFLAGS+x} ] ; then
    export LDFLAGS=
fi

if [ -z ${LIBS+x} ] ; then
    export LIBS=
fi

if [ -z ${MPICC+x} ] ; then
    export MPICC=MPICC.sh
fi

if [ -z ${MPICXX+x} ] ; then
    export MPICXX=MPICXX.sh
fi

if [ -z ${QLAT_CXX+x} ] ; then
    export QLAT_CXX=CXX.sh
fi

if [ -z ${QLAT_MPICXX+x} ] ; then
    export QLAT_MPICXX=MPICXX.sh
fi

if [ -z ${NPY_BLAS_ORDER+x} ] ; then
    export NPY_BLAS_ORDER=openblas
fi

if [ -z ${NPY_LAPACK_ORDER+x} ] ; then
    export NPY_LAPACK_ORDER=openblas
fi

if [ -z ${NPY_NUM_BUILD_JOBS+x} ] ; then
    export NPY_NUM_BUILD_JOBS=$num_proc
fi

if [ -z ${NINJA_NUM_JOBS+x} ] ; then
    export NINJA_NUM_JOBS=$num_proc
fi

if [ -z ${q_num_threads+x} ] ; then
    export q_num_threads=2
fi

export PATH="$HOME/.local/bin":"$PATH"
export PATH="$prefix/bin":"$PATH"
if [ -d "$prefix"/lib/python3/qlat-packages ] ; then
    export PYTHONPATH="$prefix"/lib/python3/qlat-packages:"$PYTHONPATH"
fi
export PYTHONPATH="$prefix/gpt/lib":"$PYTHONPATH"
export PYTHONPATH="$prefix/gpt/lib/cgpt/build":"$PYTHONPATH"
export LD_RUN_PATH="$prefix/lib":"$LD_RUN_PATH"
export LD_RUN_PATH="$prefix/lib64":"$LD_RUN_PATH"
export LD_LIBRARY_PATH="$prefix/lib":"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$prefix/lib64":"$LD_LIBRARY_PATH"
export LIBRARY_PATH="$prefix/lib":"$LIBRARY_PATH"
export LIBRARY_PATH="$prefix/lib64":"$LIBRARY_PATH"
export C_INCLUDE_PATH="$prefix/include":"$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH="$prefix/include":"$CPLUS_INCLUDE_PATH"
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig":"$PKG_CONFIG_PATH"
export PKG_CONFIG_PATH="$prefix/lib64/pkgconfig":"$PKG_CONFIG_PATH"

if which organize-colon-list.py >/dev/null 2>&1 ; then
    export PATH="$(organize-colon-list.py "$PATH")"
    export PYTHONPATH="$(organize-colon-list.py "$PYTHONPATH")"
    export LD_RUN_PATH="$(organize-colon-list.py "$LD_RUN_PATH")"
    export LD_LIBRARY_PATH="$(organize-colon-list.py "$LD_LIBRARY_PATH")"
    export LIBRARY_PATH="$(organize-colon-list.py "$LIBRARY_PATH")"
    export CPATH="$(organize-colon-list.py "$CPATH")"
    export C_INCLUDE_PATH="$(organize-colon-list.py "$C_INCLUDE_PATH")"
    export CPLUS_INCLUDE_PATH="$(organize-colon-list.py "$CPLUS_INCLUDE_PATH")"
    export PKG_CONFIG_PATH="$(organize-colon-list.py "$PKG_CONFIG_PATH")"
fi

echo
for v in \
    PATH PYTHONPATH NPY_BLAS_ORDER NPY_NUM_BUILD_JOBS NPY_LAPACK_ORDER \
    LD_PRELOAD LD_RUN_PATH LD_LIBRARY_PATH LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH PKG_CONFIG_PATH \
    CC CXX CFLAGS CXXFLAGS LDFLAGS LIBS MPICC MPICXX \
    QLAT_PREFIX QLAT_CXX QLAT_MPICXX QLAT_CXXFLAGS QLAT_LDFLAGS QLAT_LIBS \
    USE_COMPILER \
    ; do
export | grep --color=never " $v="'"' || true
done
echo

unset v
