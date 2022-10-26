echo
echo "prefix=$prefix"
echo
echo "num_proc=$num_proc"

export QLAT_PREFIX="$prefix"

add-to-colon-list () {
    local name="$1"
    local new_value="$2"
    local value="${!name}"
    if [ -z "$value" ] ; then
        export "$name"="$new_value"
    else
        IFS=':' read -a vs <<< "$value"
        local value=''
        for v in "${vs[@]}" ; do
            if [ "$new_value" != "$v" ] ; then
                value+=:"$v"
            fi
        done
        export "$name"="$new_value""$value"
    fi
}

organize-colon-list() {
    local name="$1"
    local value="${!name}"
    if [ -n "$value" ] ; then
        IFS=':' read -a vs <<< "$value"
        value=''
        for v in "${vs[@]}" ; do
            value="$v":"$value"
        done
        value="${value%:}"
        IFS=':' read -a vs <<< "$value"
        export "$name"=""
        for v in "${vs[@]}" ; do
            add-to-colon-list "$name" "$v"
        done
    fi
}

if [ "$(uname)" == "Darwin" ]; then
    echo "Setting for Mac OS X"
    export q_num_mp_processes=0
    add-to-colon-list PATH "/usr/local/opt/openssl@3/bin"
    add-to-colon-list PATH "/usr/local/opt/llvm/bin"
    add-to-colon-list LD_RUN_PATH "/usr/local/opt/llvm/lib/c++"
    add-to-colon-list LIBRARY_PATH "/usr/local/opt/openssl@3/lib"
    add-to-colon-list LIBRARY_PATH "/usr/local/opt/llvm/lib/c++"
    add-to-colon-list LIBRARY_PATH "/usr/local/opt/llvm/lib"
    add-to-colon-list CPATH "/usr/local/opt/openssl@3/include"
    add-to-colon-list CPATH "/usr/local/opt/llvm/include"
    if [ -z "$USE_COMPILER" ] ; then
        export USE_COMPILER=clang
    fi
    export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
elif [ "$(uname)" == "Linux" ]; then
    echo "Setting for Linux"
else
    echo "Setting for $(uname) as if it is a Linux"
fi

if [ -z "$USE_COMPILER" ] ; then
    export USE_COMPILER=gcc
fi

if [ -z "$CC" ] ; then
    export CC=CC.sh
fi

if [ -z "$CXX" ] ; then
    export CXX=CXX.sh
fi

if [ -z "$CFLAGS" ] ; then
    export CFLAGS=
fi

if [ -z "$CXXFLAGS" ] ; then
    export CXXFLAGS=
fi

if [ -z "$LDFLAGS" ] ; then
    export LDFLAGS=
fi

if [ -z "$LIBS" ] ; then
    export LIBS=
fi

if [ -z "$MPICC" ] ; then
    export MPICC=MPICC.sh
fi

if [ -z "$MPICXX" ] ; then
    export MPICXX=MPICXX.sh
fi

if [ -z "$QLAT_CXX" ] ; then
    export QLAT_CXX=CXX.sh
fi

if [ -z "$QLAT_MPICXX" ] ; then
    export QLAT_MPICXX=MPICXX.sh
fi

if [ -z "$NPY_BLAS_ORDER" ] ; then
    export NPY_BLAS_ORDER=openblas
fi

if [ -z "$NPY_LAPACK_ORDER" ] ; then
    export NPY_LAPACK_ORDER=openblas
fi

if [ -z "$NPY_NUM_BUILD_JOBS" ] ; then
    export NPY_NUM_BUILD_JOBS=$num_proc
fi

add-to-colon-list PATH "$prefix/bin"
for v in "$prefix"/lib/python3*/*-packages ; do
    add-to-colon-list PYTHONPATH "$v"
done
add-to-colon-list PYTHONPATH "$prefix/gpt/lib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib/cgpt/build"
add-to-colon-list LD_RUN_PATH "$prefix/lib"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib"
add-to-colon-list LIBRARY_PATH "$prefix/lib"
add-to-colon-list CPATH "$prefix/include/ncurses"
add-to-colon-list CPATH "$prefix/include"
add-to-colon-list PKG_CONFIG_PATH "$prefix/lib/pkgconfig"

organize-colon-list PATH
organize-colon-list PYTHONPATH
organize-colon-list LD_RUN_PATH
organize-colon-list LD_LIBRARY_PATH
organize-colon-list LIBRARY_PATH
organize-colon-list CPATH
organize-colon-list PKG_CONFIG_PATH

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
