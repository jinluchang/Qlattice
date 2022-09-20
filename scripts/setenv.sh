echo
echo "prefix=$prefix"
echo
echo "num_proc=$num_proc"

export QLAT_PREFIX="$prefix"

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

add-to-colon-list PATH "$prefix/bin"
for v in "$prefix"/lib/python3*/*-packages ; do
    add-to-colon-list PYTHONPATH "$v"
done
add-to-colon-list PYTHONPATH "$prefix/gpt/lib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib/cgpt/build"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib64"
add-to-colon-list LIBRARY_PATH "$prefix/lib"
add-to-colon-list LIBRARY_PATH "$prefix/lib64"
add-to-colon-list C_INCLUDE_PATH "$prefix/include/ncurses"
add-to-colon-list C_INCLUDE_PATH "$prefix/include"
add-to-colon-list CPLUS_INCLUDE_PATH "$prefix/include/ncurses"
add-to-colon-list CPLUS_INCLUDE_PATH "$prefix/include"
add-to-colon-list PKG_CONFIG_PATH "$prefix/lib/pkgconfig"
add-to-colon-list PKG_CONFIG_PATH "$prefix/lib64/pkgconfig"

organize-colon-list PATH
organize-colon-list PYTHONPATH
organize-colon-list LD_LIBRARY_PATH
organize-colon-list LIBRARY_PATH
organize-colon-list C_INCLUDE_PATH
organize-colon-list CPLUS_INCLUDE_PATH
organize-colon-list PKG_CONFIG_PATH

echo
for v in \
    PATH PYTHONPATH LD_LIBRARY_PATH LIBRARY_PATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH PKG_CONFIG_PATH \
    CC CXX CFLAGS CXXFLAGS LDFLAGS LIBS MPICC MPICXX \
    QLAT_PREFIX QLAT_MPICXX QLAT_CXXFLAGS QLAT_LDFLAGS QLAT_LIBS \
    USE_COMPILER \
    ; do
export | grep --color=never " $v="'"' || true
done
echo
