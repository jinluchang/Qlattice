add-to-colon-list() {
    local name="$1"
    local new_value="$2"
    local value="${!name}"
    local v
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

if [ "$(uname)" == "Darwin" ]; then
    echo "Setting for Mac OS X"
    export q_num_mp_processes=0
    if which brew >/dev/null 2>&1 ; then
        echo "Setting for brew in Mac OS X with prefix: $(brew --prefix)"
        if [ -e "$(brew --prefix)/opt/openssl@3/bin" ]; then
            export PATH="$(brew --prefix)/opt/openssl@3/bin""${PATH:+:$PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/findutils/libexec/gnubin" ]; then
            export PATH="$(brew --prefix)/opt/findutils/libexec/gnubin""${PATH:+:$PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/coreutils/libexec/gnubin" ]; then
            export PATH="$(brew --prefix)/opt/coreutils/libexec/gnubin""${PATH:+:$PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/gnu-sed/libexec/gnubin" ]; then
            export PATH="$(brew --prefix)/opt/gnu-sed/libexec/gnubin""${PATH:+:$PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/llvm/bin" ]; then
            export PATH="$(brew --prefix)/opt/llvm/bin""${PATH:+:$PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/llvm/lib/c++" ]; then
            export LD_RUN_PATH="$(brew --prefix)/opt/llvm/lib/c++""${LD_RUN_PATH:+:$LD_RUN_PATH}"
        fi
        if [ -e "$(brew --prefix)/Cellar/fftw/"*"/lib" ]; then
            export LIBRARY_PATH="$(brew --prefix)/Cellar/fftw/"*"/lib""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/openssl@3/lib" ]; then
            export LIBRARY_PATH="$(brew --prefix)/opt/openssl@3/lib""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/llvm/lib/c++" ]; then
            export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib/c++""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/llvm/lib" ]; then
            export LIBRARY_PATH="$(brew --prefix)/opt/llvm/lib""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/zlib/lib" ]; then
            export LIBRARY_PATH="$(brew --prefix)/opt/zlib/lib""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/lib" ]; then
            export LIBRARY_PATH="$(brew --prefix)/opt/lib""${LIBRARY_PATH:+:$LIBRARY_PATH}"
        fi
        if [ -e "$(brew --prefix)/opt/openssl@3/include" ]; then
            export CPATH="$(brew --prefix)/opt/openssl@3/include""${CPATH:+:$CPATH}"
        fi
        if [ -e "$(brew --prefix)/Cellar/fftw/"*"/include" ]; then
            export CPATH="$(brew --prefix)/Cellar/fftw/"*"/include""${CPATH:+:$CPATH}"
        fi
        if [ -e "$(brew --prefix)/opt/llvm/include" ]; then
            export CPATH="$(brew --prefix)/opt/llvm/include""${CPATH:+:$CPATH}"
        fi
        if [ -e "$(brew --prefix)/opt/zlib/include" ]; then
            export CPATH="$(brew --prefix)/opt/zlib/include""${CPATH:+:$CPATH}"
        fi
        if [ -e "$(brew --prefix)/opt/include" ]; then
            export CPATH="$(brew --prefix)/opt/include""${CPATH:+:$CPATH}"
        fi
        if [ -e "$(brew --prefix)/opt/zlib/lib/pkgconfig" ]; then
            export PKG_CONFIG_PATH="$(brew --prefix)/opt/zlib/lib/pkgconfig""${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
        fi
    fi
    if [ -z ${USE_COMPILER+x} ] ; then
        export USE_COMPILER=clang
    fi
    export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
else
    if [ "$(uname)" == "Linux" ]; then
        echo "Setting for Linux"
    else
        echo "Setting for $(uname) as if it is a Linux"
    fi
    if [ -d /usr/lib64 ] ; then
        export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}"/usr/lib64
    fi
    if [ -d /lib64 ] ; then
        export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}"/lib64
    fi
    if [ -n "$NIX_CFLAGS_COMPILE" ] ; then
        value="$NIX_CFLAGS_COMPILE"
        IFS=' ' read -a vs <<< "$value"
        local value=''
        for v in "${vs[@]}" ; do
            if [[ "$v" =~ .*"/include" ]] ; then
                add-to-colon-list CPATH "$v"
            fi
        done
    fi
    if [ -n "$NIX_LDFLAGS" ] ; then
        value="$NIX_LDFLAGS"
        IFS=' ' read -a vs <<< "$value"
        local value=''
        for v in "${vs[@]}" ; do
            if [[ "$v" =~ "-L/".*"/lib" ]] ; then
                add-to-colon-list LIBRARY_PATH "${v#-L}"
            fi
        done
    fi
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
