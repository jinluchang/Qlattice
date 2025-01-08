#!/usr/bin/env bash

name=setenv-gpu

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

    if [ -z "$NVCC_ARCH" ] ; then
        compute_cap="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | sed 's/\.//')"
        if [ -z "$compute_cap" ] ; then
            NVCC_ARCH="sm_70"
        else
            NVCC_ARCH="sm_$compute_cap"
        fi
    fi
    if [ -z "$CUDA_ROOT" ] ; then
        NVCC_PATH="$(which nvcc)"
        CUDA_ROOT="${NVCC_PATH%/bin/nvcc}"
    fi
    if [ -z "$CUDA_LIBDIR" ] ; then
        CUDA_LIBDIR="/usr/lib/x86_64-linux-gnu/stubs"
    fi
#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=2
fi
if [ -z "\$NVCC_ARCH" ] ; then
    export NVCC_ARCH=$NVCC_ARCH
fi
export NVCC_OPTIONS="-w -std=c++14 -arch=\$NVCC_ARCH --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing" # -D__DEBUG_VECUTILS__
export QLAT_CXX="NVCC.py -ccbin CXX.sh \$NVCC_OPTIONS"
export QLAT_MPICXX="NVCC.py -ccbin MPICXX.sh \$NVCC_OPTIONS"
export QLAT_CXXFLAGS="--NVCC-compile -D__QLAT_BARYON_SHARED_SMALL__" # -fPIC
export QLAT_LDFLAGS="--NVCC-link" # --shared
export CUDA_ROOT="$CUDA_ROOT"
export CUDA_LIBDIR="$CUDA_LIBDIR"
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
