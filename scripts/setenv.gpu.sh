#!/bin/bash

name=setenv-gpu

source qcore/set-prefix.sh

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh

#
cat >"$prefix/setenv.sh" <<EOF
if [ -z "\$num_proc" ] ; then
    export num_proc=2
fi
if [ -z "\$NVCC_ARCH" ] ; then
    export NVCC_ARCH="sm_86"
fi
export NVCC_OPTIONS="-w -std=c++14 -arch=\$NVCC_ARCH --expt-extended-lambda --expt-relaxed-constexpr -fopenmp -fno-strict-aliasing" # -D__DEBUG_VECUTILS__
export QLAT_CXX="NVCC.py -ccbin CXX.sh \$NVCC_OPTIONS"
export QLAT_MPICXX="NVCC.py -ccbin MPICXX.sh \$NVCC_OPTIONS"
export QLAT_CXXFLAGS="--NVCC-compile" # -fPIC 
export QLAT_LDFLAGS="--NVCC-link" # --shared
EOF

    #
    "$wd"/qcore/bin/mk-setenv-dir.sh --keep
    echo "!!!! $name build !!!!"
} } 2>&1 | tee $prefix/log.$name.txt
