#!/bin/bash

. conf.sh

name=compiler-wrappers

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix/bin"

cat - >"$prefix/bin/CC.sh" << EOF
#!/bin/bash
if [ -f "$prefix/bin/clang" ] ; then
    clang "\$@"
elif [ -f "$prefix/bin/gcc" ] ; then
    gcc "\$@"
elif which icc >/dev/null 2>&1 ; then
    icc "\$@"
else
    gcc "\$@"
fi
EOF
chmod +x "$prefix/bin/CC.sh"

cat - >"$prefix/bin/CXX.sh" << EOF
#!/bin/bash
if [ -f "$prefix/bin/clang++" ] ; then
    clang++ "\$@"
elif [ -f "$prefix/bin/g++" ] ; then
    g++ "\$@"
elif which icpc >/dev/null 2>&1 ; then
    icpc "\$@"
else
    g++ "\$@"
fi
EOF
chmod +x "$prefix/bin/CXX.sh"

cat - >"$prefix/bin/MPICC.sh" << EOF
#!/bin/bash
if which mpiicc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang" ] ; then
        mpiicc -cc=clang "\$@"
    elif [ -f "$prefix/bin/gcc" ] ; then
        mpiicc -cc=gcc "\$@"
    else
        mpiicc "\$@"
    fi
else
    mpicc "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICC.sh"

cat - >"$prefix/bin/MPICXX.sh" << EOF
#!/bin/bash
if which mpiicpc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang++" ] ; then
        mpiicpc -cxx=clang++ "\$@"
    elif [ -f "$prefix/bin/g++" ] ; then
        mpiicpc -cxx=g++ "\$@"
    else
        mpiicpc "\$@"
    fi
elif which mpicxx >/dev/null 2>&1 ; then
    mpicxx "\$@"
else
    mpic++ "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICXX.sh"

echo "!!!! $name build !!!!"

} |& tee $prefix/log.$name.txt
