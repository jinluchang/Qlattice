#!/bin/bash

. conf.sh

name=compiler-wrappers

echo "!!!! build $name !!!!"

mkdir -p "$prefix/bin"

cat - >"$prefix/bin/CC" << EOF
#!/bin/bash
if [ -f "$prefix/bin/clang" ] ; then
    clang "\$@"
elif which icc >/dev/null 2>&1 ; then
    icc "\$@"
else
    gcc "\$@"
fi
EOF
chmod +x "$prefix/bin/CC"

cat - >"$prefix/bin/CXX" << EOF
#!/bin/bash
if [ -f "$prefix/bin/clang++" ] ; then
    clang++ "\$@"
elif which icpc >/dev/null 2>&1 ; then
    icpc "\$@"
else
    g++ "\$@"
fi
EOF
chmod +x "$prefix/bin/CXX"

cat - >"$prefix/bin/MPICC" << EOF
#!/bin/bash
if which mpiicc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang" ] ; then
        mpiicc -cc=clang "\$@"
    else
        mpiicc "\$@"
    fi
else
    mpicc "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICC"

cat - >"$prefix/bin/MPICXX" << EOF
#!/bin/bash
if which mpiicpc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang++" ] ; then
        mpiicpc -cxx=clang++ "\$@"
    else
        mpiicpc "\$@"
    fi
elif which mpicxx >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang++" ] ; then
        mpicxx -cxx=clang++ "\$@"
    else
        mpicxx "\$@"
    fi
else
    mpic++ "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICXX"

echo "!!!! $name build !!!!"
