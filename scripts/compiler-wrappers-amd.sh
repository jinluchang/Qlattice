#!/bin/bash

. conf.sh

name=compiler-wrappers-amd

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix/bin"

cat - >"$prefix/bin/compiler-options.py" << EOF
#!/usr/bin/env python3
import sys
import subprocess as p
argv = sys.argv.copy()
arg_to_remove = "-xcore-avx2"
while arg_to_remove in argv:
    argv.remove(arg_to_remove)
p.run(argv)
EOF
chmod +x "$prefix/bin/compiler-options.py"

cat - >"$prefix/bin/CC.sh" << EOF
#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ -f "$prefix/bin/clang" ] ; then
    $run clang "\$@"
elif [ -f "$prefix/bin/gcc" ] ; then
    $run gcc "\$@"
elif which icc >/dev/null 2>&1 ; then
    $run icc "\$@"
else
    $run gcc "\$@"
fi
EOF
chmod +x "$prefix/bin/CC.sh"

cat - >"$prefix/bin/CXX.sh" << EOF
#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ -f "$prefix/bin/clang++" ] ; then
    $run clang++ "\$@"
elif [ -f "$prefix/bin/g++" ] ; then
    $run g++ "\$@"
elif which icpc >/dev/null 2>&1 ; then
    $run icpc "\$@"
else
    $run g++ "\$@"
fi
EOF
chmod +x "$prefix/bin/CXX.sh"

cat - >"$prefix/bin/MPICC.sh" << EOF
#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if which mpiicc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang" ] ; then
        $run mpiicc -cc=clang "\$@"
    elif [ -f "$prefix/bin/gcc" ] ; then
        $run mpiicc -cc=gcc "\$@"
    else
        $run mpiicc "\$@"
    fi
else
    mpicc "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICC.sh"

cat - >"$prefix/bin/MPICXX.sh" << EOF
#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if which mpiicpc >/dev/null 2>&1 ; then
    if [ -f "$prefix/bin/clang++" ] ; then
        $run mpiicpc -cxx=clang++ "\$@"
    elif [ -f "$prefix/bin/g++" ] ; then
        $run mpiicpc -cxx=g++ "\$@"
    else
        $run mpiicpc "\$@"
    fi
elif which mpicxx >/dev/null 2>&1 ; then
    $run mpicxx "\$@"
else
    $run mpic++ "\$@"
fi
EOF
chmod +x "$prefix/bin/MPICXX.sh"

echo "!!!! $name build !!!!"

} |& tee $prefix/log.$name.txt
