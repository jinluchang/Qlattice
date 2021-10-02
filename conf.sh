# $prefix controls the installation directory

prefix_default="$HOME/qlat-build/default"

if [ -z "$prefix" ] ; then
    prefix="$prefix_default"
fi

prefix="$(readlink -m "$prefix")"

if [ -f "$prefix/setenv.sh" ] ; then
    . "$prefix/setenv.sh"
    if which mpic++ >/dev/null 2>&1 ; then
        type mpic++
    elif which mpicxx >/dev/null 2>&1 ; then
        type mpicxx
    else
        echo "NO mpic++ or mpicxx available. Quit."
        exit
    fi
else
    echo "'$prefix/setenv.sh' does not exist."
fi

set -e

wd="$(pwd)"
distfiles="$wd/distfiles"
temp_dir="$prefix/temp"
src_dir="$temp_dir/src"
build_dir="$temp_dir/build"

rm -rf $temp_dir || true
