# $prefix controls the installation directory

prefix_default="$HOME/qlat-build/default"

if [ -z "$prefix" ] ; then
    prefix="$prefix_default"
fi

if readlink -m "$prefix" >/dev/null 2>&1 ; then
    prefix="$(readlink -m "$prefix")"
else
    echo "Cannot detect real location of '$prefix', use as it is."
fi

if [ -f "$prefix/setenv.sh" ] ; then
    . "$prefix/setenv.sh"
else
    : echo "'$prefix/setenv.sh' does not exist yet."
fi

set -e

wd="$(pwd)"
distfiles="$wd/distfiles"
temp_dir="$prefix/temp"
src_dir="$temp_dir/src"
build_dir="$temp_dir/build"

if [ -e $temp_dir ] ; then
    temp_dir_tmp=$(mktemp -d $temp_dir.tmp.XXXXX)
    mv $temp_dir $temp_dir_tmp
    rm -rf $temp_dir.* || true
fi
