# $prefix controls the installation directory

prefix_default="$HOME/qlat-build/default"

if [ -z "$prefix" ] ; then
    prefix="$prefix_default"
fi

prefix="$(readlink -m "$prefix")"

if [ -f "$prefix/setenv.sh" ] ; then
    . "$prefix/setenv.sh"
else
    echo "'$prefix/setenv.sh' does not exist."
fi

set -e

wd="$(pwd)"
distfiles="$wd/distfiles"
temp_dir="$prefix/temp"
src_dir="$temp_dir/src"
build_dir="$temp_dir/build"

if [ -e $temp_dir ] ; then
    temp_dir_tmp=$(mktemp -d $temp_dir.XXXXX)
    mv $temp_dir $temp_dir_tmp
    rm -rf $temp_dir.* || true
fi
