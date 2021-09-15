prefix_default="$(readlink -m "$HOME/qlat-build/default")"

if [ -z "$prefix" ] ; then
    prefix="$prefix_default"
fi

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

rm -rf $temp_dir || true
