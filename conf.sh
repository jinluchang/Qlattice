. setenv.sh

set -e

if [ -z $num_proc ] ; then
    num_proc=2
fi

wd=$(pwd)
distfiles=$wd/distfiles
temp_dir=$wd/temp
src_dir=$temp_dir/src
build_dir=$temp_dir/build

rm -rf $temp_dir || true
