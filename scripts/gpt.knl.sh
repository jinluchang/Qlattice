#!/bin/bash

. scripts/conf.sh

name=gpt

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix"/$name || true

rsync -av --delete $distfiles/$name/ "$prefix"/$name/

cd "$prefix/$name"

git checkout 9ad98ae25493bc9987de427970a3c3eb19b7b742

cd lib/cgpt

echo "BASIS_SIZE(4)" > lib/basis_size.h
echo "BASIS_SIZE(10)" >> lib/basis_size.h
echo "BASIS_SIZE(30)" >> lib/basis_size.h
echo "BASIS_SIZE(50)" >> lib/basis_size.h

echo "SPIN(4)" > lib/spin_color.h
echo "COLOR(3)" >> lib/spin_color.h
echo "SPIN_COLOR(4,3)" >> lib/spin_color.h

./clean

./make "$prefix"/Grid/build "$num_proc"

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} |& tee $prefix/log.$name.txt
