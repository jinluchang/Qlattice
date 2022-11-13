#!/bin/bash

. scripts/conf.sh

name=ninja-script

{

echo "!!!! build $name !!!!"

mkdir -p $prefix/bin

ninja_path="$(which ninja-backend || true)"

if [ -z "$ninja_path" ] ; then
    ninja_path="$(which ninja)"
    while [ "$ninja_path" = "$prefix/bin/ninja" ] ; do
        rm -rfv "$ninja_path"
        sleep 0.1
        ninja_path="$(which ninja || true)"
        if [ -z "$ninja_path" ] ; then
            echo "Cannot find usable ninja."
            exit 1
        fi
    done
fi

echo "ninja_path=$ninja_path"

cat - >"$prefix/bin/ninja" << EOF
#!/bin/bash

# Need to limit the number of JOBS by ninja
# See https://github.com/ninja-build/ninja/issues/1482

if [ -z "\$num_proc" ] ; then
    "$ninja_path" "\$@"
else
    "$ninja_path" -j"\$num_proc" "\$@"
fi
EOF
chmod +x "$prefix/bin/ninja"

cat "$prefix/bin/ninja"

cd $wd
echo "!!!! $name build !!!!"

rm -rf $temp_dir || true

} 2>&1 | tee $prefix/log.$name.txt
