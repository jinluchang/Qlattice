#!/usr/bin/env bash

name=ninja-script

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    rm -rf "$prefix"/bin
    mkdir -p "$prefix"/bin

    echo "Trying to find ninja-backend"

    ninja_path="$(which ninja-backend || true)"

    if [ -z "$ninja_path" ] ; then
        echo "Cannot find ninja-backend. Trying to find ninja"
        for i in {1..10} ; do
            ninja_path="$(which ninja || true)"
            if [ -z "$ninja_path" ] ; then
                echo "Cannot find usable ninja."
                exit 1
            elif [ "$ninja_path" = "$prefix/bin/ninja" ] ; then
                echo "Found $ninja_path, must be the link created by this script."
                echo "Will remove $ninja_path and find ninja again."
                rm -rfv "$ninja_path"
                sleep 0.1
            else
                echo "Found ninja at $ninja_path"
                break
            fi
        done
    fi

    echo "ninja_path=$ninja_path"

#
cat - >"$prefix/bin/ninja" << EOF
#!/usr/bin/env bash

# Need to limit the number of JOBS by ninja
# See https://github.com/ninja-build/ninja/issues/1482

if [ -z "\$NINJA_NUM_JOBS" ] ; then
    # Default should be no parallelization
    "$ninja_path" -j1 "\$@"
else
    "$ninja_path" -j"\$NINJA_NUM_JOBS" "\$@"
fi
EOF

    #
    chmod +x "$prefix/bin/ninja"

    cat "$prefix/bin/ninja"

    mk-setenv.py
    echo "!!!! $name build !!!!"
    rm -rf $temp_dir || true
} } 2>&1 | tee $prefix/log.$name.txt
