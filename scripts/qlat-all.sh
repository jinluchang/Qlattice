#!/usr/bin/env bash

source qcore/set-prefix.sh

time {
for subdir in qlat-utils qlat qlat-grid qlat-cps ; do
    ./scripts/"$subdir".sh
    if [ -f "$prefix/$subdir/build-successfully.txt" ] ; then
        echo "$subdir build successfully."
    else
        echo "$subdir failed to build."
        exit 1
    fi
done
}
