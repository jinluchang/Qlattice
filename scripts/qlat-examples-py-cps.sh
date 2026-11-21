#!/usr/bin/env bash

name=qlat-examples-py-cps

source qcore/set-prefix.sh $name

# If the CPS build/installation was skipped (marker written by
# scripts/qlat-packages.sh), skip the CPS tests explicitly instead of letting
# every test fail at `import qlat_cps`.
base_prefix="$(readlink -m "$prefix/..")"
if [ -f "$base_prefix"/cps-unavailable.txt ] ; then
    echo "WARNING: CPS is not installed for prefix '$base_prefix'; skipping CPS tests."
    cat "$base_prefix"/cps-unavailable.txt
    exit 0
fi

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run rsync -a --delete "$wd"/examples-py-cps "$prefix"/

    export mpi_options="--oversubscribe $mpi_options"

    q_verbose=1 time-run make -C "$prefix"/examples-py-cps update-sources || true
    q_verbose=1 time-run make -C "$prefix"/examples-py-cps run-cps -j "$num_test" || true

    cd "$wd"

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee $prefix/log.$name.txt
