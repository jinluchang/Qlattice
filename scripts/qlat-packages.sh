#!/usr/bin/env bash

# A failed wheel build/install must fail the step: the trailing `2>&1 | tee`
# below would otherwise mask it and the job would fail later at import time.
set -e -o pipefail

name=qlat-packages

./scripts/qlat-packages-uninstall.sh

source qcore/set-prefix.sh $name

{ time {
    echo "!!!! build $name !!!!"
    source qcore/conf.sh ..

    time-run find ~/.cache/pip/wheels -type f || true
    time-run find ~/.cache/pip/wheels -type f -name 'qlat*.whl' -delete || true

    base_prefix="$(readlink -m "$prefix/..")"

    build="$prefix/packages"
    time-run rm -rfv "$build" || true
    mkdir -p "$build"

    cd "$build"

    # opts="--verbose --no-index --no-cache-dir -f $distfiles/python-packages -f $build"
    # opts="--user --verbose --force-reinstall -f $build"
    # opts="--verbose -f $build --force-reinstall --config-settings compile-args=-j$num_proc"
    # opts="--verbose -f $build --force-reinstall --no-build-isolation --config-settings compile-args=-j$num_proc"
    opts="--verbose -f $build --no-build-isolation --config-settings compile-args=-j$num_proc"

    time-run python3 -m build -ns -o "$build" "$wd"/qlat-utils
    time-run pip3 install $opts qlat-utils
    time-run python3 -m build -ns -o "$build" "$wd"/qlat
    time-run pip3 install $opts qlat
    if [ -n "$(find-library.sh libcps.a)" ] ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-cps
        time-run pip3 install $opts qlat-cps
        rm -f "$base_prefix"/cps-unavailable.txt
    else
        # Do not silently skip: record the missing CPS installation so that
        # scripts/qlat-examples-py-cps.sh can skip the CPS tests explicitly
        # (instead of them failing at `import qlat_cps`).  CI builds CPS in an
        # earlier step and scripts/cps.sh now fails loudly when that build
        # breaks, so this normally only happens on platforms that do not
        # support CPS (e.g. macOS).
        echo "WARNING: libcps.a not found under '$prefix'; qlat-cps will NOT be installed."
        echo "WARNING: CPS tests will be skipped (marker: $base_prefix/cps-unavailable.txt)."
        printf '%s\n' \
            "CPS unavailable: libcps.a not found under '$prefix' at $(date)" \
            > "$base_prefix"/cps-unavailable.txt
    fi
    if grid-config --prefix >/dev/null 2>&1 ; then
        time-run python3 -m build -ns -o "$build" "$wd"/qlat-grid
        time-run pip3 install $opts qlat-grid
    fi

    echo "!!!! $name build !!!!"
    rm -rf "$temp_dir" || true
    touch "$prefix"/build-successfully.txt
} } 2>&1 | tee $prefix/log.$name.txt
