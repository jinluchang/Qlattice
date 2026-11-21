#!/usr/bin/env bash

# Helpers for treating CPS as optional in the non-pip build flow.
#
# CPS cannot be built on macOS (its a2a code needs GSL headers, which the
# macOS runner image does not provide) and the minimal `build.sh` configuration
# on Linux intentionally does not build CPS at all.  In those cases the
# CPS-dependent steps must be skipped *explicitly* (WARNING + the same
# `cps-unavailable.txt` marker the pip flow uses) instead of aborting the whole
# build before qlat/qlat-utils/qlat-grid are installed.  When CPS *is* built
# (the default-gpt flow on Linux), failures stay fatal.
#
# The calling build.default*.sh scripts run with `set -e -o pipefail`, and call
# these helpers as plain commands: cps_build_optional() exits the whole build on
# a fatal failure, while the optional paths return 0 and let the rest of the
# build proceed.

cps_marker_file() {
    echo "$prefix"/cps-unavailable.txt
}

cps_optional_platform() {
    # macOS does not support CPS yet.  QLAT_CPS_OPTIONAL=1 allows simulating /
    # forcing the skip on other platforms (e.g. CI or a minimal build).
    [ "$(uname)" = "Darwin" ] || [ "${QLAT_CPS_OPTIONAL:-0}" = "1" ]
}

cps_library_available() {
    # CPS built by scripts/cps.sh is installed in $prefix/CPS/lib; a system CPS
    # installation is found via find-library.sh when it is on PATH.
    [ -e "$prefix"/CPS/lib/libcps.a ] && return 0
    command -v find-library.sh >/dev/null 2>&1 && [ -n "$(find-library.sh libcps.a)" ]
}

cps_mark_unavailable() {
    local reason="$1"
    echo "WARNING: $reason"
    echo "WARNING: CPS is not available on this platform; continuing the build without CPS/qlat-cps."
    echo "WARNING: marker written to $(cps_marker_file)"
    printf '%s\n' "CPS unavailable: $reason at $(date)" > "$(cps_marker_file)"
}

cps_build_optional() {
    # Usage: cps_build_optional <label> <command> [args...]
    local step="$1"
    shift
    if "$@" ; then
        return 0
    fi
    if cps_optional_platform ; then
        cps_mark_unavailable "'$step' failed on $(uname)"
        return 0
    fi
    echo "ERROR: '$step' failed; CPS is required on this platform."
    exit 1
}

cps_python_package_optional() {
    # qlat-cps is a thin binding to CPS: build it only when the CPS library was
    # actually built in this prefix (or found as a system library).  The plain
    # `build.sh` flow does not build CPS, so a missing library is an explicit
    # skip, not a build failure.
    if cps_library_available ; then
        if ./scripts/qlat-cps.sh ; then
            rm -f "$(cps_marker_file)"
            return 0
        fi
        if cps_optional_platform ; then
            cps_mark_unavailable "'qlat-cps' failed on $(uname)"
            return 0
        fi
        echo "ERROR: 'qlat-cps' failed; CPS is required on this platform."
        exit 1
    else
        if [ -f "$(cps_marker_file)" ] ; then
            echo "WARNING: CPS library (libcps.a) was not built in '$prefix'; skipping qlat-cps."
        else
            cps_mark_unavailable "CPS library (libcps.a) was not built in '$prefix'"
        fi
    fi
}
