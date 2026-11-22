#!/bin/bash

# Shared CI debug environment/context dump.
#
# Usage: scripts/ci-debug-dump-env.sh [label]
#
# Called first in every failure()/cancelled() debug step, before
# scripts/ci-debug-collect-logs.sh / scripts/ci-debug-mac-nix-gpt.sh, so all CI
# jobs print the same operating-system / toolchain / Python / MPI context from
# one place (the inline blocks these steps used to have were duplicated across
# jobs).  The `label` argument names the job in the header.
#
# Requirements:
#   * must never fail: every probe is best effort and the script exits 0;
#   * must never leak secrets: only an allowlist of the environment is printed,
#     the full `env` is intentionally NOT dumped (GitHub tokens live there);
#   * works on both Linux (ubuntu runners) and macOS (macos runners); the
#     Darwin-only probes are guarded so running it elsewhere cannot abort.
#
# The calling workflow step keeps its own `set +e` and any job-specific
# `source .../setenv.sh` it needs (so PATH/PYTHONPATH/MPI point at the real
# toolchain), then calls this script and the log collector.

set +e
set +u

label="${1:-debug}"

section() {
    echo
    echo "===== $* ====="
}

# Run a command for its output, never aborting on failure.
probe() {
    "$@" 2>&1 || true
}

section "debug context ($label)"
echo "date (UTC): $(date -u '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || true)"
echo "date (local): $(date 2>/dev/null || true)"
echo "OS: $(uname -s 2>/dev/null || true)"
echo "uname: $(uname -a 2>/dev/null || true)"
echo "shell: ${SHELL:-<unset>}${BASH_VERSION:+ (bash $BASH_VERSION)}"
echo "pwd: $(pwd 2>/dev/null || true)"
echo "runner temp: ${RUNNER_TEMP:-<unset>}"

section "disk"
probe df -h

section "cpu / memory"
os="$(uname -s 2>/dev/null || true)"
case "$os" in
    Darwin)
        echo "--- sysctl -n hw.ncpu: $(sysctl -n hw.ncpu 2>/dev/null || true)"
        echo "--- sysctl -n hw.memsize"
        probe sysctl -n hw.memsize
        echo "--- vm_stat (first 20 lines)"
        vm_stat 2>&1 | head -n 20 || true
        ;;
    *)
        echo "--- nproc: $(nproc 2>/dev/null || true)"
        echo "--- free -h"
        probe free -h
        ;;
esac

if [ "$os" = "Darwin" ] ; then
    section "macOS"
    probe sw_vers
    probe sysctl kern.corefile
    probe ls -la /cores
fi

section "toolchain versions"
for cmd in cc c++ clang++ cmake python3 mpiexec mpirun nix git ; do
    if command -v "$cmd" >/dev/null 2>&1 ; then
        echo "--- $cmd: $(command -v "$cmd" 2>/dev/null || true)"
        case "$cmd" in
            cc|c++|clang++) probe "$cmd" --version ;;
            cmake) probe cmake --version ;;
            python3) probe python3 -VV ;;
            mpiexec|mpirun) probe "$cmd" --version ;;
            nix) probe nix --version ;;
            git) probe git --version ;;
        esac
    else
        echo "--- $cmd: not found"
    fi
done

section "python / MPI"
echo "python3: $(command -v python3 2>/dev/null || true)"
python3 -VV 2>&1 || true
python3 -c "import sys ; print(sys.version) ; print(sys.path)" 2>&1 || true
python3 -c "import multiprocessing as mp ; print('mp.get_start_method() =', mp.get_start_method()) ; print('mp.get_all_start_methods() =', mp.get_all_start_methods())" 2>&1 || true
mpiexec --version 2>&1 || true

section "ulimit"
echo "--- ulimit -c: $(ulimit -c 2>/dev/null || true)"
probe ulimit -a

section "environment (allowlist; the full env is intentionally not dumped)"
# Only names that are useful for debugging the build/test environment and that
# cannot carry GitHub secrets.  Kept as a prefix/suffix allowlist: never print
# all of `env`.
env 2>/dev/null | sort | grep -E '^(QLAT|NIX|MPI|OMPI|OPENMPI|PYTHON|PYTHONPATH|LD_|DYLD_|LIBRARY_PATH|PATH|CXX|CC|CFLAGS|CXXFLAGS|CPPFLAGS|LDFLAGS|HOME|RUNNER_)' || echo "(no matching environment variables)"

section "end of debug context ($label)"
exit 0
