#!/usr/bin/env bash

# Collect full, untruncated logs of the example tests that failed in a CI run.
#
# Usage: scripts/ci-debug-collect-logs.sh <build-root> <out-dir>
#
# <build-root> is the build prefix used by the direct (non-Nix) CI jobs, e.g.
# $HOME/qlat-build/default; it contains qlat-examples-*/examples-* trees.
# <out-dir> is where the logs are copied for the artifact upload.
#
# The example Makefiles detect a failure by diffing the newly produced CHECK
# output against the committed reference and only print `tail -n 100` of the
# full log, so a CI failure hides the real traceback.  This script finds every
# test whose CHECK output differs from the reference (or was not produced at
# all, e.g. because the test was killed), prints its complete log.full.txt and
# copies the logs into <out-dir>.

set -u

root="${1:?usage: scripts/ci-debug-collect-logs.sh <build-root> <out-dir>}"
out="${2:?usage: scripts/ci-debug-collect-logs.sh <build-root> <out-dir>}"

mkdir -p "$out"

n_checked=0
n_failed=0

show_test() {
    local dir="$1"
    local full="$2"
    local ref="$3"
    local new="$4"
    local tag="$5"
    local reason="$6"
    echo
    echo "########################################################################"
    echo "## $reason"
    echo "## test dir: $dir"
    echo "## reference: $ref"
    echo "## new:       $new"
    echo "########################################################################"
    echo "## directory listing:"
    ls -la "$dir" 2>&1 || true
    echo "## log.check.txt.new ($(wc -c < "$new" 2>/dev/null || echo 0) bytes):"
    if [ -s "$new" ] ; then
        cat "$new"
    else
        echo "<EMPTY/none: the test did not reach its final 'CHECK: finished successfully.' line>"
    fi
    echo "## FULL log (untruncated): $full"
    if [ -f "$full" ] ; then
        cat "$full"
    else
        echo "<missing: $full>"
    fi
    echo "## end of full log: $tag"
    mkdir -p "$out/$tag"
    local f
    for f in "$full" "$new" "$dir/log.txt" ; do
        if [ -f "$f" ] ; then
            cp -a "$f" "$out/$tag/" || true
        fi
    done
}

echo "===== collect failed/incomplete logs; build root = '$root' ====="

# Python / GPT / CPS tests: <test>.py.p/log.full.txt
while IFS= read -r -d '' new ; do
    dir="${new%/log.check.txt.new}"
    ref="$dir/log.check.txt"
    full="$dir/log.full.txt"
    n_checked=$((n_checked + 1))
    if [ -s "$new" ] && [ -f "$ref" ] && diff -q "$ref" "$new" >/dev/null 2>&1 ; then
        continue
    fi
    n_failed=$((n_failed + 1))
    rel="${dir#"$root"/}"
    tag="${rel//\//__}"
    show_test "$dir" "$full" "$ref" "$new" "$tag" "FAILED Python test (CHECK output differs from the reference)"
done < <(find "$root" -name 'log.check.txt.new' -path '*examples-py*' -print0 2>/dev/null | sort -z)

# Tests that were killed before the Makefile could write log.check.txt.new
# (e.g. the CI step hit its timeout-minutes while the test was still running).
while IFS= read -r -d '' full ; do
    dir="${full%/log.full.txt}"
    new="$dir/log.check.txt.new"
    if [ -e "$new" ] ; then
        continue
    fi
    n_checked=$((n_checked + 1))
    n_failed=$((n_failed + 1))
    rel="${dir#"$root"/}"
    tag="${rel//\//__}"
    show_test "$dir" "$full" "$dir/log.check.txt" "$new" "$tag" "KILLED/INCOMPLETE Python test (no log.check.txt.new was produced)"
done < <(find "$root" -name 'log.full.txt' -path '*.py.p/*' -print0 2>/dev/null | sort -z)

# C++ tests: <test>/build/log.check.txt.new and <test>/log.full
while IFS= read -r -d '' new ; do
    build_dir="${new%/log.check.txt.new}"
    dir="${build_dir%/build}"
    ref="$build_dir/log.check.txt"
    full="$dir/log.full"
    n_checked=$((n_checked + 1))
    if [ -s "$new" ] && [ -f "$ref" ] && diff -q "$ref" "$new" >/dev/null 2>&1 ; then
        continue
    fi
    n_failed=$((n_failed + 1))
    rel="${dir#"$root"/}"
    tag="${rel//\//__}"
    show_test "$dir" "$full" "$ref" "$new" "$tag" "FAILED C++ test (CHECK output differs from the reference)"
done < <(find "$root" -name 'log.check.txt.new' -path '*examples-cpp*' -print0 2>/dev/null | sort -z)

# Keep the per-script logs (they contain the pass/fail summary of a whole
# examples-py group) in the artifact as well.
for f in "$root"/qlat-examples-*/log.*.txt ; do
    if [ -f "$f" ] ; then
        mkdir -p "$out/script-logs"
        cp -a "$f" "$out/script-logs/" || true
    fi
done

echo
echo "===== summary: checked $n_checked test log(s), $n_failed failed/incomplete ====="
