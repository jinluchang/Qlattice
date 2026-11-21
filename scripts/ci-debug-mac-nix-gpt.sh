#!/usr/bin/env bash

# macOS-Nix SIGTRAP diagnosis for the build_nix_mac_clang CI job.
#
# Usage: scripts/ci-debug-mac-nix-gpt.sh <out-dir> [test-name]
#
# Must run from the repository root on a macOS runner with nix installed.  It
# collects the information a failed run cannot currently show:
#
#   1. macOS / toolchain / debugger environment,
#   2. the full per-test logs preserved by `nix-build --keep-failed`
#      (defeating the Makefile's `tail -n 100` truncation),
#   3. a rerun of the first failing GPT example with Grid `--debug-signals`, so
#      Grid_fatal_signal_handler turns the SIGTRAP into a backtrace,
#   4. a best-effort `lldb --batch -o run -o 'thread backtrace all'` rerun of a
#      minimal qlat->GPT copy-plan reproducer (the known crash site).
#
# The reruns are guarded: a failure here must not fail the job, and anything
# that cannot be determined is printed as unavailable.

set -u

out="${1:?usage: scripts/ci-debug-mac-nix-gpt.sh <out-dir> [test-name]}"
mkdir -p "$out"

env_attr="pkgs-cpsless-clang-ucxless.qlat-env"

cat_file() {
    local f="$1"
    if [ -r "$f" ] ; then
        cat "$f" 2>&1
    elif sudo -n test -r "$f" >/dev/null 2>&1 ; then
        sudo -n cat "$f" 2>&1
    else
        echo "<not readable: $f>"
    fi
}

copy_file() {
    local f="$1"
    local dest="$2"
    mkdir -p "$(dirname "$dest")" 2>/dev/null || true
    if cp -a "$f" "$dest" 2>/dev/null ; then
        :
    elif sudo -n cp -a "$f" "$dest" 2>/dev/null ; then
        :
    fi
}

test_name="${2:-}"
if [ -z "$test_name" ] ; then
    # The first entry of examples-py-gpt/Makefile's tests_gpt list is the first
    # test the CI run executes (and the first one that SIGTRAPs).
    test_name="$(sed -n '/^tests_gpt =/,$p' examples-py-gpt/Makefile | grep -o '[A-Za-z0-9_-]*\.log' | head -n 1 | sed 's/\.log$//')"
fi
echo "first failing GPT test to rerun: '${test_name}'"

echo
echo "==================== macOS / toolchain ===================="
uname -a || true
sw_vers || true
ulimit -a || true
ulimit -c unlimited 2>/dev/null || true
nix --version || true
nix-store --version || true
if command -v lldb >/dev/null 2>&1 ; then
    lldb --version || true
else
    echo "lldb: not found"
fi
sysctl kern.corefile 2>/dev/null || true
ls -la /cores 2>/dev/null || true
echo "system python3: $(command -v python3 || true)"
python3 -VV 2>&1 || true
mpiexec --version 2>&1 || true

echo
echo "==================== preserved nix build dirs ===================="
echo "(nix-build --keep-failed keeps failed build directories; the full logs"
echo " of the failed example tests are inside them)"
shopt -s nullglob
kept_dirs=(/private/tmp/nix-build-* /tmp/nix-build-* "${TMPDIR:-/tmp}"/nix-build-*)
if [ "${#kept_dirs[@]}" -gt 0 ] ; then
    kept_dirs=($(printf '%s\n' "${kept_dirs[@]}" | sort -u))
fi
if [ "${#kept_dirs[@]}" -eq 0 ] ; then
    echo "no nix-build-* directory found; check that the nix-build steps use --keep-failed"
fi
for d in ${kept_dirs[@]+"${kept_dirs[@]}"} ; do
    echo
    echo "--- kept build dir: $d ---"
    ls -la "$d" 2>&1 || true
    tag="$(basename "$d")"
    # Directory listing of the failed tests' .p dirs and their log files.
    while IFS= read -r pdir ; do
        echo "--- test dir: $pdir ---"
        ls -la "$pdir" 2>&1 || true
    done < <(find "$d" -maxdepth 8 -type d -name '*.py.p' 2>/dev/null | sort)
    echo "--- *.log* files under the kept build dir ---"
    find "$d" -maxdepth 8 -type f -name '*.log*' 2>/dev/null | sort | while IFS= read -r f ; do
        ls -la "$f" 2>&1 || true
    done
    # Full untruncated per-test logs (the Makefile only prints `tail -n 100`).
    i=0
    while IFS= read -r f ; do
        echo "===== FULL LOG: $f ====="
        cat_file "$f"
        echo "===== END LOG: $f ====="
        i=$((i + 1))
        copy_file "$f" "$out/kept-nix-logs/$tag/$(printf '%02d' "$i")-$(basename "$f")"
    done < <(find "$d" -maxdepth 8 -type f \( -name 'log.full.txt' -o -name 'log.full' -o -name 'log.check.txt.new' \) 2>/dev/null | sort)
    # The failing command lines and the full nix build log of the derivation.
    drv_path="/nix/store/$(basename "$d" | sed -e 's/^nix-build-//' -e 's/-[0-9][0-9]*$//')"
    if [ -f "$drv_path" ] ; then
        nix_log="$out/nix-build-log-$(basename "$drv_path").txt"
        if nix-store -l "$drv_path" > "$nix_log" 2>/dev/null || nix log "$drv_path" > "$nix_log" 2>/dev/null ; then
            echo "--- failing command lines from '$drv_path' ---"
            grep -a -E 'time \(|mpiexec|lldb|Trace/BPT|SIGTRAP' "$nix_log" | head -n 40 || true
            echo "(full nix build log: $nix_log)"
        else
            echo "could not read the nix build log of '$drv_path'"
        fi
    fi
done

echo
echo "==================== nix GPT test environment ===================="
# Keep the nix result symlink and the rerun work dir OUTSIDE <out-dir>: the
# <out-dir> is uploaded as an artifact and must not follow the env symlink
# into the nix store or contain test caches.
base_dir="$(dirname "$out")"
env_link="$base_dir/qlat-debug-result-gpt-env"
rm -rf "$env_link"
if ! nix-build nixpkgs/q-pkgs.nix -A "$env_attr" -o "$env_link" -j 1 --cores 3 ; then
    echo "could not build '$env_attr' (the earlier job step probably failed); skipping reruns"
    exit 0
fi
env_dir="$env_link"
export PATH="$env_dir/bin:$PATH"
echo "env_dir = $env_dir"

# Run a command with a hard time limit; the nix environment does not always
# provide GNU `timeout` and macOS has no timeout(1) by default.
run_with_timeout() {
    local seconds="$1"
    shift
    if [ -x "$env_dir/bin/timeout" ] ; then
        "$env_dir/bin/timeout" -s KILL "$seconds" "$@"
    elif command -v timeout >/dev/null 2>&1 ; then
        timeout -s KILL "$seconds" "$@"
    elif command -v gtimeout >/dev/null 2>&1 ; then
        gtimeout -s KILL "$seconds" "$@"
    elif command -v perl >/dev/null 2>&1 ; then
        perl -e 'alarm shift; exec @ARGV' "$seconds" "$@"
    else
        "$@"
    fi
}
python3 -VV 2>&1 || true
python3 -c 'import multiprocessing as mp ; print("python mp.get_start_method() =", mp.get_start_method())' 2>&1 || true
mpiexec --version 2>&1 | head -n 5 || true
echo
echo "--- binary / deployment-target info (best effort) ---"
echo "python executable: $(python3 -c 'import sys ; print(sys.executable)' 2>/dev/null || true)"
python_exe="$(python3 -c 'import sys ; print(sys.executable)' 2>/dev/null || true)"
if [ -n "$python_exe" ] && [ -f "$python_exe" ] ; then
    file "$python_exe" 2>&1 || true
    otool -l "$python_exe" 2>/dev/null | grep -A4 'LC_BUILD_VERSION\|LC_VERSION_MIN_MACOSX' | head -n 20 || true
    otool -L "$python_exe" 2>/dev/null | head -n 20 || true
fi
while IFS= read -r f ; do
    echo "--- $f"
    file "$f" 2>&1 || true
    otool -L "$f" 2>/dev/null | head -n 25 || true
done < <(find "$env_dir" -maxdepth 10 \( -name 'cgpt*.so' -o -name 'libGrid*.a' -o -name 'libgpt*.a' \) 2>/dev/null | head -n 5)

work="$base_dir/qlat-debug-gpt-rerun"
rm -rf "$work"
mkdir -p "$work"
if [ ! -f "examples-py-gpt/$test_name.py" ] ; then
    echo "test 'examples-py-gpt/$test_name.py' not found; skipping reruns"
    exit 0
fi
cp "examples-py-gpt/$test_name.py" "$work/" 2>/dev/null || true
cp "examples-py-gpt/$test_name.log.json" "$work/" 2>/dev/null || true

# Grid installs its SIGTRAP/SIGSEGV handler only when --debug-signals is in
# argv before gpt initializes (gpt calls cgpt.init(sys.argv) at import time),
# so the flag has to be added before the test script runs.  qlat accepts extra
# argv entries, so running the original script through runpy is safe.
cat > "$work/mpi4py_debug_wrapper.py" <<'PYEOF'
import runpy
import sys

sys.argv += ["--debug-signals"]
test_script = sys.argv[1]
sys.argv[1:2] = []
runpy.run_path(test_script, run_name="__main__")
PYEOF

echo
echo "==================== rerun '$test_name' with Grid --debug-signals ===================="
run_log="$out/gpt-rerun-debug-signals.log"
(
    cd "$work" || exit 0
    ulimit -c unlimited 2>/dev/null || true
    export PYTHONFAULTHANDLER=1
    run_with_timeout 900 mpiexec -n 2 --oversubscribe --bind-to none bash "$env_dir/bin/bind-gpu-qlat.sh" \
        python3 -m mpi4py ./mpi4py_debug_wrapper.py "./$test_name.py" \
        --test --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1
) > "$run_log" 2>&1
run_status=$?
echo "rerun exit status: $run_status"
echo "--- interesting lines from the rerun (full log: $run_log) ---"
grep -a -i -E 'GRID_ASSERT|Grid_fatal_signal_handler|blame|SIGTRAP|Trace/BPT|signal|backtrace|ERR:|qassert|assert' "$run_log" | head -n 80 || true
echo "--- last 300 lines of the rerun ---"
tail -n 300 "$run_log" || true

echo
echo "==================== lldb reproducer for the copy-plan SIGTRAP ===================="
real_python="$(python3 -c 'import sys ; print(sys.executable)' 2>/dev/null || true)"
python_path="$(python3 -c 'import os, sys ; print(os.pathsep.join(p for p in sys.path if p.endswith("site-packages")))' 2>/dev/null || true)"
if command -v lldb >/dev/null 2>&1 && [ -n "$real_python" ] && [ -f "$real_python" ] ; then
    cat > "$work/repro_gpt_copy_plan.py" <<'PYEOF'
# Minimal reproducer of the macOS-Nix SIGTRAP: the first qlat->GPT copy plan
# (qlat/qlat_gpt.py -> cgpt.copy_create_plan).  --debug-signals makes Grid
# print a backtrace if a fatal signal is raised before lldb sees it.
import sys

sys.argv += ["--grid", "8.8.8.8", "--debug-signals"]
import qlat_gpt as qg
import qlat as q

qg.begin_with_gpt()
gf = q.GaugeField(q.Geometry(q.Coordinate([8, 8, 8, 8])))
gf.set_unit()
plan = qg.gpt_from_qlat(gf)
print("copy plan created:", type(plan), flush=True)
qg.end_with_gpt()
PYEOF
    lldb_log="$out/lldb-repro.log"
    (
        cd "$work" || exit 0
        ulimit -c unlimited 2>/dev/null || true
        export PYTHONFAULTHANDLER=1
        export PYTHONPATH="${python_path}${PYTHONPATH:+:$PYTHONPATH}"
        run_with_timeout 600 lldb --batch \
            -o "settings set target.env-vars PYTHONFAULTHANDLER=1" \
            -o run \
            -o "thread backtrace all" \
            -o quit \
            -- "$real_python" "$work/repro_gpt_copy_plan.py"
    ) > "$lldb_log" 2>&1
    lldb_status=$?
    echo "lldb exit status: $lldb_status"
    echo "--- lldb output (full log: $lldb_log) ---"
    tail -n 400 "$lldb_log" || true
else
    echo "lldb not usable (lldb missing or nix python not found); skipped"
fi

echo
echo "==================== debug run finished ===================="
echo "collected logs are under: $out"
