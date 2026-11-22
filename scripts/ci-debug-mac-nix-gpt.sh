#!/usr/bin/env bash

# macOS-Nix SIGTRAP diagnosis for the build_nix_mac_clang CI job.
#
# Usage: scripts/ci-debug-mac-nix-gpt.sh <out-dir> [test-name]
#
# Set CI_DEBUG_LIGHT=1 for the cheap collection only: steps 1-2 above plus the
# kept-build-dir / nix-store-log / crash-report collection, without the nix env
# build and the two reruns (3-4).  The workflow uses that mode for the
# cancel/failure artifact so a cancelled job does not spend 5-25 min in lldb;
# the heavy repro runs as a separate, gated step after the artifact upload.
#
# The light mode is deliberately build-free: it only reads files under
# $RUNNER_TEMP and queries nix (`nix config show`, `nix-store -l`, which
# fetches an already-recorded build log and never builds) -- so cancelling a
# job can never trigger additional nix compilation.  The heavy nix-build in
# step 3 is reachable only from the non-cancelled "tests actually failed"
# step, where the environment it names was already built by the preceding
# gating build step (the nix-build is a pure store lookup there, no compiler
# is started for a cancelled/failed derivation).
#
# Must run from the repository root on a macOS runner with nix installed.  It
# collects the information a failed run cannot show:
#
#   1. macOS / toolchain / debugger environment,
#   2. the full per-test logs preserved by `nix-build --keep-failed`: on nix >=
#      2.19 the kept build directory lives under the nix `build-dir` setting
#      (e.g. /nix/var/nix/builds/nix-...), not under $TMPDIR, so it is located
#      by parsing `note: keeping build directory "..."` from the
#      `$RUNNER_TEMP/nix-build-logs/*.log` tees of the workflow's nix-build
#      steps (with `nix config show build-dir` and the legacy globs as
#      fallbacks), together with the full `nix-store -l <drv>` build log,
#   3. a rerun of the first failing GPT example with Grid `--debug-signals`, so
#      Grid_fatal_signal_handler turns the SIGTRAP into a backtrace,
#   4. a rerun of a minimal qlat->GPT copy-plan reproducer under lldb.  The
#      target is the *real* Mach-O interpreter resolved from
#      `sys.base_prefix` (the nix `$env/bin/python3` is an exec wrapper, which
#      made the previous attempt stop at `exec` in `_dyld_start` and never run
#      the reproducer), `stop-on-exec` is disabled, SIGTRAP is handled by lldb
#      (the reproducer does not pass Grid's `--debug-signals`, which installs
#      Grid's handler first), the environment is exported in the shell
#      (inherited by lldb) instead of `target.env-vars`, and both
#      `thread backtrace all` and `disassemble --frame` are captured,
#   5. the macOS crash reports (*.ips), which give the faulting thread/stack
#      without a debugger.
#
# The reruns are guarded: a failure here must not fail the job, and anything
# that cannot be determined is printed as unavailable.

set -u

out="${1:?usage: scripts/ci-debug-mac-nix-gpt.sh <out-dir> [test-name]}"
mkdir -p "$out"

env_attr="pkgs-clang-ucxless.qlat-env"

# The workflow tees every nix-build step into this directory.
nix_build_logs_dir="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/nix-build-logs"

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

list_dir() {
    local d="$1"
    if ls -la "$d" 2>&1 ; then
        :
    elif sudo -n ls -la "$d" 2>&1 ; then
        :
    else
        echo "<not readable: $d>"
    fi
}

list_file() {
    local f="$1"
    if ls -la "$f" 2>&1 ; then
        :
    elif sudo -n ls -la "$f" 2>&1 ; then
        :
    else
        echo "<not readable: $f>"
    fi
}

# Run a command with a hard time limit; the nix environment does not always
# provide GNU `timeout` and macOS has no timeout(1) by default.
run_with_timeout() {
    local seconds="$1"
    shift
    if [ -n "${env_dir:-}" ] && [ -x "$env_dir/bin/timeout" ] ; then
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

# Resolve the genuine Mach-O python interpreter, not the nix exec wrapper.
# `sys.base_prefix` is the real store python; `sys.executable` may be the tiny
# wrapper that execs it (which is what made lldb stop at `exec`).
resolve_real_python() {
    if [ -n "${CI_DEBUG_REAL_PYTHON:-}" ] ; then
        echo "$CI_DEBUG_REAL_PYTHON"
        return 0
    fi
    local cand
    cand="$(python3 - <<'PYEOF' 2>/dev/null
import os
import sys

cand = os.path.join(sys.base_prefix, "bin", "python%d.%d" % sys.version_info[:2])
if os.path.exists(cand):
    print(cand)
elif sys.base_prefix != sys.prefix:
    cand2 = os.path.join(sys.base_prefix, "bin", "python3")
    if os.path.exists(cand2):
        print(cand2)
elif os.path.realpath(sys.executable) != sys.executable:
    print(os.path.realpath(sys.executable))
PYEOF
)"
    if [ -z "$cand" ] ; then
        cand="$(python3 -c 'import os, sys ; print(os.path.realpath(sys.executable))' 2>/dev/null || true)"
    fi
    printf '%s\n' "$cand"
}

# Copy the raw nix-build step logs into the artifact, so the next run keeps the
# `note: keeping build directory "..."` lines even if the kept dir is gone.
copy_nix_build_step_logs() {
    local dest="$out/nix-build-step-logs"
    local f
    shopt -s nullglob
    for f in "$nix_build_logs_dir"/*.log ; do
        copy_file "$f" "$dest/$(basename "$f")"
    done
}

# Copy the macOS crash reports of the trapped process (no debugger needed).
collect_crash_reports() {
    echo
    echo "==================== macOS crash reports ===================="
    local dest="$out/crash-reports"
    local src f
    local copied=0
    mkdir -p "$dest" 2>/dev/null || true
    for src in "$HOME/Library/Logs/DiagnosticReports" "/Library/Logs/DiagnosticReports" ; do
        [ -d "$src" ] || continue
        while IFS= read -r f ; do
            if copy_file "$f" "$dest/$(basename "$f")" ; then
                copied=$((copied + 1))
            fi
        done < <(find "$src" -maxdepth 1 -type f -name '*.ips*' 2>/dev/null | sort)
    done
    echo "copied $copied crash report(s) to $dest"
}

# Print the full untruncated logs of one kept build directory and copy them
# into the artifact.
process_kept_build_dir() {
    local d="$1"
    local tag
    tag="$(basename "$d")"
    echo
    echo "--- kept build dir: $d ---"
    list_dir "$d"
    echo "--- *.log* files under the kept build dir ---"
    while IFS= read -r f ; do
        list_file "$f"
    done < <(find "$d" -maxdepth 8 -type f -name '*.log*' 2>/dev/null | sort)
    # Full untruncated per-test logs (the Makefile only prints `tail -n 100`).
    local i=0
    local f
    while IFS= read -r f ; do
        echo "===== FULL LOG: $f ====="
        cat_file "$f"
        echo "===== END LOG: $f ====="
        i=$((i + 1))
        copy_file "$f" "$out/kept-nix-logs/$tag/$(printf '%02d' "$i")-$(basename "$f")"
    done < <(find "$d" -maxdepth 10 -type f \( -name 'log.full.txt' -o -name 'log.full' -o -name 'log.check.txt.new' \) 2>/dev/null | sort)
    # Directory listings of the failed *.py.p dirs (contain log.json + logs).
    local pdir
    while IFS= read -r pdir ; do
        echo "--- test dir: $pdir ---"
        list_dir "$pdir"
    done < <(find "$d" -maxdepth 10 -type d -name '*.py.p' 2>/dev/null | sort)
    # The full nix build log of the failing derivation: use the derivation named
    # in the step log's `error: Cannot build '<drv>'` line, else reconstruct it
    # from the legacy `nix-build-<drv-name>/` directory name.
    local drv=""
    local logf
    shopt -s nullglob
    for logf in "$nix_build_logs_dir"/*.log ; do
        drv="$(sed -n "s/^error: Cannot build '\(.*\.drv\)'.*/\1/p" "$logf" 2>/dev/null | head -n 1)"
        if [ -n "$drv" ] ; then
            break
        fi
    done
    if [ -z "$drv" ] ; then
        local drv_path="/nix/store/$(basename "$d" | sed -e 's/^nix-build-//' -e 's/-[0-9][0-9]*$//')"
        if [ -f "$drv_path" ] ; then
            drv="$drv_path"
        fi
    fi
    if [ -n "$drv" ] ; then
        local nix_log="$out/nix-build-log-$(basename "$drv").txt"
        # `nix-store -l` / `--print-log` only reads the log nix already recorded
        # for this (failed) derivation; it never builds or retries it.  (The
        # newer `nix log` command is not used here on purpose: the debug path
        # must never be able to start a nix build after a cancellation.)
        if nix-store -l "$drv" > "$nix_log" 2>/dev/null ; then
            echo "--- failing command lines from '$drv' ---"
            grep -a -E 'time \(|mpiexec|lldb|Trace/BPT|SIGTRAP' "$nix_log" | head -n 40 || true
            echo "(full nix build log: $nix_log)"
        else
            echo "could not read the nix build log of '$drv'"
        fi
    fi
}

# Locate and report every kept build directory, parsing the `note: keeping
# build directory "..."` lines from the teed nix-build step logs first (that is
# the authoritative source for the current run), then falling back to the nix
# `build-dir` setting and the legacy $TMPDIR/nix-build-* locations.
collect_kept_build_dirs() {
    echo
    echo "==================== preserved nix build dirs ===================="
    echo "(nix-build --keep-failed keeps failed build directories; the full logs"
    echo " of the failed example tests are inside them)"
    shopt -s nullglob
    local kept_dirs=()
    local d f
    local note_re='^note: keeping build directory "\(.*\)"$'
    for f in "$nix_build_logs_dir"/*.log ; do
        [ -f "$f" ] || continue
        while IFS= read -r d ; do
            [ -n "$d" ] || continue
            kept_dirs+=("$d")
            echo "from $(basename "$f"): keeping build directory $d"
        done < <(sed -n "s/${note_re}/\1/p" "$f" 2>/dev/null)
    done
    local build_dir
    build_dir="$(nix config show 2>/dev/null | awk -F' = ' '$1=="build-dir"{print $2}')"
    if [ -n "$build_dir" ] ; then
        for d in "$build_dir"/nix-* ; do
            [ -d "$d" ] || continue
            kept_dirs+=("$d")
            echo "from nix config show build-dir: keeping build directory $d"
        done
    fi
    for d in /private/tmp/nix-build-* /tmp/nix-build-* "${TMPDIR:-/tmp}"/nix-build-* ; do
        [ -d "$d" ] || continue
        kept_dirs+=("$d")
    done
    if [ "${#kept_dirs[@]}" -gt 0 ] ; then
        kept_dirs=($(printf '%s\n' "${kept_dirs[@]}" | sort -u))
    fi
    if [ "${#kept_dirs[@]}" -eq 0 ] ; then
        echo "no kept build directory found; check that the nix-build steps still"
        echo "use --keep-failed and tee into '$nix_build_logs_dir'"
        return 0
    fi
    for d in ${kept_dirs[@]+"${kept_dirs[@]}"} ; do
        process_kept_build_dir "$d"
    done
}

# Rerun the minimal copy-plan reproducer under lldb, with SIGTRAP handled by
# lldb and the real Mach-O python as the target.  One rank is enough (the trap
# is per-rank and rank-local), and the reproducer must not use --debug-signals
# (Grid would install its own SIGTRAP handler and could swallow the trap before
# lldb sees it).
run_lldb_repro() {
    echo
    echo "==================== lldb reproducer for the copy-plan SIGTRAP ===================="
    if ! command -v lldb >/dev/null 2>&1 ; then
        echo "lldb: not found; skipped"
        return 0
    fi
    local real_python
    real_python="$(resolve_real_python)"
    if [ -z "$real_python" ] || [ ! -f "$real_python" ] ; then
        echo "lldb not usable (no python interpreter found); skipped"
        return 0
    fi
    if file "$real_python" 2>/dev/null | grep -q 'Mach-O' ; then
        echo "lldb target (Mach-O): $real_python"
    else
        echo "WARNING: '$real_python' does not look like a Mach-O binary;"
        echo "         lldb may still stop early (expected: nix exec wrapper)"
    fi
    cat > "$work/repro_gpt_copy_plan.py" <<'PYEOF'
# Minimal reproducer of the macOS-Nix SIGTRAP: the first qlat->GPT copy plan
# (qlat/qlat_gpt.py -> cgpt.copy_create_plan).  Do NOT pass --debug-signals
# here: Grid's SIGTRAP handler would be installed first and could swallow the
# trap before lldb sees it.
import sys

sys.argv += ["--grid", "8.8.8.8"]
import qlat_gpt as qg
import qlat as q

qg.begin_with_gpt()
gf = q.GaugeField(q.Geometry(q.Coordinate([8, 8, 8, 8])))
gf.set_unit()
plan = qg.gpt_from_qlat(gf)  # dies here today (qlat_gpt.py:113 -> cgpt copy_create_plan)
print("copy plan created:", type(plan), flush=True)
qg.end_with_gpt()
PYEOF
    local python_path
    python_path="$(python3 -c 'import sys ; print(":".join(p for p in sys.path if p))' 2>/dev/null || true)"
    local lldb_log="$out/lldb-repro.log"
    echo "lldb command:"
    echo "  lldb --batch -o 'settings set target.process.stop-on-exec false' \\"
    echo "       -o 'process handle SIGTRAP -s true -n false -p false' \\"
    echo "       -o run -o 'thread backtrace all' -o 'disassemble --frame' \\"
    echo "       -- $real_python -m mpi4py $work/repro_gpt_copy_plan.py"
    (
        cd "$work" || exit 0
        ulimit -c unlimited 2>/dev/null || true
        # The environment is exported here (inherited by lldb and the inferior);
        # do not use `settings set target.env-vars`, which replaces it.
        export PYTHONFAULTHANDLER=1
        export OMP_NUM_THREADS=2
        export PYTHONPATH="${python_path}${PYTHONPATH:+:$PYTHONPATH}"
        run_with_timeout 600 lldb --batch \
            -o "settings set target.process.stop-on-exec false" \
            -o "process handle SIGTRAP -s true -n false -p false" \
            -o "process handle SIGABRT -s true -n false -p false" \
            -o "run" \
            -o "thread list" \
            -o "thread backtrace all" \
            -o "frame select 0" \
            -o "register read pc sp fp x0 x1" \
            -o "disassemble --frame --count 8" \
            -o "process status" \
            -o "quit" \
            -- "$real_python" -m mpi4py "$work/repro_gpt_copy_plan.py"
    ) > "$lldb_log" 2>&1
    local lldb_status=$?
    echo "lldb exit status: $lldb_status"
    echo "--- interesting lines from the lldb log (full log: $lldb_log) ---"
    grep -a -E 'stop reason|SIGTRAP|SIGABRT|thread #|frame #|brk|__builtin_trap|cgpt|Grid_' "$lldb_log" | head -n 60 || true
    echo "--- last 400 lines of the lldb log ---"
    tail -n 400 "$lldb_log" || true
}

# Testing hook: `CI_DEBUG_SOURCE_ONLY=1 source scripts/ci-debug-mac-nix-gpt.sh
# <out-dir>` defines the functions above without running the collection (used
# by tmp/task31-stub-test.sh on Linux).
if [ "${CI_DEBUG_SOURCE_ONLY:-0}" = "1" ] ; then
    return 0
fi

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

copy_nix_build_step_logs
collect_kept_build_dirs
collect_crash_reports

debug_light="${CI_DEBUG_LIGHT:-0}"
{
    echo "collection mode: $([ "$debug_light" = "1" ] && echo light || echo full)"
    echo "QLAT_GPT_TESTS_FAILED=${QLAT_GPT_TESTS_FAILED:-<unset>}"
    echo "heavy nix env build / --debug-signals rerun / lldb repro: $([ "$debug_light" = "1" ] && echo skipped-light-mode || echo enabled)"
} > "$out/debug-mode.txt" 2>/dev/null || true

if [ "$debug_light" = "1" ] ; then
    echo
    echo "==================== light mode: skipping nix env build, --debug-signals rerun and lldb repro ===================="
    echo "no nix build was triggered (light mode only reads the recorded logs and kept build dirs)"
    echo "collected environment, nix build step logs, kept build dirs, nix-store logs and crash reports are under: $out"
    exit 0
fi

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

run_lldb_repro

echo
echo "==================== debug run finished ===================="
echo "collected logs are under: $out"
