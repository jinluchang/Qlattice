#!/usr/bin/env bash

# Run a single qlat Python example test using a pre-built nix environment.
#
# Usage:
#   ./run-one-example-py.sh <test-name> [options]
#
# Options:
#   --cuda          Use CUDA-enabled build (result-py-local-cuda)
#   --cudasupport   Use CUDA support build (result-py-local-cudasupport)
#   --cu            Use CUDA utilities build (result-py-local-cu)
#   --clang         Use clang build (result-py-local-clang)
#   --n <num>       Number of MPI processes (default: 2)
#   --timeout <dur> Timeout duration (default: 60m)
#   --help          Show this help message
#
# Examples:
#   ./run-one-example-py.sh utils
#   ./run-one-example-py.sh auto-contract-01 --cuda
#   ./run-one-example-py.sh qtopo-measure --cudasupport --n 4

set -euo pipefail

script_path="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 ; pwd -P )"

# Default options
nprocs=2
timeout_dur="60m"
build_variant=""
test_name=""

show_help() {
    sed -n '/^# Usage:/,/^$/p' "$0" | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cuda)
            build_variant="-cuda"
            shift
            ;;
        --cudasupport)
            build_variant="-cudasupport"
            shift
            ;;
        --cu)
            build_variant="-cu"
            shift
            ;;
        --clang)
            build_variant="-clang"
            shift
            ;;
        --n)
            nprocs="$2"
            shift 2
            ;;
        --timeout)
            timeout_dur="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            if [ -z "$test_name" ]; then
                test_name="$1"
            else
                echo "Unexpected argument: $1" >&2
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$test_name" ]; then
    echo "Error: test name required" >&2
    echo "Usage: $0 <test-name> [--cuda|--cudasupport|--cu|--clang]" >&2
    exit 1
fi

# Locate the result directory (in project root, next to nixpkgs/)
project_root="$(dirname "$script_path")"
result_dir="$project_root/result-py-local${build_variant}"
if [ ! -d "$result_dir" ]; then
    echo "Error: result directory not found: $result_dir" >&2
    echo "Run nix-build first to create the environment." >&2
    exit 1
fi

python3="$result_dir/bin/python3"
if [ ! -x "$python3" ]; then
    echo "Error: python3 not found in $result_dir" >&2
    exit 1
fi

# Set up environment
export SHELL="$result_dir/bin/bash"
export PATH="$result_dir/bin:/run/current-system/sw/bin:$PATH"
export PKG_CONFIG_PATH="$result_dir/lib/pkgconfig:$result_dir/share/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="/run/opengl-driver/lib:$result_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LIBRARY_PATH="$result_dir/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export PYTHONPATH=""
export CUBACORES=0
export OMP_NUM_THREADS=2
export JAX_ENABLE_X64=True

# CPATH: CUDA headers conflict with glibc headers (noexcept mismatch)
# For CUDA builds, clear CPATH and remove conflicting glibc math headers
if [[ "$build_variant" == *cuda* ]]; then
    export CPATH=""
    # Remove glibc math headers that conflict with CUDA
    if [ -d "$result_dir/include/bits" ]; then
        rm -f "$result_dir/include/bits/mathcalls.h" "$result_dir/include/bits/mathcalls-macros.h" 2>/dev/null || true
    fi
fi
export CPATH="$result_dir/include${CPATH:+:$CPATH}"

if [ "$build_variant" != "-cudasupport" ]; then
    export JAX_PLATFORMS=cpu
fi

# Source CUDA compilation environment if using CUDA variant
cuda_mpi_qlat_sh="$result_dir/bin/cuda-mpi-qlat.sh"
if [[ "$build_variant" == *cuda* ]] && [ -e "$cuda_mpi_qlat_sh" ]; then
    # The symlinks may point to binary wrappers (ELF). Extract the real shell script.
    real_cuda_mpi_qlat_sh=$(strings "$cuda_mpi_qlat_sh" 2>/dev/null | grep -m1 "cuda-mpi-qlat.sh$" || true)
    if [ -z "$real_cuda_mpi_qlat_sh" ]; then
        real_cuda_mpi_qlat_sh="$cuda_mpi_qlat_sh"
    fi
    echo "Sourcing CUDA compilation environment from: $real_cuda_mpi_qlat_sh"
    source "$real_cuda_mpi_qlat_sh" echo
    # Re-export PKG_CONFIG_PATH after sourcing (it may override env)
    export PKG_CONFIG_PATH="$result_dir/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    echo "CXX=${CXX:-unset}"
    echo "MPICXX=${MPICXX:-unset}"
    echo "OMPI_CXX=${OMPI_CXX:-unset}"
    echo "PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-unset}"
fi

# Locate the test script
test_script="$project_root/examples-py/${test_name}.py"
if [ ! -f "$test_script" ]; then
    echo "Error: test script not found: $test_script" >&2
    exit 1
fi

# Set up mpi options
mpi_options="--oversubscribe --bind-to none"
if [[ "$build_variant" == *cuda* ]]; then
    mpi_options="$mpi_options --mca pml ^ucx"
fi

# Run the test
work_dir="$project_root/examples-py/${test_name}.py.p"
rm -rf "$work_dir"
mkdir -p "$work_dir"

echo "Running test: $test_name"
echo "Build variant: ${build_variant:-default (CPU)}"
echo "MPI processes: $nprocs"
echo "Working directory: $work_dir"
echo

cd "$work_dir"
cp "$test_script" .
json_file="${test_name}.log.json"
if [ -f "$project_root/examples-py/$json_file" ]; then
    cp "$project_root/examples-py/$json_file" .
fi

q_verbose=1 timeout -s KILL "$timeout_dur" \
    mpiexec -n "$nprocs" $mpi_options \
    "$python3" -m mpi4py "./${test_name}.py" --test \
    -qmp-geom 1 1 1 2 \
    --mpi 1.1.1.2 --mpi_split 1.1.1.1 \
    --mpi 1.1.2 --mpi_split 1.1.1 \
    > log.full.txt 2>&1 || true

# Extract and show CHECK lines
echo "--- CHECK output ---"
grep -a "CHECK: " log.full.txt || true
echo "--- end CHECK output ---"

# Compare with reference log if it exists
ref_log="$project_root/examples-py/${test_name}.log"
if [ -f "$ref_log" ]; then
    grep -a "CHECK: " log.full.txt > log.check.txt.new
    grep -a "CHECK: " "$ref_log" > log.check.txt
    if diff log.check.txt log.check.txt.new; then
        echo
        echo "PASSED: $test_name"
    else
        echo
        echo "FAILED: $test_name (output differs from reference)"
        echo "Last 50 lines of output:"
        tail -n 50 log.full.txt
        exit 1
    fi
else
    echo
    echo "WARNING: no reference log found at $ref_log"
    echo "Test ran but cannot verify output."
fi
