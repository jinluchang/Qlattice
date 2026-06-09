#!/usr/bin/env python3

"""
Run a single qlat Python example test using a pre-built nix environment.

Usage:
  ./run-one-example-py.py <test-name> [options]

Options:
  --cuda          Use CUDA-enabled build (result-py-local-cuda)
  --cudasupport   Use CUDA support build (result-py-local-cudasupport)
  --cu            Use CUDA utilities build (result-py-local-cu)
  --clang         Use clang build (result-py-local-clang)
  --pypi          Use PyPI build (result-py-local-pypi)
  --n <num>       Number of MPI processes (default: 2)
  --timeout <dur> Timeout duration (default: 60m)
  --help          Show this help message

Examples:
  ./run-one-example-py.py utils
  ./run-one-example-py.py auto-contract-01 --cuda
  ./run-one-example-py.py qtopo-measure --cudasupport --n 4
"""

import argparse
import os
import shutil
import subprocess
import sys


def show_help():
    print(__doc__.strip())
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("test_name", nargs="?", default="")
    parser.add_argument("--cuda", action="store_true")
    parser.add_argument("--cudasupport", action="store_true")
    parser.add_argument("--cu", action="store_true")
    parser.add_argument("--clang", action="store_true")
    parser.add_argument("--pypi", action="store_true")
    parser.add_argument("--n", type=int, default=2)
    parser.add_argument("--timeout", default="60m")
    parser.add_argument("--help", "-h", action="store_true")
    args = parser.parse_args()

    if args.help:
        show_help()

    if not args.test_name:
        script_path = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_path)
        examples_dir = os.path.join(project_root, "examples-py")
        tests = sorted(
            f[:-3] for f in os.listdir(examples_dir)
            if f.endswith(".py")
        )
        print(f"Usage: {sys.argv[0]} <test-name> [--cuda|--cudasupport|--cu|--clang|--pypi]")
        print()
        print("Available tests:")
        for t in tests:
            print(f"  {t}")
        sys.exit(0)

    build_variants = []
    if args.cuda:
        build_variants.append("-cuda")
    if args.cudasupport:
        build_variants.append("-cudasupport")
    if args.cu:
        build_variants.append("-cu")
    if args.clang:
        build_variants.append("-clang")
    if args.pypi:
        build_variants.append("-pypi")

    if len(build_variants) > 1:
        print("Error: only one of --cuda, --cudasupport, --cu, --clang, --pypi may be specified", file=sys.stderr)
        sys.exit(1)

    build_variant = build_variants[0] if build_variants else ""

    script_path = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_path)
    result_dir = os.path.join(project_root, f"result-py-local{build_variant}")

    if not os.path.isdir(result_dir):
        print(f"Error: result directory not found: {result_dir}", file=sys.stderr)
        print("Run nix-build first to create the environment.", file=sys.stderr)
        sys.exit(1)

    python3 = os.path.join(result_dir, "bin", "python3")
    if not os.path.isfile(python3) or not os.access(python3, os.X_OK):
        print(f"Error: python3 not found in {result_dir}", file=sys.stderr)
        sys.exit(1)

    # Set up environment via setenv-qlat.sh
    setenv_sh = os.path.join(result_dir, "bin", "setenv-qlat.sh")
    if not os.path.isfile(setenv_sh):
        print(f"Error: setenv-qlat.sh not found: {setenv_sh}", file=sys.stderr)
        sys.exit(1)

    result = subprocess.run(
        ["bash", "-c", f"source {setenv_sh} && env"],
        capture_output=True, text=True
    )
    env = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            env[key] = value



    # Locate the test script
    test_script = os.path.join(project_root, "examples-py", f"{args.test_name}.py")
    if not os.path.isfile(test_script):
        print(f"Error: test script not found: {test_script}", file=sys.stderr)
        sys.exit(1)

    # Set up mpi options
    mpi_options = ["--oversubscribe", "--bind-to", "none"]
    if "cuda" in build_variant:
        mpi_options.extend(["--mca", "pml", "^ucx"])

    # Run the test
    work_dir = os.path.join(project_root, "tmp", "examples-py", f"{args.test_name}.py.p")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    print(f"Running test: {args.test_name}")
    print(f"Build variant: {build_variant if build_variant else 'default (CPU)'}")
    print(f"MPI processes: {args.n}")
    print(f"Working directory: {work_dir}")
    print()

    shutil.copy2(test_script, work_dir)
    json_file = f"{args.test_name}.log.json"
    json_src = os.path.join(project_root, "examples-py", json_file)
    if os.path.isfile(json_src):
        shutil.copy2(json_src, work_dir)

    test_script_in_dir = os.path.join(work_dir, os.path.basename(test_script))
    log_full = os.path.join(work_dir, "log.full.txt")

    cmd = [
        "timeout", "-s", "KILL", args.timeout,
        "mpiexec", "-n", str(args.n),
        *mpi_options,
        python3, "-m", "mpi4py", test_script_in_dir, "--test",
        "-qmp-geom", "1", "1", "1", "2",
        "--mpi", "1.1.1.2", "--mpi_split", "1.1.1.1",
        "--mpi", "1.1.2", "--mpi_split", "1.1.1",
    ]

    env["q_verbose"] = "1"

    with open(log_full, "w") as f:
        subprocess.run(cmd, cwd=work_dir, env=env, stdout=f, stderr=subprocess.STDOUT)

    # Extract and show CHECK lines
    print("--- CHECK output ---")
    with open(log_full, "rb") as f:
        for line in f:
            if b"CHECK: " in line:
                print(line.decode("utf-8", errors="replace").rstrip())
    print("--- end CHECK output ---")

    # Compare with reference log if it exists
    ref_log = os.path.join(project_root, "examples-py", f"{args.test_name}.log")
    if os.path.isfile(ref_log):
        with open(log_full, "rb") as f:
            new_checks = [l for l in f if b"CHECK: " in l]
        with open(ref_log, "rb") as f:
            ref_checks = [l for l in f if b"CHECK: " in l]

        if new_checks == ref_checks:
            print()
            print(f"PASSED: {args.test_name}")
        else:
            print()
            print(f"FAILED: {args.test_name} (output differs from reference)")
            print("Last 50 lines of output:")
            with open(log_full, "rb") as f:
                lines = f.readlines()
                for line in lines[-50:]:
                    print(line.decode("utf-8", errors="replace").rstrip())
            sys.exit(1)
    else:
        print()
        print(f"WARNING: no reference log found at {ref_log}")
        print("Test ran but cannot verify output.")


if __name__ == "__main__":
    main()
