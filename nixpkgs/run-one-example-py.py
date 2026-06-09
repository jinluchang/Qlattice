#!/usr/bin/env python3

"""
Run a single qlat Python example test using a pre-built nix environment.\n
Usage:
  ./run-one-example-py.py <test-name> [options]\n
Options:
  --cuda          Use CUDA-enabled build (result-py-local-cuda)
  --cudasupport   Use CUDA support build (result-py-local-cudasupport)
  --cu            Use CUDA utilities build (result-py-local-cu)
  --clang         Use clang build (result-py-local-clang)
  --pypi          Use PyPI build (result-py-local-pypi)
  --n <num>       Number of MPI processes (default: 2)
  --timeout <dur> Timeout duration (default: 60m)
  --help          Show this help message\n
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

def parse_args():
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
    return parser.parse_args()

def get_project_root():
    script_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(script_path)

def list_available_tests(project_root):
    examples_dir = os.path.join(project_root, "examples-py")
    tests = sorted(f[:-3] for f in os.listdir(examples_dir) if f.endswith(".py"))
    print(
        f"Usage: {sys.argv[0]} <test-name> [--cuda|--cudasupport|--cu|--clang|--pypi]"
    )
    print()
    print("Available tests:")
    for t in tests:
        print(f"  {t}")
    sys.exit(0)

def resolve_build_variant(args):
    variants = []
    if args.cuda:
        variants.append("-cuda")
    if args.cudasupport:
        variants.append("-cudasupport")
    if args.cu:
        variants.append("-cu")
    if args.clang:
        variants.append("-clang")
    if args.pypi:
        variants.append("-pypi")
    #
    if len(variants) > 1:
        print(
            "Error: only one of --cuda, --cudasupport, --cu, --clang, --pypi may be specified",
            file=sys.stderr,
        )
        sys.exit(1)
    #
    return variants[0] if variants else ""

def validate_result_dir(result_dir):
    if not os.path.isdir(result_dir):
        print(f"Error: result directory not found: {result_dir}", file=sys.stderr)
        print("Run nix-build first to create the environment.", file=sys.stderr)
        sys.exit(1)

def validate_python3(result_dir):
    python3 = os.path.join(result_dir, "bin", "python3")
    if not os.path.isfile(python3) or not os.access(python3, os.X_OK):
        print(f"Error: python3 not found in {result_dir}", file=sys.stderr)
        sys.exit(1)
    return python3

def load_environment(result_dir):
    setenv_sh = os.path.join(result_dir, "bin", "setenv-qlat.sh")
    if not os.path.isfile(setenv_sh):
        print(f"Error: setenv-qlat.sh not found: {setenv_sh}", file=sys.stderr)
        sys.exit(1)
    #
    result = subprocess.run(
        ["bash", "-c", f"source {setenv_sh} && env"], capture_output=True, text=True
    )
    env = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            env[key] = value
    return env

def validate_test_script(project_root, test_name):
    test_script = os.path.join(project_root, "examples-py", f"{test_name}.py")
    if not os.path.isfile(test_script):
        print(f"Error: test script not found: {test_script}", file=sys.stderr)
        sys.exit(1)
    return test_script

def get_mpi_options(build_variant):
    mpi_options = ["--oversubscribe", "--bind-to", "none"]
    return mpi_options

def prepare_work_dir(project_root, test_name, test_script):
    work_dir = os.path.join(project_root, "tmp", "examples-py", f"{test_name}.py.p")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    #
    shutil.copy2(test_script, work_dir)
    json_file = f"{test_name}.log.json"
    json_src = os.path.join(project_root, "examples-py", json_file)
    if os.path.isfile(json_src):
        shutil.copy2(json_src, work_dir)
    #
    return work_dir

def build_command(python3, work_dir, test_name, n_procs, timeout, mpi_options):
    test_script_in_dir = os.path.join(work_dir, f"{test_name}.py")
    return [
        "timeout",
        "-s",
        "KILL",
        timeout,
        "mpiexec",
        "-n",
        str(n_procs),
        *mpi_options,
        python3,
        "-m",
        "mpi4py",
        test_script_in_dir,
        "--test",
        "-qmp-geom",
        "1",
        "1",
        "1",
        "2",
        "--mpi",
        "1.1.1.2",
        "--mpi_split",
        "1.1.1.1",
        "--mpi",
        "1.1.2",
        "--mpi_split",
        "1.1.1",
    ]

def run_test(cmd, work_dir, env):
    env["q_verbose"] = "1"
    log_full = os.path.join(work_dir, "log.full.txt")
    with open(log_full, "w") as f:
        subprocess.run(cmd, cwd=work_dir, env=env, stdout=f, stderr=subprocess.STDOUT)
    return log_full

def show_check_output(log_full):
    print("--- CHECK output ---")
    with open(log_full, "rb") as f:
        for line in f:
            if b"CHECK: " in line:
                print(line.decode("utf-8", errors="replace").rstrip())
    print("--- end CHECK output ---")

def extract_checks(filepath):
    with open(filepath, "rb") as f:
        return [l for l in f if b"CHECK: " in l]

def compare_with_reference(project_root, test_name, log_full):
    ref_log = os.path.join(project_root, "examples-py", f"{test_name}.log")
    if not os.path.isfile(ref_log):
        print()
        print(f"WARNING: no reference log found at {ref_log}")
        print("Test ran but cannot verify output.")
        return
    #
    new_checks = extract_checks(log_full)
    ref_checks = extract_checks(ref_log)
    #
    if new_checks == ref_checks:
        print()
        print(f"PASSED: {test_name}")
    else:
        print()
        print(f"FAILED: {test_name} (output differs from reference)")
        print("Last 50 lines of output:")
        with open(log_full, "rb") as f:
            lines = f.readlines()
            for line in lines[-50:]:
                print(line.decode("utf-8", errors="replace").rstrip())
        sys.exit(1)

def main():
    args = parse_args()
    #
    if args.help:
        show_help()
    #
    project_root = get_project_root()
    #
    if not args.test_name:
        list_available_tests(project_root)
    #
    build_variant = resolve_build_variant(args)
    result_dir = os.path.join(project_root, f"result-py-local{build_variant}")
    #
    validate_result_dir(result_dir)
    python3 = validate_python3(result_dir)
    env = load_environment(result_dir)
    test_script = validate_test_script(project_root, args.test_name)
    #
    work_dir = prepare_work_dir(project_root, args.test_name, test_script)
    mpi_options = get_mpi_options(build_variant)
    cmd = build_command(
        python3, work_dir, args.test_name, args.n, args.timeout, mpi_options
    )
    #
    print(f"Running test: {args.test_name}")
    print(f"Build variant: {build_variant if build_variant else 'default (CPU)'}")
    print(f"MPI processes: {args.n}")
    print(f"Working directory: {work_dir}")
    print()
    #
    log_full = run_test(cmd, work_dir, env)
    show_check_output(log_full)
    compare_with_reference(project_root, args.test_name, log_full)

if __name__ == "__main__":
    main()
