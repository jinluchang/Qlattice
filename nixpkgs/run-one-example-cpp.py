#!/usr/bin/env python3

"""
Run a single qlat C++ example test using a pre-built nix environment.

This script sets up the environment from a qlat build produced by nix and then
runs the test via make with the examples-cpp/Makefile.  The Makefile itself is
standalone — it runs tests using whatever qlat is available in the current
environment.

The nix build is created by nixpkgs/install-py-local-kernel-with-nix.sh, which
runs nix-build and creates a ./result-py-local symlink (or ./result-py-local-*
for variant builds) pointing to the nix store path.  This script sources the
setenv-qlat.sh from that result directory to configure PATH, PYTHONPATH,
LD_LIBRARY_PATH, etc., then delegates to make.

The test is run inside ./tmp/examples-cpp/<test-name>/build/.  After the run,
check that directory for log.full (full output), log.txt (filtered), and
log.check.txt (CHECK: lines only).  The top-level
./tmp/examples-cpp/<test-name>/log and log.full are also available.

Usage:
  ./nixpkgs/run-one-example-cpp.py <test-name> [options]

Options:
  --cuda          Use CUDA-enabled build (result-py-local-cuda)
  --cudasupport   Use CUDA support build (result-py-local-cudasupport)
  --cu            Use CUDA utilities build (result-py-local-cu)
  --clang         Use clang build (result-py-local-clang)
  --pypi          Use PyPI build (result-py-local-pypi)
  --help          Show this help message

Examples:
  ./nixpkgs/run-one-example-cpp.py simple-1
  ./nixpkgs/run-one-example-cpp.py hmc --cuda
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
    parser.add_argument("--help", "-h", action="store_true")
    return parser.parse_args()

def get_project_root():
    script_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(script_path)

def list_available_tests(project_root):
    examples_dir = os.path.join(project_root, "examples-cpp")
    tests = sorted(f for f in os.listdir(examples_dir) if os.path.isdir(os.path.join(examples_dir, f)) and not f.startswith('.'))
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

def validate_test_dir(project_root, test_name):
    test_dir = os.path.join(project_root, "examples-cpp", test_name)
    if not os.path.isdir(test_dir):
        print(f"Error: test directory not found: {test_dir}", file=sys.stderr)
        sys.exit(1)
    return test_dir

def prepare_work_dir(project_root, test_name):
    work_dir = os.path.join(project_root, "tmp", "examples-cpp")
    os.makedirs(work_dir, exist_ok=True)
    #
    build_dir = os.path.join(work_dir, test_name, "build")
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    #
    src_dir = os.path.join(project_root, "examples-cpp")
    shutil.copytree(
        os.path.join(src_dir, test_name),
        os.path.join(work_dir, test_name),
        dirs_exist_ok=True,
    )
    shutil.copy2(os.path.join(src_dir, "Makefile"), work_dir)
    #
    return work_dir

def run_test(work_dir, test_name, env):
    cmd = ["make", "-f", "../Makefile", "-C", test_name, "run-proj"]
    print(f"Running: make -f ../Makefile -C {test_name} run-proj")
    print(f"In: {work_dir}")
    print()
    result = subprocess.run(cmd, cwd=work_dir, env=env)
    return result.returncode

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
    env = load_environment(result_dir)
    validate_test_dir(project_root, args.test_name)
    work_dir = prepare_work_dir(project_root, args.test_name)
    #
    print(f"Running test: {args.test_name}")
    print(f"Build variant: {build_variant if build_variant else 'default (CPU)'}")
    print(f"Work directory: {work_dir}")
    print()
    #
    rc = run_test(work_dir, args.test_name, env)
    p_dir = os.path.join(work_dir, args.test_name, "build")
    print()
    print("NOTE: This run did NOT modify the original source files")
    print(f"      in {os.path.join(project_root, 'examples-cpp')}/")
    print()
    print(f"Regenerated files are in: {p_dir}/")
    sys.exit(rc)

if __name__ == "__main__":
    main()
