#!/usr/bin/env python3

"""
Build and install a local Qlattice Jupyter kernel with nix.\n
This replaces the former 'install-py-local-kernel-with-nix.sh' and
'install-py-local-kernel-with-nix-cuda.sh'.  Everything that used to be
configured through environment variables ('name', 'use_nom', 'kernel_prefix')
is a command line option now, and '--all-variants' builds every variant in a
single nix-build call, like the old '-cuda' script did.\n
Nothing outside of the repo is written to when the usual user level directories
are not writable: the nix-build out link always goes into the repo, and the
Jupyter kernel spec, the nix cache and the nom state directory fall back to
'<repo>/tmp/' as needed.
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# '<repo>/nixpkgs', with symlinks resolved like 'pwd -P' does.
SCRIPT_PATH = Path(__file__).resolve().parent
REPO_PATH = SCRIPT_PATH.parent
NIXPKGS_FILE = SCRIPT_PATH / "q-pkgs.nix"

# Build variants, keyed by the name used on the command line.  The suffix is
# used for the nix attribute ('pkgs<suffix>.qlat-jhub-env'), the out link
# ('result-py-local<suffix>') and the kernel name ('py-local<suffix>').
VARIANT_SUFFIX = {
    "default": "",
    "cu": "-cu",
    "cuda": "-cuda",
    "cudasupport": "-cudasupport",
    "clang": "-clang",
    "pypi": "-pypi",
}
VARIANT_NAMES = ", ".join(VARIANT_SUFFIX)
SUFFIX_VARIANT = {suffix: name for name, suffix in VARIANT_SUFFIX.items()}

DEFAULT_JOBS = 6
DEFAULT_CORES = 15

EPILOG = f"""
examples:
  ./nixpkgs/install-py-local-kernel-with-nix.py
  ./nixpkgs/install-py-local-kernel-with-nix.py --variant cuda
  ./nixpkgs/install-py-local-kernel-with-nix.py --variant clang --variant pypi
  ./nixpkgs/install-py-local-kernel-with-nix.py --all-variants
  ./nixpkgs/install-py-local-kernel-with-nix.py --no-nom -j 4 --cores 15

variants (exact names, without a leading '-'):
  default       Qlattice and Grid/GPT.  This is the default when '--variant' is omitted.
  cu            some CUDA utilities.  Qlattice and Grid/GPT are NOT compiled with CUDA.
  cuda          some CUDA utilities.  Qlattice and Grid/GPT are compiled with CUDA.
  cudasupport   full CUDA support when possible (JAX GPU platforms are left enabled).
  clang         built with clang.
  pypi          the latest PyPI version.

a variant builds the nix attribute 'pkgs<suffix>.qlat-jhub-env', creates the out
link './result-py-local<suffix>' and installs the kernel 'py-local<suffix>', where
the suffix is '' for 'default' and '-<name>' otherwise: '--variant cuda' builds
'pkgs-cuda.qlat-jhub-env' and installs the kernel 'py-local-cuda'.

arguments that are not recognised here are passed on to nix-build.

outputs (all inside the repo):
  ./result-py-local<variant>              nix-build out link
  ./tmp/jupyter/share/jupyter/kernels/    kernel spec, used when it cannot be
                                          installed for the current user
  ./tmp/nix-cache, ./tmp/nix-state        nix and nom caches, used when the
                                          user level directories are not writable

The kernel spec is installed for the current user when the user level kernel
directory ('$JUPYTER_DATA_DIR/kernels', by default
'$HOME/.local/share/jupyter/kernels') is writable.  Otherwise it is installed in
'./tmp/jupyter' (or '--kernel-prefix') and has to be added to JUPYTER_PATH; the
run prints the exact line, e.g.:
  export JUPYTER_PATH="{REPO_PATH}/tmp/jupyter/share/jupyter${{JUPYTER_PATH:+:$JUPYTER_PATH}}"
"""

def parse_variant(value):
    name = value.strip()
    if name not in VARIANT_SUFFIX:
        raise argparse.ArgumentTypeError(
            f"unknown variant '{value}' (valid: {VARIANT_NAMES})"
        )
    return VARIANT_SUFFIX[name]

def check_variant_args(argv):
    """Reject the old '-cu' style names before argparse sees them.\n
    argparse would only report 'expected one argument' for '--variant -cu', which
    does not explain that the leading '-' has to be dropped.
    """
    for index, arg in enumerate(argv):
        if arg != "--variant" or index + 1 >= len(argv):
            continue
        value = argv[index + 1]
        if value.startswith("-") and not value.startswith("--"):
            die(
                f"'--variant' does not take '{value}': variant names are bare "
                f"names, without a leading '-' (valid: {VARIANT_NAMES})",
                2,
            )

def die(message, code=1):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(code)

def parse_args(argv):
    argv = sys.argv[1:] if argv is None else list(argv)
    check_variant_args(argv)
    parser = argparse.ArgumentParser(
        description="Build and install a local Qlattice Jupyter kernel with nix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    parser.add_argument(
        "--variant",
        action="append",
        type=parse_variant,
        metavar="NAME",
        help=f"build variant, may be repeated (default: default; valid: {VARIANT_NAMES})",
    )
    parser.add_argument(
        "--all-variants",
        action="store_true",
        help="build and install every variant",
    )
    parser.add_argument(
        "--kernel-prefix",
        type=Path,
        metavar="DIR",
        help="install the repo local kernel spec into DIR/share/jupyter/kernels "
        f"(default: {REPO_PATH / 'tmp' / 'jupyter'}); only used when the kernel "
        "spec cannot be installed for the current user",
    )
    parser.add_argument(
        "--no-nom",
        action="store_true",
        help="do not monitor nix-build with nom",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=DEFAULT_JOBS,
        help=f"nix-build -j (default: {DEFAULT_JOBS})",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=DEFAULT_CORES,
        help=f"nix-build --cores (default: {DEFAULT_CORES})",
    )
    args, nix_build_args = parser.parse_known_args(argv)
    args.nix_build_args = nix_build_args
    #
    if args.all_variants and args.variant:
        parser.error("--all-variants cannot be combined with --variant")
    if args.all_variants:
        args.variants = list(VARIANT_SUFFIX.values())
    elif args.variant:
        args.variants = list(dict.fromkeys(args.variant))
    else:
        args.variants = [""]
    #
    if args.kernel_prefix is None:
        args.kernel_prefix = REPO_PATH / "tmp" / "jupyter"
    return args

def run(cmd, env):
    print(f"+ {shlex.join(str(c) for c in cmd)}", flush=True)
    try:
        return subprocess.run(cmd, env=env, check=True)
    except FileNotFoundError as e:
        die(f"command not found: '{e.filename}'", 127)
    except subprocess.CalledProcessError as e:
        die(f"'{e.cmd[0]}' failed with exit code {e.returncode}", e.returncode)

def run_with_nom(cmd, env):
    """Run 'cmd' with its output piped into 'nom --json', like `cmd |& nom --json`."""
    print(f"+ {shlex.join(str(c) for c in cmd)} |& nom --json", flush=True)
    try:
        nix_build = subprocess.Popen(
            cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        nom = subprocess.Popen(["nom", "--json"], env=env, stdin=nix_build.stdout)
    except FileNotFoundError as e:
        die(f"command not found: '{e.filename}'", 127)
    nix_build.stdout.close()  # let nix-build see the pipe close when nom exits
    nom_returncode = nom.wait()
    nix_build_returncode = nix_build.wait()
    # Unlike a shell pipeline, both exit codes are checked here.
    if nix_build_returncode != 0:
        die(
            f"'nix-build' failed with exit code {nix_build_returncode}",
            nix_build_returncode,
        )
    if nom_returncode != 0:
        die(f"'nom' failed with exit code {nom_returncode}", nom_returncode)

def can_write_dir(path):
    """Return True when a file can be created in 'path', creating it if needed.\n
    A read-only mount is not reported by 'os.access(path, os.W_OK)', so the
    check has to actually write something.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError:
        return False
    probe = path / f".write-probe-{os.getpid()}"
    try:
        with open(probe, "w"):
            pass
    except OSError:
        return False
    try:
        probe.unlink()
    except OSError:
        pass
    return True

def get_user_jupyter_data():
    """The '--user' data directory, as computed by 'jupyter_core.paths.jupyter_data_dir'."""
    if os.environ.get("JUPYTER_DATA_DIR"):
        return Path(os.environ["JUPYTER_DATA_DIR"])
    xdg_data_home = os.environ.get("XDG_DATA_HOME") or str(
        Path.home() / ".local" / "share"
    )
    return Path(xdg_data_home) / "jupyter"

def seed_nix_cache(user_nix_cache, repo_cache):
    """Copy the readable user git cache into the repo cache to avoid re-downloads."""
    src = user_nix_cache / "gitv3"
    dst = repo_cache / "nix" / "gitv3"
    if dst.exists() or not src.is_dir():
        return
    print(f"Seeding the repo local nix cache from '{src}'.", flush=True)
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst, symlinks=True)
    except OSError as e:
        print(f"WARNING: Could not seed the repo local nix cache: {e}", file=sys.stderr)
        print("WARNING: It may be re-downloaded.", file=sys.stderr)
        shutil.rmtree(dst, ignore_errors=True)

def setup_cache_fallbacks(env):
    """Point the nix and nom cache directories into the repo when they are not writable.\n
    nix keeps its evaluation and git cache in '$XDG_CACHE_HOME/nix' (by default
    '$HOME/.cache/nix') and nom keeps its build reports in '$XDG_STATE_HOME' (by
    default '$HOME/.local/state').
    """
    user_cache_home = Path(env.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    if not can_write_dir(user_cache_home):
        repo_cache = REPO_PATH / "tmp" / "nix-cache"
        print(
            f"The user nix cache is not writable. Using 'XDG_CACHE_HOME={repo_cache}'."
        )
        env["XDG_CACHE_HOME"] = str(repo_cache)
        seed_nix_cache(user_cache_home / "nix", repo_cache)
    #
    user_state_home = Path(
        env.get("XDG_STATE_HOME") or Path.home() / ".local" / "state"
    )
    if not can_write_dir(user_state_home):
        repo_state = REPO_PATH / "tmp" / "nix-state"
        print(
            f"The user state directory is not writable. Using 'XDG_STATE_HOME={repo_state}'."
        )
        env["XDG_STATE_HOME"] = str(repo_state)

def setup_kernel_install(kernel_prefix):
    """Return True when the kernel spec is to be installed for the current user."""
    user_kernels = get_user_jupyter_data() / "kernels"
    if can_write_dir(user_kernels):
        print(
            f"Installing the Jupyter kernel spec for the current user in '{user_kernels}'."
        )
        return True
    print(f"The user Jupyter kernel directory '{user_kernels}' is not writable.")
    print(
        f"Installing the Jupyter kernel spec in the repo instead, in '{kernel_prefix}'."
    )
    kernel_prefix.mkdir(parents=True, exist_ok=True)
    return False

def nix_build(attrs, out_link, args, env, use_nom):
    cmd = ["nix-build", str(NIXPKGS_FILE)]
    for attr in attrs:
        cmd += ["-A", attr]
    if out_link is None:
        cmd.append("--no-out-link")
    else:
        cmd += ["-o", str(out_link)]
    if use_nom:
        cmd += ["--log-format", "internal-json", "-v"]
    cmd += ["-j", str(args.jobs), "--cores", str(args.cores)]
    cmd += args.nix_build_args
    #
    if use_nom:
        run_with_nom(cmd, env)
    else:
        run(cmd, env)

def remove_stale_out_link(out_link):
    """Remove an out link that nix-build would refuse to replace.\n
    nix-build fails with 'cannot create symlink ...; already exists' when the
    out link points somewhere outside of the nix store, which is what older
    versions of this script left behind.
    """
    if out_link.is_symlink():
        print(f"Removing stale out link '{out_link}'.")
        out_link.unlink()

def ipykernel_env_args(variant, out_link):
    args = [
        "--env",
        "SHELL",
        str(out_link / "bin" / "bash"),
        "--env",
        "LOCALE_ARCHIVE",
        "/run/current-system/sw/lib/locale/locale-archive",
        "--env",
        "PATH",
        f"{out_link}/bin:/run/current-system/sw/bin",
        "--env",
        "PKG_CONFIG_PATH",
        str(out_link / "lib" / "pkgconfig"),
        "--env",
        "LD_LIBRARY_PATH",
        f"/run/opengl-driver/lib:{out_link}/lib",
        "--env",
        "LIBRARY_PATH",
        str(out_link / "lib"),
        "--env",
        "CPATH",
        str(out_link / "include"),
        "--env",
        "PYTHONPATH",
        "",
        "--env",
        "CUBACORES",
        "0",
        "--env",
        "OMP_NUM_THREADS",
        "2",
        "--env",
        "JAX_ENABLE_X64",
        "True",
    ]
    if variant != "-cudasupport":
        args += ["--env", "JAX_PLATFORMS", "cpu"]
    return args

def install_kernel(variant, out_link, kernel_prefix, use_user_install, env):
    kernel_name = f"py-local{variant}"
    cmd = [str(out_link / "bin" / "python3"), "-m", "ipykernel", "install"]
    if use_user_install:
        cmd.append("--user")
    else:
        cmd += ["--prefix", str(kernel_prefix)]
    cmd += ipykernel_env_args(variant, out_link)
    cmd.append(f"--name={kernel_name}")
    run(cmd, env)
    #
    print()
    if use_user_install:
        print(f"Kernel '{kernel_name}' installed for the current user.")
    else:
        jupyter_data_dir = kernel_prefix / "share" / "jupyter"
        print(
            f"Kernel spec installed in '{jupyter_data_dir / 'kernels' / kernel_name}'."
        )
        print("To use it, add it to JUPYTER_PATH, e.g.:")
        print(
            f'    export JUPYTER_PATH="{jupyter_data_dir}${{JUPYTER_PATH:+:$JUPYTER_PATH}}"'
        )

def install_one_variant(variant, args, env, use_user_install, bulk_build):
    kernel_name = f"py-local{variant}"
    out_link = REPO_PATH / f"result-{kernel_name}"
    #
    print()
    print(f"Building 'pkgs{variant}.qlat-jhub-env' for kernel '{kernel_name}'.")
    remove_stale_out_link(out_link)
    use_nom = (not args.no_nom) and (not bulk_build)
    nix_build([f"pkgs{variant}.qlat-jhub-env"], out_link, args, env, use_nom=use_nom)
    if not (out_link / "bin" / "python3").exists():
        die(f"'{out_link}/bin/python3' not found.")
    print(f"'{out_link}' -> '{os.path.realpath(out_link)}'")
    install_kernel(variant, out_link, args.kernel_prefix, use_user_install, env)

def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    args = parse_args(argv)
    #
    if not NIXPKGS_FILE.is_file():
        die(f"'{NIXPKGS_FILE}' not found.  This script belongs in 'nixpkgs/'.")
    #
    env = dict(os.environ)
    setup_cache_fallbacks(env)
    use_user_install = setup_kernel_install(args.kernel_prefix)
    #
    bulk_build = len(args.variants) > 1
    if bulk_build:
        attrs = [f"pkgs{variant}.qlat-jhub-env" for variant in args.variants]
        names = ", ".join(SUFFIX_VARIANT[variant] for variant in args.variants)
        print(f"Building {len(attrs)} variants: {names}.")
        nix_build(attrs, None, args, env, use_nom=not args.no_nom)
    #
    for variant in args.variants:
        install_one_variant(variant, args, env, use_user_install, bulk_build)
    print()
    print("Finished successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
