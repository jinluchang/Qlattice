#!/usr/bin/env python3

"""
Build and install a local Qlattice Jupyter kernel with nix.\n
This replaces the former 'install-py-local-kernel-with-nix.sh' and
'install-py-local-kernel-with-nix-cuda.sh'.  Everything that used to be
configured through environment variables ('name', 'use_nom', 'kernel_prefix')
is a command line option now, and '--all-variants' builds every variant that
nix lists for this machine in a single nix-build call, like the old '-cuda'
script did.  The list comes from 'nixpkgs/q-pkgs.nix' (the 'qlat-name-list'
attribute, queried with 'nix-instantiate --eval'), so it agrees with the
attributes that actually exist: without a usable NVIDIA GPU the '-cu', '-cuda'
and '-cudasupport' names are not listed at all.  The query needs the nix cache,
so it fails with an actionable error instead of guessing.\n
The nix-build out link is installed into the current user's kernel directory,
'$JUPYTER_DATA_DIR/kernels/nix-build-py-local<variant>/result' (by default
'$HOME/.local/share/jupyter/kernels/...'), when that directory is writable.
That is what the old bash script did, and it keeps the nix store path alive
through the indirect gc root nix registers for the out link.  The repo link
'./result-py-local<variant>' points at that out link, so it stays usable.\n
When the user level kernel directory is not writable, the Jupyter kernel spec,
the nix cache, the nom state directory and the temporary directory of the nix
and ipykernel subprocesses ('$TMPDIR') fall back to '<repo>/tmp/', so that
nothing outside of the repo is written to.  The out link still lives under
'<repo>/tmp/jupyter' because nix registers its indirect gc root there, but the
repo link './result-py-local<variant>' points at the nix store path directly
instead of going through that out link.  Files that the subprocesses would put
in a temporary directory ('$TMPDIR', by default '/tmp') then go to
'<repo>/tmp/tmp', and the IPython config directory of the kernel install goes to
'<repo>/tmp/jupyter/ipython', because some runners discard '/tmp'.\n
The nixpkgs used for the build can be selected with '--version' (a
NixOS release or channel such as '26.05' or 'unstable', fetched from
channels.nixos.org) or with '--nixpkgs' (a local checkout or a rev-pinned
tarball URL).  The default is the '<nixpkgs>' of NIX_PATH, and '--version ""'
selects that same default explicitly.  26.05 is the tested release for the
kernel environment; 24.11 and 25.11 currently fail to evaluate 'qlat-jhub-env'
(nixfmt-tree / xvfb are gone), so use '--nixpkgs' with a patched tree for those.
"""

import argparse
import json
import os
import re
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
  ./nixpkgs/install-py-local-kernel-with-nix.py --version 26.05
  ./nixpkgs/install-py-local-kernel-with-nix.py --variant clang --version nixos-26.05
  ./nixpkgs/install-py-local-kernel-with-nix.py --version ""
  ./nixpkgs/install-py-local-kernel-with-nix.py --nixpkgs /path/to/nixpkgs

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

'--all-variants' builds the variants that 'nixpkgs/q-pkgs.nix' lists, queried
from nix with

  nix-instantiate --read-write-mode --eval --strict --json q-pkgs.nix \
      -A qlat-name-list <the same --version/--nixpkgs/--arg/--argstr arguments>

'qlat-name-list' holds exactly the suffixes of the 'pkgs<suffix>' attributes, so
the cuda variants are simply not listed when there is no usable NVIDIA GPU (the
nix side collapses them to the plain CPU build when 'ngpu' is "0") and the build
can never ask for an attribute that does not exist.  The query needs a writable
nix cache; if it cannot run (no 'nix-instantiate', a non-zero exit, invalid
output, or a nixpkgs that cannot be fetched) the script stops with an error that
shows the command, the exit code and the last stderr lines, and suggests passing
explicit '--variant' options -- it never guesses the variant set from '/dev'.
'QLAT_NGPU=<int>' is forwarded to nix as '--argstr ngpu <int>' (an explicit
'--argstr ngpu ...' still wins), which both the name list and the build see.  An
explicit '--variant' is never filtered.

nixpkgs:
  --version VERSION          build with the nixpkgs of the NixOS release or
                             channel VERSION ('26.05', 'nixos-25.11',
                             'unstable', ...), fetched from
                             https://channels.nixos.org/nixos-<VERSION>/nixexprs.tar.xz
  --nixpkgs PATH-OR-URL      build with a local nixpkgs checkout or a
                             rev-pinned tarball such as
                             https://github.com/NixOS/nixpkgs/archive/<rev>.tar.gz
  Without either option the '<nixpkgs>' of NIX_PATH is used; '--version ""'
  selects that same default explicitly.  26.05 is the tested release for the
  kernel environment; 24.11 and 25.11 currently fail to evaluate 'qlat-jhub-env'
  (nixfmt-tree / xvfb are gone), so use '--nixpkgs' with a patched tree for
  those.

arguments that are not recognised here are passed on to nix-build, after the
'--version' / '--nixpkgs' arguments; an explicit '--argstr version ...'
or '--arg nixpkgs ...' therefore still wins (nix uses the last value).

outputs:
  $JUPYTER_DATA_DIR/kernels/nix-build-py-local<variant>/result
                                          nix-build out link and its gc root,
                                          used when the user level kernel
                                          directory (by default
                                          '$HOME/.local/share/jupyter/kernels')
                                          is writable
  ./result-py-local<variant>              link to that out link, or directly to
                                          the nix store path when the out link
                                          is installed in './tmp/jupyter'
                                          because the user level kernel
                                          directory is not writable
  ./tmp/jupyter/share/jupyter/kernels/    out link and kernel spec, used when
                                          the user level kernel directory is
                                          not writable
  ./tmp/nix-cache, ./tmp/nix-state        nix and nom caches, used when the
                                          user level directories are not writable
  ./tmp/tmp                               '$TMPDIR' of the nix and ipykernel
                                          subprocesses, used when the user level
                                          directories are not writable

The out link and the kernel spec are installed for the current user when the
user level kernel directory ('$JUPYTER_DATA_DIR/kernels', by default
'$HOME/.local/share/jupyter/kernels') is writable.  Otherwise they are installed
in './tmp/jupyter' (or '--kernel-prefix') and it has to be added to
JUPYTER_PATH; the run prints the exact line, e.g.:
  export JUPYTER_PATH="{REPO_PATH}/tmp/jupyter/share/jupyter${{JUPYTER_PATH:+:$JUPYTER_PATH}}"
"""

def parse_variant(value):
    name = value.strip()
    if name not in VARIANT_SUFFIX:
        raise argparse.ArgumentTypeError(
            f"unknown variant '{value}' (valid: {VARIANT_NAMES})"
        )
    return VARIANT_SUFFIX[name]

NIXPKGS_VERSION_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

def parse_nixpkgs_version(value):
    """Validate a NixOS release or channel; it is substituted into the channel URL.\n
    The literal empty string is valid and selects the default '<nixpkgs>' of NIX_PATH.
    """
    if value == "":
        return ""
    version = value.strip()
    if version.startswith("nixos-"):
        version = version[len("nixos-") :]
    if not NIXPKGS_VERSION_RE.match(version):
        raise argparse.ArgumentTypeError(
            f"expected a nixpkgs release or channel such as '26.05', "
            f"'nixos-25.11' or 'unstable', got '{value}'"
        )
    return version

def parse_nixpkgs_target(value):
    """Return ('dir', absolute path) or ('url', url) for a nixpkgs source."""
    target = value.strip()
    if target.startswith(("http://", "https://")):
        return ("url", target)
    path = Path(target).expanduser()
    if path.is_dir():
        return ("dir", str(path.resolve()))
    raise argparse.ArgumentTypeError(
        f"expected an existing nixpkgs directory or an http(s) tarball URL, "
        f"got '{value}'"
    )

def nix_str(value):
    """Quote 'value' as a nix string literal."""
    escaped = (
        value.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("$", "\\$")
        .replace("\n", "\\n")
        .replace("\r", "\\r")
        .replace("\t", "\\t")
    )
    return f'"{escaped}"'

def build_nixpkgs_args(version, target):
    """nix-build arguments for '--version' / '--nixpkgs' (empty by default)."""
    if version is not None:
        return ["--argstr", "version", version]
    if target is None:
        return []
    kind, value = target
    if kind == "dir":
        return ["--argstr", "nixpkgs", value]
    return ["--arg", "nixpkgs", f"(builtins.fetchTarball {nix_str(value)})"]

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
        help="build and install every variant that 'nix qlat-name-list' lists "
        "for this machine (the cuda variants are missing without a usable GPU; "
        "'QLAT_NGPU' is forwarded to nix as '--argstr ngpu')",
    )
    parser.add_argument(
        "--kernel-prefix",
        type=Path,
        metavar="DIR",
        help="install the repo local out link and kernel spec into "
        "DIR/share/jupyter/kernels "
        f"(default: {REPO_PATH / 'tmp' / 'jupyter'}); only used when they "
        "cannot be installed for the current user",
    )
    parser.add_argument(
        "--version",
        type=parse_nixpkgs_version,
        metavar="VERSION",
        help="nixpkgs NixOS release or channel to build with, e.g. '26.05' or "
        "'nixos-unstable'; fetches "
        "https://channels.nixos.org/nixos-<VERSION>/nixexprs.tar.xz "
        "(default: the '<nixpkgs>' of NIX_PATH; an empty value selects that "
        "same default)",
    )
    parser.add_argument(
        "--nixpkgs",
        type=parse_nixpkgs_target,
        metavar="PATH-OR-URL",
        help="nixpkgs source to build with: a local checkout directory or a "
        "rev-pinned tarball URL such as "
        "https://github.com/NixOS/nixpkgs/archive/<rev>.tar.gz "
        "(default: the '<nixpkgs>' of NIX_PATH)",
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
    if args.version is not None and args.nixpkgs is not None:
        parser.error("--nixpkgs cannot be combined with --version")
    args.nixpkgs_args = build_nixpkgs_args(args.version, args.nixpkgs)
    #
    ngpu = ngpu_from_env()
    if ngpu is not None:
        # 'QLAT_NGPU' is only an input to nix; a later explicit '--argstr ngpu ...'
        # still wins, because nix uses the last value.
        args.nix_build_args = ["--argstr", "ngpu", str(ngpu)] + args.nix_build_args
    #
    if args.all_variants and args.variant:
        parser.error("--all-variants cannot be combined with --variant")
    if args.all_variants:
        # Resolved in 'main()' from nix's 'qlat-name-list', which needs a writable
        # nix cache (set up there, after 'parse_args()' has run).
        args.variants = []
        args.skipped_variants = []
    elif args.variant:
        args.variants = list(dict.fromkeys(args.variant))
        args.skipped_variants = []
    else:
        args.variants = [""]
        args.skipped_variants = []
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

def ngpu_from_env():
    """Return the 'QLAT_NGPU' override as an int, or None when it is not set.\n
    It is only an INPUT to nix: 'nixpkgs/options.nix' turns 'ngpu' into the
    'use-cuda*' options, and with 'ngpu' equal to "0" the '-cu', '-cuda' and
    '-cudasupport' names disappear from 'qlat-name-list' altogether.  The value is
    therefore forwarded to nix as '--argstr ngpu <int>' (and both the name list and
    the build see it), instead of being interpreted by this script.
    """
    override = os.environ.get("QLAT_NGPU")
    if override is None:
        return None
    value = override.strip()
    if not value.isdigit():
        die(f"QLAT_NGPU must be a non-negative integer, got '{override}'")
    return int(value)

def arg_pairs(argv):
    """Keep only the arguments that 'nix-instantiate --eval' accepts.\n
    Arguments that are not recognised by the script are passed on to nix-build, so
    they can also contain nix-build only flags ('--no-out-link', '-o',
    '--keep-going', '-j', '--cores'), which 'nix-instantiate' rejects.  Only the
    ones that can change the evaluated attributes are kept: the '--arg',
    '--argstr' and '--option' triples and the '-I' / '--include' pairs.
    """
    result = []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg in ("--arg", "--argstr", "--option") and index + 2 < len(argv):
            result.extend(argv[index : index + 3])
            index += 3
        elif arg in ("-I", "--include") and index + 1 < len(argv):
            result.extend(argv[index : index + 2])
            index += 2
        elif arg.startswith("-I") and len(arg) > 2:
            result.append(arg)
            index += 1
        elif arg.startswith("--include=") or arg.startswith("--argstr="):
            result.append(arg)
            index += 1
        else:
            index += 1
    return result

def nix_query_error(reason, cmd, result=None):
    """Report why the 'qlat-name-list' query failed and exit.\n
    There is no fallback to a '/dev/nvidia*' count: a guess can select variants
    whose nix attributes do not exist, which is exactly the failure this query
    exists to avoid.
    """
    lines = [
        f"cannot query nix for the available variants: {reason}.",
        f"command: {shlex.join(str(c) for c in cmd)}",
    ]
    if result is not None:
        lines.append(f"exit code: {result.returncode}")
        stderr = (result.stderr or "").strip()
        if stderr:
            # nix's own 'error: ...' lines say what went wrong; the end of the
            # stderr can be the tail of an HTML error page from a failed fetch.
            stderr_lines = stderr.splitlines()
            error_indexes = [
                index
                for index, line in enumerate(stderr_lines)
                if line.startswith("error:")
            ]
            if error_indexes:
                # nix wraps its message, so keep the continuation lines too.
                tail = stderr_lines[error_indexes[-1] :][:6]
            else:
                tail = stderr_lines[-5:]
            lines.append("stderr (last lines):")
            lines.extend(f"  {line}" for line in tail)
    lines.append(
        "Pass explicit '--variant' options instead of '--all-variants' "
        f"(valid: {VARIANT_NAMES}) to build without the query."
    )
    die("\n".join(lines))

def nix_listed_variant_suffixes(args, env):
    """Suffixes of the variants that 'nixpkgs/q-pkgs.nix' actually lists.\n
    'qlat-name-list' is computed by nix itself and holds exactly the suffixes of
    the 'pkgs<suffix>' attributes, because it applies the same rules that decide
    whether an attribute exists at all: without a usable GPU the '-cu', '-cuda'
    and '-cudasupport' names are simply not emitted, so 'pkgs-cu' and friends do
    not exist.  Asking nix is therefore authoritative and can never disagree with
    the build.\n
    The query needs a writable nix cache, so it is called from 'main()' after the
    cache fallbacks are set up, not from 'parse_args()'.  A failure is fatal: see
    'nix_query_error()'.
    """
    cmd = [
        "nix-instantiate",
        "--read-write-mode",
        "--eval",
        "--strict",
        "--json",
        str(NIXPKGS_FILE),
        "-A",
        "qlat-name-list",
    ]
    cmd += list(args.nixpkgs_args) + arg_pairs(args.nix_build_args)
    print(f"+ {shlex.join(str(c) for c in cmd)}", flush=True)
    try:
        result = subprocess.run(
            cmd, env=env, capture_output=True, text=True, timeout=300
        )
    except FileNotFoundError as e:
        nix_query_error(f"command not found: '{e.filename}'", cmd)
    except subprocess.TimeoutExpired:
        nix_query_error("the query timed out after 300 s", cmd)
    if result.returncode != 0:
        nix_query_error("'nix-instantiate' failed", cmd, result)
    try:
        names = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        nix_query_error(f"the output is not valid JSON ({e})", cmd, result)
    if not isinstance(names, list):
        nix_query_error("the output is not a JSON list", cmd, result)
    listed = [name for name in names if isinstance(name, str)]
    suffixes = [suffix for suffix in VARIANT_SUFFIX.values() if suffix in listed]
    if not suffixes:
        nix_query_error(
            "'qlat-name-list' contains none of the known variant names", cmd, result
        )
    return suffixes

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
    default '$HOME/.local/state').\n
    Return True when at least one of the two directories fell back into the repo.
    """
    used_repo_fallback = False
    user_cache_home = Path(env.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    if not can_write_dir(user_cache_home):
        repo_cache = REPO_PATH / "tmp" / "nix-cache"
        print(
            f"The user nix cache is not writable. Using 'XDG_CACHE_HOME={repo_cache}'."
        )
        env["XDG_CACHE_HOME"] = str(repo_cache)
        seed_nix_cache(user_cache_home / "nix", repo_cache)
        used_repo_fallback = True
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
        used_repo_fallback = True
    #
    return used_repo_fallback

def setup_tmpdir(env):
    """Point '$TMPDIR' into the repo so that the subprocesses write their temp files there.\n
    nix, nom and the ipykernel install write their temporary files to '$TMPDIR'
    (by default '/tmp').  When the user level directories are not writable the
    script keeps everything else in '<repo>/tmp', so this directory has to
    follow; otherwise the writes land in '/tmp', which some runners discard
    together with the process that made them.  A '$TMPDIR' that is already set
    is left alone.
    """
    if env.get("TMPDIR"):
        return
    tmp_dir = REPO_PATH / "tmp" / "tmp"
    try:
        tmp_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"WARNING: Could not create '{tmp_dir}': {e}", file=sys.stderr)
        return
    env["TMPDIR"] = str(tmp_dir)
    print(f"Using 'TMPDIR={tmp_dir}' for the subprocesses.")

def setup_kernel_install(kernel_prefix):
    """Return '(kernels_dir, use_user_install)' for the out link and kernel spec.\n
    The current user's kernel directory
    ('$JUPYTER_DATA_DIR/kernels', by default '$HOME/.local/share/jupyter/kernels')
    is used when it is writable; otherwise everything is installed in the repo,
    under '<kernel_prefix>/share/jupyter/kernels'.
    """
    user_kernels = get_user_jupyter_data() / "kernels"
    if can_write_dir(user_kernels):
        print(f"Installing for the current user in '{user_kernels}'.")
        return user_kernels, True
    print(f"The user Jupyter kernel directory '{user_kernels}' is not writable.")
    repo_kernels = kernel_prefix / "share" / "jupyter" / "kernels"
    print(f"Installing in the repo instead, in '{repo_kernels}'.")
    repo_kernels.mkdir(parents=True, exist_ok=True)
    return repo_kernels, False

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
    cmd += args.nixpkgs_args
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

def kernel_out_link(kernels_dir, kernel_name):
    """Return the nix-build out link for 'kernel_name' in its 'nix-build-*' directory.\n
    nix registers an indirect gc root for the out link, so keeping it inside the
    kernel directory (as the old bash script did) is what keeps the store path
    alive.
    """
    out_dir = kernels_dir / f"nix-build-{kernel_name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / "result"

def link_repo_out_link(repo_link, target):
    """Make '<repo>/result-<kernel>' point at 'target', replacing any old link."""
    if repo_link.is_symlink() or repo_link.exists():
        repo_link.unlink()
    repo_link.symlink_to(target)
    print(f"'{repo_link}' -> '{target}'")

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
        # '$HOME' is not writable, so IPython would put its config directory in a
        # temporary directory outside of the repo; keep it in the repo instead.
        ipython_dir = (kernel_prefix / "ipython").resolve()
        ipython_dir.mkdir(parents=True, exist_ok=True)
        env = dict(env, IPYTHONDIR=str(ipython_dir))
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

def install_one_variant(variant, args, env, kernels_dir, use_user_install, bulk_build):
    kernel_name = f"py-local{variant}"
    out_link = kernel_out_link(kernels_dir, kernel_name)
    repo_link = REPO_PATH / f"result-{kernel_name}"
    #
    print()
    print(f"Building 'pkgs{variant}.qlat-jhub-env' for kernel '{kernel_name}'.")
    remove_stale_out_link(out_link)
    use_nom = (not args.no_nom) and (not bulk_build)
    nix_build([f"pkgs{variant}.qlat-jhub-env"], out_link, args, env, use_nom=use_nom)
    if not (out_link / "bin" / "python3").exists():
        die(f"'{out_link}/bin/python3' not found.")
    store_path = Path(os.path.realpath(out_link))
    print(f"'{out_link}' -> '{store_path}'")
    # The out link is kept for its indirect gc root, but when it is installed in
    # the repo (because the user level kernel directory is not writable), the
    # repo link points at the store path directly, so that it does not depend on
    # '<kernel_prefix>/share/jupyter/kernels/nix-build-*/result'.
    repo_target = out_link if use_user_install else store_path
    link_repo_out_link(repo_link, repo_target)
    install_kernel(variant, out_link, args.kernel_prefix, use_user_install, env)

def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    args = parse_args(argv)
    #
    if not NIXPKGS_FILE.is_file():
        die(f"'{NIXPKGS_FILE}' not found.  This script belongs in 'nixpkgs/'.")
    #
    env = dict(os.environ)
    used_repo_cache = setup_cache_fallbacks(env)
    kernels_dir, use_user_install = setup_kernel_install(args.kernel_prefix)
    if used_repo_cache or not use_user_install:
        setup_tmpdir(env)
    #
    if args.all_variants:
        listed = nix_listed_variant_suffixes(args, env)
        args.variants = listed
        args.skipped_variants = [
            suffix for suffix in VARIANT_SUFFIX.values() if suffix not in listed
        ]
    #
    if args.skipped_variants:
        names = ", ".join(SUFFIX_VARIANT[variant] for variant in args.skipped_variants)
        print(f"Skipping variants not listed by 'nix qlat-name-list': {names}.")
    #
    bulk_build = len(args.variants) > 1
    if bulk_build:
        attrs = [f"pkgs{variant}.qlat-jhub-env" for variant in args.variants]
        names = ", ".join(SUFFIX_VARIANT[variant] for variant in args.variants)
        print(f"Building {len(attrs)} variants: {names}.")
        nix_build(attrs, None, args, env, use_nom=not args.no_nom)
    #
    for variant in args.variants:
        install_one_variant(
            variant, args, env, kernels_dir, use_user_install, bulk_build
        )
    print()
    print("Finished successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
