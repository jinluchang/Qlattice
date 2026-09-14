#!/usr/bin/env python3

"""
Build many qlat nix packages at once.\n
This replaces the former seven 'build-many-qlat-pkgs-*.sh' scripts.  Each of them
is a '--group' here, and the arguments they pass to 'many-qlat-pkgs.nix' can be
overridden on the command line, so custom package sets are possible too.\n
The result symlink is '$HOME/qlat-build/nix/<group>/result', as before.  When
'$HOME' is not writable it falls back to
'<repo>/tmp/qlat-build/nix/<group>/result', and the nix and nom cache
directories fall back to '<repo>/tmp/' in the same way, so nothing outside of
the repo has to be written to.
"""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# '<repo>/nixpkgs', with symlinks resolved like 'pwd -P' does.
SCRIPT_PATH = Path(__file__).resolve().parent
REPO_PATH = SCRIPT_PATH.parent
NIXPKGS_FILE = SCRIPT_PATH / "many-qlat-pkgs.nix"

DEFAULT_OUT_DIR = Path.home() / "qlat-build" / "nix"
REPO_OUT_DIR = REPO_PATH / "tmp" / "qlat-build" / "nix"

# A nix (or JSON) double quoted string, with the escapes left in place.
NIX_STR_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')

# The values many-qlat-pkgs.nix falls back to when the matching argument is not
# passed, spelled out so that GROUPS states the actual arguments of every group.
# 'ALL_VERSIONS' is its built in version list and 'ALL_NAMES' is every name of
# 'q-pkgs.qlat-name-list' (defined through options.nix).  Keep 'ALL_NAMES' in
# sync when a variant is added there, otherwise the groups that build every name
# would silently miss it.
ALL_VERSIONS = ["", "26.05", "25.11"]
ALL_NAMES = [
    "",
    "-std",
    "-cpsless",
    "-gridless",
    "-cu",
    "-cuda",
    "-cudasupport",
    "-ucxless",
    "-clang",
    "-pypi",
    "-clang-ucxless",
    "-cuda-ucxless",
    "-gridless-cubaquadless",
    "-gridless-clang",
    "-std-ucxless",
    "-std-clang-ucxless",
    "-cpsless-ucxless",
    "-cpsless-clang-ucxless",
    "-cpsless-clang",
    "-std-cu",
    "-std-cuda",
    "-std-cudasupport",
]

# The package sets the old scripts built.  Every group spells out every field,
# so reading this table is enough to know what a group builds.  'qlat_tests' is
# 'all' (no test filter), 'none' (no test packages), 'none-cuda' (every test
# package but the '-cuda*' ones), 'only' (only test packages) or 'only-cuda'
# (only the '-cuda*' test packages), and 'jobs'/'cores' are the '-j'/'--cores'
# defaults.
GROUPS = {
    "all": {
        "description": "all qlat packages (skip cuda tests)",
        "qlat_tests": "none-cuda",
        "version_list": ALL_VERSIONS,
        "qlat_name_list": ALL_NAMES,
        "jobs": 8,
        "cores": 7,
    },
    "all-cuda-tests": {
        "description": "all qlat packages (cuda tests only)",
        "qlat_tests": "only-cuda",
        "version_list": ALL_VERSIONS,
        "qlat_name_list": ALL_NAMES,
        "jobs": 2,
        "cores": 15,
    },
    "core": {
        "description": "core packages only",
        "qlat_tests": "all",
        "version_list": [""],
        "qlat_name_list": ["", "-pypi"],
        "jobs": 4,
        "cores": 15,
    },
    "tests": {
        "description": "all tests (default version and name)",
        "qlat_tests": "only",
        "version_list": [""],
        "qlat_name_list": [""],
        "jobs": 4,
        "cores": 31,
    },
    "cuda-core": {
        "description": "cuda core packages",
        "qlat_tests": "all",
        "version_list": [""],
        "qlat_name_list": ["", "-clang", "-cudasupport", "-pypi"],
        "jobs": 6,
        "cores": 15,
    },
    "cuda": {
        "description": "cuda packages (skip cuda tests)",
        "qlat_tests": "none-cuda",
        "version_list": ALL_VERSIONS,
        "qlat_name_list": ["", "-clang", "-cu", "-cudasupport", "-ucxless", "-pypi"],
        "jobs": 4,
        "cores": 15,
    },
    "cuda-tests": {
        "description": "cuda packages (cuda tests only)",
        "qlat_tests": "only-cuda",
        "version_list": ALL_VERSIONS,
        "qlat_name_list": ["", "-clang", "-cu", "-cudasupport", "-ucxless", "-pypi"],
        "jobs": 2,
        "cores": 15,
    },
    "small": {
        "description": "small package set",
        "qlat_tests": "all",
        "version_list": ALL_VERSIONS,
        "qlat_name_list": ["", "-clang", "-ucxless", "-pypi"],
        "jobs": 4,
        "cores": 15,
    },
}

# The group that is built when '--group' is not given.
DEFAULT_GROUP = "core"

EPILOG = f"""
groups ({", ".join(GROUPS)}), '{DEFAULT_GROUP}' is built when '--group' is omitted:
  --group all             {GROUPS["all"]["description"]}
  --group all-cuda-tests  {GROUPS["all-cuda-tests"]["description"]}
  --group core            {GROUPS["core"]["description"]}
  --group tests           {GROUPS["tests"]["description"]}
  --group cuda-core       {GROUPS["cuda-core"]["description"]}
  --group cuda            {GROUPS["cuda"]["description"]}
  --group cuda-tests      {GROUPS["cuda-tests"]["description"]}
  --group small           {GROUPS["small"]["description"]}

'--group' may be repeated to build several sets one after another, and
'--version-list', '--qlat-name-list' and '--qlat-tests' override what the group
passes to many-qlat-pkgs.nix:
  --group small --qlat-name-list '["" "-cuda"]' --qlat-tests only
'--qlat-tests all' omits the argument, which is the many-qlat-pkgs.nix default
of building the env and test packages of every selected name.  The other values
select test packages: 'none' (none), 'none-cuda' (all but the '-cuda*' ones),
'only' (all, and no env packages) and 'only-cuda' (just the '-cuda*' ones, and
no env packages).  Other arguments of many-qlat-pkgs.nix ('ngpu',
'cudaCapability', ...) can be passed on to nix-build directly, since arguments
that are not recognised here are forwarded:
  --group small --argstr ngpu 2

outputs:
  $HOME/qlat-build/nix/<group>/result      the nix-build out link; when '$HOME'
                                           is not writable this is
                                           '<repo>/tmp/qlat-build/nix/<group>/result'
  <repo>/tmp/nix-cache, <repo>/tmp/nix-state
                                           nix and nom caches, used when the
                                           user level directories are not writable

'--out-dir' changes the base of the out link and '--out-link' sets it exactly;
both are used as given, without the '$HOME' fallback.
"""

def die(message, code=1):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(code)

def check_groups():
    """Fail loudly when a GROUPS entry has a malformed field.\n
    A wrong value here would silently change what is built, so the table is
    validated before anything is passed to nix-build.
    """
    if DEFAULT_GROUP not in GROUPS:
        die(f"DEFAULT_GROUP '{DEFAULT_GROUP}' is not one of {', '.join(GROUPS)}")
    for group, preset in GROUPS.items():
        for name in ("version_list", "qlat_name_list"):
            value = preset[name]
            if not isinstance(value, list) or not all(
                isinstance(item, str) for item in value
            ):
                die(
                    f"GROUPS['{group}']['{name}'] must be a list of strings, "
                    f"got: {value!r}"
                )
        value = preset["qlat_tests"]
        if value not in ("all", "none", "none-cuda", "only", "only-cuda"):
            die(
                f"GROUPS['{group}']['qlat_tests'] must be 'all', 'none', "
                f"'none-cuda', 'only' or 'only-cuda', got: {value!r}"
            )
        for name in ("jobs", "cores"):
            if not isinstance(preset[name], int):
                die(
                    f"GROUPS['{group}']['{name}'] must be an int, got: {preset[name]!r}"
                )

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

def nix_list(values):
    """Serialize 'values' as a nix list literal, e.g. '["" "-pypi"]'.\n
    nix lists are whitespace separated; a comma is a syntax error there, so
    'json.dumps' output must not be used.
    """
    return "[" + " ".join(nix_str(value) for value in values) + "]"

def unescape_nix_str(value):
    result = []
    index = 0
    while index < len(value):
        char = value[index]
        if char == "\\" and index + 1 < len(value):
            index += 1
            nxt = value[index]
            result.append({"n": "\n", "r": "\r", "t": "\t"}.get(nxt, nxt))
        else:
            result.append(char)
        index += 1
    return "".join(result)

def parse_nix_list(value):
    """Parse a list of quoted strings, e.g. '["" "-pypi"]'.\n
    Commas are accepted as separators too, so the JSON spelling
    '["", "-pypi"]' works as well.
    """
    matches = NIX_STR_RE.findall(value)
    remainder = NIX_STR_RE.sub("", value).strip("[], \t\n")
    if not matches or remainder:
        raise argparse.ArgumentTypeError(
            f'expected a list of quoted strings, e.g. \'["" "-pypi"]\', got: {value}'
        )
    return [unescape_nix_str(match) for match in matches]

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Build many qlat nix packages at once.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EPILOG,
    )
    parser.add_argument(
        "--group",
        action="append",
        choices=sorted(GROUPS),
        metavar="NAME",
        help=f"package set to build, may be repeated (default: {DEFAULT_GROUP}; "
        "see the list below)",
    )
    parser.add_argument(
        "--list-groups",
        action="store_true",
        help="list the package sets and exit",
    )
    parser.add_argument(
        "--version-list",
        type=parse_nix_list,
        metavar="LIST",
        help="override the group's '--arg version-list', e.g. '[\"\" \"26.05\"]'",
    )
    parser.add_argument(
        "--qlat-name-list",
        type=parse_nix_list,
        metavar="LIST",
        help="override the group's '--arg qlat-name-list', e.g. '[\"\" \"-cuda\"]'",
    )
    parser.add_argument(
        "--qlat-tests",
        choices=["all", "none", "none-cuda", "only", "only-cuda"],
        help="override the group's '--argstr qlat-tests'; 'all' omits the "
        "argument, which is the many-qlat-pkgs.nix default of building the env "
        "and test packages of every selected name",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        metavar="DIR",
        help=f"base of the out link, i.e. DIR/<group>/result (default: {DEFAULT_OUT_DIR})",
    )
    parser.add_argument(
        "--out-link",
        type=Path,
        metavar="PATH",
        help="exact path of the result symlink (default: <out dir>/<group>/result)",
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
        help="nix-build -j (default: the group's value, see --list-groups)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        help="nix-build --cores (default: the group's value, see --list-groups)",
    )
    args, nix_build_args = parser.parse_known_args(argv)
    args.nix_build_args = nix_build_args
    #
    if args.list_groups:
        return args
    if not args.group:
        args.group = [DEFAULT_GROUP]
    if args.out_link is not None and args.out_dir is not None:
        parser.error("--out-link and --out-dir cannot be combined")
    if args.out_link is not None and len(args.group) > 1:
        parser.error("--out-link cannot be combined with several --group values")
    return args

def describe_group(group):
    """Every nix argument of a group, as the values that are actually used."""
    preset = GROUPS[group]
    return [
        ("version-list", nix_list(preset["version_list"])),
        ("qlat-name-list", nix_list(preset["qlat_name_list"])),
        ("qlat-tests", preset["qlat_tests"]),
        ("jobs/cores", f"-j {preset['jobs']} --cores {preset['cores']}"),
    ]

def print_groups():
    print("groups (use --group NAME, which may be repeated):")
    print("'qlat-tests: all' means no test filter, in which case that argument is")
    print("not passed (many-qlat-pkgs.nix accepts null, 'none', 'none-cuda', 'only'")
    print("or 'only-cuda').  'none'/'none-cuda' leave out test packages,")
    print("'only'/'only-cuda' leave out the env packages and keep only tests.")
    for group in GROUPS:
        print()
        print(f"  {group}" + (" (default)" if group == DEFAULT_GROUP else ""))
        print(f"      {GROUPS[group]['description']}")
        for name, value in describe_group(group):
            print(f"      {name + ':':<16} {value}")
        print(f"      {'out link:':<16} {DEFAULT_OUT_DIR / group / 'result'}")

def resolve_out_link(args, group):
    if args.out_link is not None:
        out_link = args.out_link
    else:
        if args.out_dir is not None:
            out_dir = args.out_dir
        else:
            out_dir = DEFAULT_OUT_DIR
            if not can_write_dir(out_dir):
                print(
                    f"The user build directory '{out_dir}' is not writable. "
                    f"Using '{REPO_OUT_DIR}'."
                )
                out_dir = REPO_OUT_DIR
        out_link = out_dir / group / "result"
    #
    try:
        out_link.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        die(f"cannot create '{out_link.parent}': {e}")
    return out_link

def remove_stale_out_link(out_link):
    """Remove an out link that nix-build would refuse to replace.\n
    nix-build only refuses to replace a symlink that does not point into the nix
    store, so a valid result link from an earlier build is left alone.
    """
    if not out_link.is_symlink():
        return
    target = os.path.realpath(out_link)
    if target.startswith("/nix/store/") and os.path.exists(target):
        return
    print(f"Removing stale out link '{out_link}'.")
    out_link.unlink()

def group_nix_args(args, group):
    """The '--arg'/'--argstr' arguments for a group, with the overrides applied.\n
    'None' on the parsed arguments means the option was not given on the command
    line, in which case the group's own value is used.  A 'qlat_tests' of 'all'
    means no test filter, and many-qlat-pkgs.nix accepts only null, 'none',
    'none-cuda', 'only' or 'only-cuda', so the argument is left out for it.
    """
    preset = GROUPS[group]
    #
    version_list = args.version_list
    if version_list is None:
        version_list = preset["version_list"]
    #
    qlat_name_list = args.qlat_name_list
    if qlat_name_list is None:
        qlat_name_list = preset["qlat_name_list"]
    #
    qlat_tests = args.qlat_tests
    if qlat_tests is None:
        qlat_tests = preset["qlat_tests"]
    #
    result = [
        "--arg",
        "version-list",
        nix_list(version_list),
        "--arg",
        "qlat-name-list",
        nix_list(qlat_name_list),
    ]
    if qlat_tests != "all":
        result += ["--argstr", "qlat-tests", qlat_tests]
    return result

def build_group(args, group, env):
    preset = GROUPS[group]
    out_link = resolve_out_link(args, group)
    jobs = args.jobs if args.jobs is not None else preset["jobs"]
    cores = args.cores if args.cores is not None else preset["cores"]
    #
    print()
    print(f"Building group '{group}': {preset['description']}.")
    print(f"Out link: '{out_link}'.  -j {jobs} --cores {cores}.")
    remove_stale_out_link(out_link)
    #
    cmd = ["nix-build", str(NIXPKGS_FILE), "-o", str(out_link)]
    cmd += group_nix_args(args, group)
    if not args.no_nom:
        cmd += ["--log-format", "internal-json", "-v"]
    cmd += ["-j", str(jobs), "--cores", str(cores)]
    cmd += args.nix_build_args
    #
    if args.no_nom:
        run(cmd, env)
    else:
        run_with_nom(cmd, env)

def main(argv=None):
    sys.stdout.reconfigure(line_buffering=True)
    start_time = time.monotonic()
    args = parse_args(argv)
    #
    check_groups()
    #
    if args.list_groups:
        print_groups()
        return 0
    #
    if not NIXPKGS_FILE.is_file():
        die(f"'{NIXPKGS_FILE}' not found.  This script belongs in 'nixpkgs/'.")
    #
    env = dict(os.environ)
    setup_cache_fallbacks(env)
    for group in args.group:
        build_group(args, group, env)
    #
    print()
    print(f"Finished successfully in {time.monotonic() - start_time:.1f} s.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
