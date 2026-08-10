#!/usr/bin/env python3
"""
Check and upgrade pinned external package versions in nixpkgs/*.nix.

How to upgrade:
  1. Check which packages have updates available:
       ./nixpkgs/upgrade.py
     Each row shows: name, .nix file, current version, latest version, status.
     Status "UPDATE" means an upgrade is available.
  2. Apply the upgrades:
       ./nixpkgs/upgrade.py --update              # upgrade everything outdated
       ./nixpkgs/upgrade.py --update gvar lsqfit  # upgrade only selected packages
     Only the version pins in the .nix files are edited; nothing is built.
  3. Review the changes:
       git diff nixpkgs/
  4. Rebuild qlat with the new versions (this refetches the upgraded sources):
       name='' ./nixpkgs/install-py-local-kernel-with-nix.sh
  5. Run the test suite to make sure the new versions work:
       nix-build nixpkgs/q-pkgs.nix -A pkgs.qlat-tests -j 4 --cores 31
  6. Commit the updated .nix files.

Usage:
  ./nixpkgs/upgrade.py                 # check all packages, print status table
  ./nixpkgs/upgrade.py --update        # apply available upgrades to the .nix files
  ./nixpkgs/upgrade.py gpt grid gvar   # check (or with --update, upgrade) selected packages
                                       # (unique name prefixes are accepted)

Requirements:
  - Network access to github.com (`git ls-remote` and the GitHub API are used
    to find the latest commits, tags, and commit dates).
  - `nix` on PATH, only needed for tag-hash packages (qcd_ml) to recompute
    the src hash.

Package kinds:
  rev       pin a git commit of a branch (builtins.fetchGit rev). Upgrading moves
            the pin to the branch HEAD; the old pin is kept as a comment with its
            date, following the existing history style in the .nix files.
  tag       pin a release tag. Upgrading bumps `version = "..."` (and the literal
            tag in `ref = ...` for usqcd-style tags such as c-lime1-3-2).
  tag-hash  like tag, but the src hash must be recomputed (requires `nix`).
  manual    tarball hosted in the Qlattice-distfiles repo; reported only.
            To upgrade these, add the new tarball to the Qlattice-distfiles repo
            and edit version + hash in the .nix file by hand.
"""

import argparse
import datetime
import json
import os
import re
import subprocess
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))

PKGS = [
    dict(name="grid-lehner", file="grid-lehner.nix", kind="rev",
         repo="jinluchang/Grid-clehner", branch="feature/gpt", field="version"),
    dict(name="gpt-lehner", file="gpt-lehner.nix", kind="rev",
         repo="jinluchang/gpt", branch="master", field="version"),
    dict(name="cps", file="cps.nix", kind="rev",
         repo="RBC-UKQCD/CPS_public", branch="master", field="version"),
    dict(name="gvar", file="gvar.nix", kind="tag",
         repo="gplepage/gvar", tag_prefix="v"),
    dict(name="lsqfit", file="lsqfit.nix", kind="tag",
         repo="gplepage/lsqfit", tag_prefix="v"),
    dict(name="corrfitter", file="corrfitter.nix", kind="tag",
         repo="gplepage/corrfitter", tag_prefix="v"),
    dict(name="vegas", file="vegas.nix", kind="tag",
         repo="gplepage/vegas", tag_prefix="v"),
    dict(name="c-lime", file="c-lime.nix", kind="tag",
         repo="usqcd-software/c-lime", tag_prefix="c-lime", tag_sep="-"),
    dict(name="qmp", file="qmp.nix", kind="tag",
         repo="usqcd-software/qmp", tag_prefix="qmp", tag_sep="-"),
    dict(name="qcd_ml", file="qcd_ml.nix", kind="tag-hash",
         repo="daknuett/qcd_ml", tag_prefix="v"),
    dict(name="nixgl", file="options.nix", kind="rev",
         repo="jinluchang/nixGL", branch="main", field="rev"),
    dict(name="qio", file="qio.nix", kind="manual"),
    dict(name="cuba-quad", file="cuba-quad.nix", kind="manual"),
]

### -------------------------------------------------------------------


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def git_url(repo):
    return f"https://github.com/{repo}"


def pkg_path(pkg):
    return os.path.join(HERE, pkg["file"])


def read_pkg(pkg):
    with open(pkg_path(pkg)) as f:
        return f.read()


def current_version(pkg):
    """First non-commented `field = "...";` line in the .nix file."""
    field = pkg.get("field", "version")
    m = re.search(rf'^(\s*){field} = "([^"]+)";', read_pkg(pkg), re.M)
    assert m, f"{pkg['file']}: cannot find current {field} = ..."
    return m.group(2)


def ls_remote_head(repo, branch):
    out = run(["git", "ls-remote", git_url(repo), f"refs/heads/{branch}"])
    return out.split()[0]


def ls_remote_tags(repo):
    out = run(["git", "ls-remote", "--tags", git_url(repo)])
    tags = set()
    for line in out.splitlines():
        ref = line.split("\t")[1]
        tag = ref[len("refs/tags/"):]
        if tag.endswith("^{}"):
            continue
        tags.add(tag)
    return sorted(tags)


def parse_version(s):
    return tuple(int(x) for x in re.findall(r"\d+", s))


def latest_tag(pkg):
    """Return (version_str, tag) of the highest release tag, or None."""
    prefix = pkg["tag_prefix"]
    sep = pkg.get("tag_sep", ".")
    cands = []
    for tag in ls_remote_tags(pkg["repo"]):
        if not tag.startswith(prefix):
            continue
        ver = tag[len(prefix):]
        if sep == "-":
            ver = ver.replace("-", ".")
        if not re.fullmatch(r"\d+(\.\d+)*", ver):
            continue
        cands.append((parse_version(ver), ver, tag))
    if not cands:
        return None
    cands.sort()
    _, ver, tag = cands[-1]
    return ver, tag


def commit_date(repo, sha):
    url = f"https://api.github.com/repos/{repo}/commits/{sha}"
    try:
        with urllib.request.urlopen(url, timeout=20) as f:
            data = json.load(f)
        return data["commit"]["committer"]["date"][:10].replace("-", "/")
    except Exception:
        return datetime.date.today().strftime("%Y/%m/%d")


def prefetch_unpack_hash(url):
    """SRI sha256 of the unpacked tarball, as used by fetchzip/fetchFromGitHub."""
    try:
        out = run(["nix", "--extra-experimental-features", "nix-command",
                   "store", "prefetch-file", "--json", "--unpack", url])
        return json.loads(out)["hash"]
    except Exception:
        pass
    try:
        base32 = run(["nix-prefetch-url", "--unpack", url]).strip()
        out = run(["nix", "--extra-experimental-features", "nix-command",
                   "hash", "convert", "--hash-algo", "sha256", "--to", "sri", base32])
        return out.strip()
    except Exception:
        raise RuntimeError(f"cannot prefetch hash for {url}; "
                           f"install nix or run `nix store prefetch-file --unpack {url}` manually")


### -------------------------------------------------------------------


def edit_file(pkg, transform):
    path = pkg_path(pkg)
    with open(path) as f:
        text = f.read()
    new_text = transform(text)
    assert new_text != text, f"{pkg['file']}: edit made no change"
    with open(path, "w") as f:
        f.write(new_text)


def update_rev(pkg, new_sha, date):
    field = pkg.get("field", "version")
    old = current_version(pkg)

    def transform(text):
        m = re.search(rf'^(\s*){field} = "{re.escape(old)}";([^\n]*)$', text, re.M)
        assert m, f"{pkg['file']}: cannot locate active {field} line"
        indent, comment = m.group(1), m.group(2).rstrip()
        old_line = f"{indent}# {field} = \"{old}\";{comment}"
        new_line = f"{indent}{field} = \"{new_sha}\"; # {date}"
        return text[:m.start()] + old_line + "\n" + new_line + text[m.end():]

    edit_file(pkg, transform)


def update_tag(pkg, new_ver, new_tag):
    old = current_version(pkg)

    def transform(text):
        text, n = re.subn(rf'^(\s*)version = "{re.escape(old)}";',
                          rf'\g<1>version = "{new_ver}";', text, count=1, flags=re.M)
        assert n == 1, f"{pkg['file']}: cannot locate version line"
        if pkg.get("tag_sep") == "-":
            old_tag = pkg["tag_prefix"] + old.replace(".", "-")
            text, n = re.subn(rf'ref = "refs/tags/{re.escape(old_tag)}";',
                              f'ref = "refs/tags/{new_tag}";', text, count=1)
            assert n == 1, f"{pkg['file']}: cannot locate ref tag line"
        return text

    edit_file(pkg, transform)


def update_tag_hash(pkg, new_ver, new_tag):
    update_tag(pkg, new_ver, new_tag)
    url = f"https://github.com/{pkg['repo']}/archive/refs/tags/{new_tag}.tar.gz"
    new_hash = prefetch_unpack_hash(url)

    def transform(text):
        text, n = re.subn(r'hash = "sha256-[^"]+";', f'hash = "{new_hash}";',
                          text, count=1)
        assert n == 1, f"{pkg['file']}: cannot locate hash line"
        return text

    edit_file(pkg, transform)


### -------------------------------------------------------------------


def check_pkg(pkg):
    """Return dict with current, latest, status, and data needed for update."""
    res = dict(name=pkg["name"], file=pkg["file"], kind=pkg["kind"],
               current=None, latest=None, status=None, data=None)
    try:
        res["current"] = current_version(pkg)
        if pkg["kind"] == "manual":
            res["latest"] = "-"
            res["status"] = "manual (distfiles tarball)"
        elif pkg["kind"] == "rev":
            sha = ls_remote_head(pkg["repo"], pkg["branch"])
            res["latest"] = sha
            res["data"] = dict(sha=sha)
            res["status"] = "up-to-date" if sha == res["current"] else "UPDATE"
        else:
            r = latest_tag(pkg)
            if r is None:
                res["latest"] = "?"
                res["status"] = "no tags found"
            else:
                ver, tag = r
                res["latest"] = ver
                res["data"] = dict(ver=ver, tag=tag)
                res["status"] = ("up-to-date" if parse_version(ver) == parse_version(res["current"])
                                 else "UPDATE")
    except Exception as e:
        res["status"] = f"error: {e}"
    return res


def apply_update(pkg, res):
    if pkg["kind"] == "rev":
        sha = res["data"]["sha"]
        date = commit_date(pkg["repo"], sha)
        update_rev(pkg, sha, date)
    elif pkg["kind"] == "tag":
        update_tag(pkg, res["data"]["ver"], res["data"]["tag"])
    elif pkg["kind"] == "tag-hash":
        update_tag_hash(pkg, res["data"]["ver"], res["data"]["tag"])
    else:
        raise RuntimeError(f"cannot auto-update kind {pkg['kind']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pkgs", nargs="*",
                        help="package names to check (default: all); "
                             "unique prefixes such as 'gpt' or 'grid' are accepted")
    parser.add_argument("--update", action="store_true",
                        help="apply available upgrades to the .nix files")
    args = parser.parse_args()

    selected = []
    for name in args.pkgs:
        matches = [p for p in PKGS if p["name"] == name or p["name"].startswith(name)]
        if not matches:
            parser.error(f"unknown package: {name}")
        for m in matches:
            if m not in selected:
                selected.append(m)
    if not selected:
        selected = PKGS

    results = [check_pkg(p) for p in selected]

    width = max(len(r["name"]) for r in results)
    for r in results:
        print(f"{r['name']:<{width}}  {r['file']:<18}  {str(r['current']):<44}  "
              f"{str(r['latest']):<44}  {r['status']}")

    if not args.update:
        print()
        print(__doc__.strip())
        return

    updated = []
    for pkg, res in zip(selected, results):
        if res["status"] != "UPDATE":
            continue
        print(f"updating {pkg['name']} ({pkg['file']}): {res['current']} -> {res['latest']}")
        apply_update(pkg, res)
        updated.append(pkg["name"])

    if updated:
        print()
        print(f"updated: {', '.join(updated)}")
        print("rebuild with: name='' ./nixpkgs/install-py-local-kernel-with-nix.sh")
    else:
        print("nothing to update")


if __name__ == "__main__":
    main()
