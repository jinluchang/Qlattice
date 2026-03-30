#!/usr/bin/env python3
"""
Qlattice Release Automation Script.

A single-file Python script that automates the entire Qlattice release procedure.
Uses only Python standard library (no external dependencies).

The release process consists of 9 sequential steps:
    1. Environment checks (VERSION file exists, git root, clean worktree)
    2. Build-pull loop (build tarballs, verify artifacts, pull remote changes, rebuild if needed)
    3. Tag the current commit with the version from VERSION file
    4. Push commit and tags to both origin (GitHub) and gitee (mirror)
    5. Upload built packages to PyPI
    6. Bump version in VERSION file, update-sources.sh, and nix config
    7. Verify new build compiles successfully
    8. Commit the version bump changes
    9. Push the version bump to both remotes

All steps use hard-fail error handling: any failure causes the script to exit
with a non-zero status. The script is fully automated with no interactive prompts.

Usage:
    python3 release.py

Requirements:
    - Script must be run from the project root (same directory as VERSION file)
    - Git working tree must be clean (no uncommitted changes)
    - ./VERSION file must exist and contain version string (e.g., "v0.90")
    - Remote "origin" (GitHub) and "gitee" (mirror) must be configured
    - Build dependencies: nix, twine, git
    - PyPI credentials configured via ~/.pypirc or environment variables

Design goals:
    - Prevent releasing from wrong directory or with uncommitted changes
    - Handle race conditions where remote commits arrive during 15-40 min build
    - Mark release commit with git tag for traceability
    - Retry unreliable gitee pushes up to 3 times with 30s delays
    - Safe re-runs via twine --skip-existing
"""

from pathlib import Path
from typing import Optional
import subprocess
import glob
import sys
import re
import time


def print_step(n: int, description: str) -> None:
    """Print step number and description for user feedback."""
    print(f"Step {n}: {description}")


def run(cmd, check: bool = True, capture_output: bool = True, text: bool = True, cwd: Optional[Path] = None):
    """Run a subprocess command with configurable options.

    Args:
        cmd: Command to run (string or list).
        check: If True, raise CalledProcessError on non-zero exit. Default True.
        capture_output: If True, capture stdout/stderr. Default True.
        text: If True, decode output as text. Default True.
        cwd: Optional working directory for the command.

    Returns:
        subprocess.CompletedProcess result.
    """
    result = subprocess.run(
        cmd,
        check=False if not check else True,
        stdout=subprocess.PIPE if capture_output else None,
        stderr=subprocess.STDOUT,
        text=text,
        cwd=str(cwd) if cwd else None,
    )
    return result


def check_environment() -> None:
    """Step 1: Verify the environment is suitable for release.

    Checks:
        - VERSION file exists at project root
        - Script is running from git repository root
        - Working tree is clean (no uncommitted changes)

    Exits with error message if any check fails.
    """
    print_step(1, "Environment sanity checks: VERSION existence, git root, and clean worktree")
    version_path = Path("VERSION")
    if not version_path.exists():
        print("ERROR: VERSION file does not exist at project root.")
        sys.exit(1)

    try:
        toplevel = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=True,
        ).stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"ERROR: git rev-parse failed: {e}")
        sys.exit(1)

    script_parent = str(Path(__file__).resolve().parent)
    if toplevel != script_parent:
        print(
            f"ERROR: git top-level '{toplevel}' does not match script dir '{script_parent}'"
        )
        sys.exit(1)

    status = subprocess.run(
        ["git", "status", "--porcelain"], stdout=subprocess.PIPE, text=True
    ).stdout.strip()
    if status:
        print("ERROR: Working tree is not clean. Commit or stash changes before proceeding.")
        sys.exit(1)

    print("Environment checks passed.")


def build_pull_loop() -> None:
    """Step 2: Build tarballs with pull-to-rebuild loop.

    Runs the build script, verifies artifacts exist, then pulls from remote.
    If new commits arrived during the build, rebuilds to include them.
    Repeats until no new commits are pulled (converges to stable state).

    This handles the race condition where someone pushes to the repo during
    the 15-40 minute build — ensures the release always contains the latest code.

    Exits on build failure, missing artifacts, or merge conflicts.
    """
    print_step(2, "Build loop: build tarballs, verify, pull and re-build on updates")
    while True:
        # Step A: build
        try:
            subprocess.run(["./nixpkgs/build-many-qlat-pkgs-core.sh"], check=True)
        except subprocess.CalledProcessError:
            print("ERROR: Build script failed.")
            sys.exit(1)

        # Step B: verify artifacts exist
        tar_pattern = str(Path.home() / "qlat-build" / "nix" / "result-*-q-pkgs" /
                          "share" / "qlat-pypi" / "qlat*.tar.gz")
        tarballs = glob.glob(tar_pattern)
        if not tarballs:
            print(f"ERROR: Build artifacts not found matching {tar_pattern}")
            sys.exit(1)

        # Step C: record HEAD and pull
        head_before = subprocess.run(["git", "rev-parse", "HEAD"], stdout=subprocess.PIPE, text=True).stdout.strip()
        pull = subprocess.run(["git", "pull", "origin"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if pull.returncode != 0:
            out = pull.stdout or ""
            if "CONFLICT" in out or "conflict" in out:
                subprocess.run(["git", "merge", "--abort"], check=True)
                print("ERROR: Merge conflict during git pull. Aborting release.")
                sys.exit(1)
            else:
                print("ERROR: git pull failed.")
                sys.exit(1)

        head_after = subprocess.run(["git", "rev-parse", "HEAD"], stdout=subprocess.PIPE, text=True).stdout.strip()
        if head_after != head_before:
            # New commits were pulled; re-run build
            continue
        else:
            # No updates from remote; finish loop
            break
    print("Build loop completed.")


def tag_version() -> None:
    """Step 3: Create a git tag for the current version.

    Reads version from VERSION file and creates a git tag.
    Fails fast if the tag already exists to prevent accidental overwrites.
    """
    print_step(3, "Tag version from VERSION file")
    version = Path("VERSION").read_text().strip()
    if not version:
        print("ERROR: VERSION file is empty.")
        sys.exit(1)
    # Ensure the tag doesn't already exist
    existing = subprocess.run(["git", "tag", "--list", version], stdout=subprocess.PIPE, text=True).stdout.strip()
    if existing:
        print(f"ERROR: tag {version} already exists — delete it or bump VERSION")
        sys.exit(1)
    subprocess.run(["git", "tag", version], check=True)
    print(f"Tagged {version}.")


def push_contents_and_tags() -> None:
    """Step 4: Push commit and tags to both remotes.

    Pushes branch contents and git tags to origin (GitHub) and gitee (mirror).
    Includes retry logic for gitee which can be unreliable:
        - Up to 3 attempts per push operation
        - 30 second sleep between attempts
    """
    print_step(4, "Push contents and tags to remotes (origin, gitee)")
    subprocess.run(["git", "push", "origin"], check=True)
    subprocess.run(["git", "push", "origin", "--tags"], check=True)

    # Push to gitee with retries
    for remote in ["gitee"]:
        ok = False
        for attempt in range(3):
            res = subprocess.run(["git", "push", remote], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if res.returncode == 0:
                ok = True
                break
            time.sleep(30)
        if not ok:
            print(f"ERROR: git push {remote} failed after retries.")
            sys.exit(1)

    # Push tags to gitee as well
    for remote in ["gitee"]:
        ok = False
        for attempt in range(3):
            res = subprocess.run(["git", "push", remote, "--tags"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if res.returncode == 0:
                ok = True
                break
            time.sleep(30)
        if not ok:
            print(f"ERROR: git push {remote} --tags failed after retries.")
            sys.exit(1)

    print("Pushes to origins completed.")


def upload_pypi() -> None:
    """Step 5: Upload built packages to PyPI.

    Uploads all tarballs matching qlat*.tar.gz from the build output directory.
    Uses twine with --skip-existing to safely handle re-runs without failing
    on already-uploaded versions.
    """
    print_step(5, "Upload tarballs to PyPI (twine)")
    tar_pattern = str(Path.home() / "qlat-build" / "nix" / "result-*-q-pkgs" /
                          "share" / "qlat-pypi" / "qlat*.tar.gz")
    tarballs = glob.glob(tar_pattern)
    if not tarballs:
        print(f"ERROR: No tarballs found for PyPI upload matching {tar_pattern}")
        sys.exit(1)
    cmd = ["twine", "upload", "--skip-existing"] + tarballs
    subprocess.run(cmd, check=True)
    print("PyPI upload completed.")


def bump_version() -> None:
    """Step 6: Bump version for next development cycle.

    Increments the minor version number in three locations:
        - ./VERSION file (format: vX.Y)
        - Source files via ./scripts/update-sources.sh
        - version-pypi in nixpkgs/q-pkgs.nix (format: X.Y, no 'v' prefix)

    Exits with error if VERSION format is not recognized.
    """
    print_step(6, "Bump minor version, update sources and nix q-pkgs")
    version = Path("VERSION").read_text().strip()
    m = re.match(r"^v?(\d+)\.(\d+.*)$", version)
    if not m:
        print("ERROR: VERSION format not recognized (expected vX.Y).")
        sys.exit(1)
    major = int(m.group(1))
    minor_str = m.group(2)
    # Increment the numeric minor component (treat non-numeric tail as part of minor)
    try:
        minor_num = int(minor_str.split(".")[0]) if "." in minor_str else int(minor_str)
    except ValueError:
        minor_num = 0
    minor_num += 1
    new_version = f"v{major}.{minor_num}"
    Path("VERSION").write_text(new_version + "\n")

    # Update sources
    subprocess.run(["./scripts/update-sources.sh"], check=True)

    # Update version-pypi in nix file (bare number, no 'v')
    q_pkgs_nix = Path("nixpkgs/q-pkgs.nix")
    if not q_pkgs_nix.exists():
        print("ERROR: nixpkgs/q-pkgs.nix not found.")
        sys.exit(1)
    text = q_pkgs_nix.read_text()
    new_pypi_version = f"{major}.{minor_num}"
    # Use \g<1> to avoid ambiguity when replacement starts with digits
    text_new = re.sub(r'(version-pypi\s*=\s*")[^"]+(")', r"\g<1>" + new_pypi_version + r"\g<2>", text)
    if text == text_new:
        print("WARNING: version-pypi pattern not found in q-pkgs.nix — manual update needed.")
    q_pkgs_nix.write_text(text_new)
    print(f"Version bumped to {new_version}, version-pypi updated to {new_pypi_version}.")


def verify_new_build() -> None:
    """Step 7: Verify the new version builds successfully.

    Re-runs the build script with the bumped version.
    Catches version-related build breakage immediately before committing.
    """
    print_step(7, "Verify new build by rebuilding packages")
    subprocess.run(["./nixpkgs/build-many-qlat-pkgs-core.sh"], check=True)
    print("New build verification completed.")


def commit_bump() -> None:
    """Step 8: Commit the version bump changes.

    Creates a single atomic commit recording all version bump changes
    (VERSION file, source files, nix config) with message "Upgrade version".
    """
    print_step(8, "Commit the version bump changes")
    subprocess.run(["git", "commit", "-am", "Upgrade version"], check=True)
    print("Commit created.")


def push_bump() -> None:
    """Step 9: Push the version bump to both remotes.

    Pushes the version bump commit to origin (GitHub) and gitee (mirror).
    Includes retry logic for gitee (up to 3 attempts, 30s sleep).
    """
    print_step(9, "Push the bump to remotes (origin and gitee)")
    subprocess.run(["git", "push", "origin"], check=True)
    # Push gitee with retries
    for attempt in range(3):
        res = subprocess.run(["git", "push", "gitee"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if res.returncode == 0:
            break
        time.sleep(30)
    else:
        print("ERROR: git push gitee failed after retries.")
        sys.exit(1)
    print("Push bump completed.")


def main() -> None:
    """Main entry point: runs all 9 release steps sequentially.

    Executes the complete release procedure:
        1. check_environment
        2. build_pull_loop
        3. tag_version
        4. push_contents_and_tags
        5. upload_pypi
        6. bump_version
        7. verify_new_build
        8. commit_bump
        9. push_bump

    Prints "SUCCESS" on completion, exits with non-zero on any failure.
    """
    check_environment()
    version = Path("VERSION").read_text().strip()
    print(f"Current version: {version}")
    build_pull_loop()
    tag_version()
    push_contents_and_tags()
    upload_pypi()
    bump_version()
    verify_new_build()
    commit_bump()
    push_bump()
    print("SUCCESS")


if __name__ == "__main__":
    main()
