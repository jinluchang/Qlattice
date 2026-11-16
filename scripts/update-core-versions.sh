#!/usr/bin/env bash

# Check and advance the pins in ./scripts/core-versions.sh.
#
# How to upgrade:
#   1. Check which repositories have a newer upstream branch tip:
#        ./scripts/update-core-versions.sh
#   2. Apply the upgrades:
#        ./scripts/update-core-versions.sh --update
#        ./scripts/update-core-versions.sh --update grid gpt   # selected only
#   3. Review the changes:
#        git diff scripts/core-versions.sh
#   4. Download/build with the new pins (this reuses the existing checkouts):
#        ./scripts/download-core.sh
#   5. Optionally advance the Nix pins of the same sources:
#        ./nixpkgs/upgrade.py --update gpt grid cps
#   6. Run the tests and commit scripts/core-versions.sh (and nixpkgs/*.nix).
#
# "Latest" means the tip of the branch recorded in scripts/core-versions.sh,
# resolved with `git ls-remote` (the same rule ./nixpkgs/upgrade.py uses). Pins
# are never moved automatically: without --update this is a dry run.
#
# Usage:
#   ./scripts/update-core-versions.sh [--update] [NAME...]
#
# Requirements: git and network access to github.com.

set -euo pipefail

wd="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
versions_file="$wd/scripts/core-versions.sh"

source "$versions_file"

update=false
selected=()
for arg in "$@" ; do
    case "$arg" in
        --update) update=true ;;
        -h|--help) sed -n '2,26p' "$0" ; exit 0 ;;
        -*) echo "Unknown option: $arg" >&2 ; exit 2 ;;
        *) selected+=("$arg") ;;
    esac
done

want() {
    # Case-insensitive prefix match (like ./nixpkgs/upgrade.py).
    local name lower
    name="$1"
    lower="$(echo "$name" | tr '[:upper:]' '[:lower:]')"
    if [ ${#selected[@]} -eq 0 ] ; then return 0 ; fi
    local s
    for s in "${selected[@]}" ; do
        s="$(echo "$s" | tr '[:upper:]' '[:lower:]')"
        if [ "${lower#"$s"}" != "$lower" ] ; then return 0 ; fi
    done
    return 1
}

resolve_tip() {
    local url="$1" branch="$2" out
    if ! out="$(timeout 120 git ls-remote "$url" "refs/heads/$branch" 2>/dev/null)" ; then
        return 1
    fi
    echo "$out" | cut -f1
}

# Committer date of a commit, like the dated pin comments in nixpkgs/*.nix.
# Falls back to today's date when the GitHub API is unavailable (rate limit).
commit_date() {
    local repo="$1" rev="$2"
    python3 - "$repo" "$rev" <<'PY' 2>/dev/null || date +%Y/%m/%d
import json, sys, urllib.request
repo, rev = sys.argv[1], sys.argv[2]
url = f"https://api.github.com/repos/{repo}/commits/{rev}"
with urllib.request.urlopen(url, timeout=20) as f:
    d = json.load(f)
print(d["commit"]["committer"]["date"][:10].replace("-", "/"))
PY
}

github_repo() {
    local repo="$1"
    repo="${repo#https://github.com/}"
    echo "${repo%.git}"
}

upd_rec=()
failed=0

for entry in "${CORE_REPOS[@]}" ; do
    IFS='|' read -r name url branch rev date <<< "$entry"
    if ! want "$name" ; then
        continue
    fi
    if latest="$(resolve_tip "$url" "$branch")" && [ -n "$latest" ] ; then
        if [ "$latest" = "$rev" ] ; then
            status="up-to-date"
        else
            status="UPDATE"
        fi
    else
        latest="?"
        status="error: cannot resolve refs/heads/$branch"
        failed=1
    fi
    printf '%-20s %-14s %s -> %s  %s\n' "$name" "$branch" "$rev" "$latest" "$status"
    if [ "$status" = "UPDATE" ] ; then
        upd_rec+=("$name|$url|$branch|$latest|$(commit_date "$(github_repo "$url")" "$latest")")
    fi
done

if [ "$failed" -ne 0 ] ; then
    echo
    echo "ERROR: some upstream branches could not be resolved; nothing changed." >&2
    exit 1
fi

if ! $update ; then
    echo
    echo "Dry run: nothing changed. Use --update to advance the pins."
    exit 0
fi

if [ ${#upd_rec[@]} -eq 0 ] ; then
    echo
    echo "All pins up-to-date; scripts/core-versions.sh unchanged."
    exit 0
fi

# Rewrite only the affected entry lines; comments and layout are preserved.
tmp="$(mktemp "$versions_file".tmp.XXXXXX)"
trap 'rm -f "$tmp"' EXIT
in_array=false
while IFS= read -r line || [ -n "$line" ] ; do
    if [ "$line" = "CORE_REPOS=(" ] ; then
        in_array=true
        echo "$line" >> "$tmp"
        continue
    fi
    if $in_array && [ "$line" = ")" ] ; then
        in_array=false
        echo "$line" >> "$tmp"
        continue
    fi
    if $in_array && [[ "$line" =~ ^([[:space:]]*\")([^|]+)\|([^|]+)\|([^|]+)\|([0-9a-f]{40})\|([^\"]*)(\")[[:space:]]*$ ]] ; then
        entry_name="${BASH_REMATCH[2]}"
        replaced=false
        for rec in "${upd_rec[@]}" ; do
            if [ "${rec%%|*}" = "$entry_name" ] ; then
                echo "${BASH_REMATCH[1]}$rec${BASH_REMATCH[7]}" >> "$tmp"
                replaced=true
                break
            fi
        done
        if $replaced ; then continue ; fi
    fi
    echo "$line" >> "$tmp"
done < "$versions_file"
chmod 644 "$tmp"
mv "$tmp" "$versions_file"
trap - EXIT

echo
echo "Updated ${#upd_rec[@]} pin(s) in scripts/core-versions.sh:"
git -C "$wd" diff -- scripts/core-versions.sh || true
