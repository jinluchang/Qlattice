#!/usr/bin/env bash

source qcore/setenv.sh
source qcore/conf.sh

source "$wd"/scripts/core-versions.sh

mkdir -p "$distfiles"

cd "$distfiles"

dget() {
    name="$1"
    url="$2"
    if [ -f "$name" ] ; then
        echo "$name is downloaded"
    else
        wget --no-check-certificate -O "$name" -c "$url"
    fi
}

aget() {
    url="$1"
    name="${1##*/}"
    dget "$name" "$url"
}

aget "http://usqcd-software.github.io/downloads/c-lime/lime-1.3.2.tar.gz"

aget "http://usqcd-software.github.io/downloads/qmp/qmp-2.5.4.tar.gz"

aget "http://usqcd-software.github.io/downloads/qio/qio-3.0.0.tar.gz"

aget "https://downloads.sourceforge.net/libtirpc/libtirpc-1.3.4.tar.bz2"

aget "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2"

aget "https://fftw.org/pub/fftw/fftw-3.3.10.tar.gz"

# aget "https://feynarts.de/cuba/Cuba-4.2.2.tar.gz"
aget "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/Cuba-4.2.2.tar.gz"

# aget "https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-1.14.5.tar.gz"
aget "https://github.com/jinluchang/Qlattice-distfiles/raw/refs/heads/main/distfiles/hdf5-1.14.5.tar.gz"

### -------------------------------------------------------------------
# The git repositories are checked out at the exact commits recorded in
# scripts/core-versions.sh (no "latest"; update the pins with
# ./scripts/update-core-versions.sh --update).

pin_repo() {
    local name="$1"
    local entry=""
    local e
    for e in "${CORE_REPOS[@]}" ; do
        if [ "${e%%|*}" = "$name" ] ; then
            entry="$e"
            break
        fi
    done
    if [ -z "$entry" ] ; then
        echo "ERROR: '$name' is not pinned in 'scripts/core-versions.sh'." >&2
        exit 1
    fi
    local url="" branch="" rev="" date=""
    IFS='|' read -r name url branch rev date <<< "$entry"
    local dir="$distfiles/$name"
    if [ -d "$dir" ] ; then
        local origin=""
        origin="$(git -C "$dir" remote get-url origin 2>/dev/null || true)"
        if [ "$origin" != "$url" ] ; then
            echo "'$name': repoint origin from '$origin' to '$url'"
            git -C "$dir" remote set-url origin "$url"
        fi
        if [ "$(git -C "$dir" rev-parse HEAD 2>/dev/null || true)" = "$rev" ] ; then
            echo "'$name': already at $rev"
            return 0
        fi
        if ! git -C "$dir" cat-file -e "$rev^{commit}" 2>/dev/null ; then
            if ! git -C "$dir" diff --quiet 2>/dev/null ; then
                echo "'$name': WARNING: discarding local changes (distfiles/ is a cache)."
            fi
            echo "'$name': fetching '$url'"
            git -C "$dir" fetch --all --prune
        fi
    else
        echo "'$name': cloning '$url' (branch '$branch')"
        git clone "$url" "$dir"
    fi
    if ! git -C "$dir" cat-file -e "$rev^{commit}" 2>/dev/null ; then
        echo "'$name': commit '$rev' is not in the local clone; fetching it from origin"
        if ! git -C "$dir" fetch origin "$rev" ; then
            echo "ERROR: cannot fetch commit '$rev' of '$name' from '$url'." >&2
            echo "ERROR: check the '$name' entry in 'scripts/core-versions.sh'." >&2
            exit 1
        fi
    fi
    echo "'$name': checking out $rev"
    git -C "$dir" checkout --detach --force "$rev"
    local got=""
    got="$(git -C "$dir" rev-parse HEAD)"
    if [ "$got" != "$rev" ] ; then
        echo "ERROR: '$name' is at '$got' but 'scripts/core-versions.sh' pins '$rev'." >&2
        exit 1
    fi
    # 'configure' is generated (untracked); drop it so the bootstrap below
    # regenerates it for the newly checked out revision.
    rm -f "$dir/configure"
    echo "'$name': pinned at $got"
}

pin_repo Grid

if [ -e Grid/configure ] ; then
    echo "Grid bootstrapped."
else
    ( cd Grid ; git clean -f ; ./bootstrap.sh || true ; ls -l Eigen || true )
fi

pin_repo Grid-clehner

if [ -e Grid-clehner/configure ] ; then
    echo "Grid-clehner bootstrapped."
else
    ( cd Grid-clehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

pin_repo Grid-ljin

if [ -e Grid-ljin/configure ] ; then
    echo "Grid-ljin bootstrapped."
else
    ( cd Grid-ljin ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

pin_repo gpt

pin_repo gpt-ljin

pin_repo Hadrons

if [ -e Hadrons/configure ] ; then
    echo "Hadrons bootstrapped."
else
    ( cd Hadrons ; git clean -f ; ./bootstrap.sh || true )
fi

pin_repo Grid-tblum

if [ -e Grid-tblum/configure ] ; then
    echo "Grid-tblum bootstrapped."
else
    ( cd Grid-tblum ; git clean -f ; if command -v gsha256sum ; then echo OK ; else sed -i 's/gsha256sum/sha256sum/g' bootstrap.sh || true ; fi ; ./bootstrap.sh || true ; ls -l Eigen || true ; git checkout bootstrap.sh || true )
fi

pin_repo Hadrons-tblum

if [ -e Hadrons-tblum/configure ] ; then
    echo "Hadrons-tblum bootstrapped."
else
    ( cd Hadrons-tblum ; git clean -f ; ./bootstrap.sh || true )
fi

pin_repo Hadrons-vaishakhiN

if [ -e Hadrons-vaishakhiN/configure ] ; then
    echo "Hadrons-vaishakhiN bootstrapped."
else
    ( cd Hadrons-vaishakhiN ; git clean -f ; ./bootstrap.sh || true )
fi

pin_repo CPS

### -------------------------------------------------------------------
# Summary (also handy in CI logs).

echo "Pinned core repositories:"
for entry in "${CORE_REPOS[@]}" ; do
    IFS='|' read -r name url branch rev date <<< "$entry"
    if [ -d "$distfiles/$name/.git" ] ; then
        echo "  $name: $(git -C "$distfiles/$name" rev-parse HEAD) ($branch)"
    fi
done

cd "$wd"
