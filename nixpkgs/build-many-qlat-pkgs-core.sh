#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

set -ex

src="$script_path"
dst="$HOME/qlat-build/nix"
mkdir -p "$dst"
cd "$dst"

nix_version_list=(
    # "25.11"
    "25.05"
    # "24.11"
)

name_list=(
    q-pkgs
    # q-pkgs-clang
    # q-pkgs-ucxless
    q-pkgs-pypi
)

time (

    for nix_version in "${nix_version_list[@]}" ; do
        build_list=()
        for name in "${name_list[@]}" ; do
            build_list+=("-A" "$name.qlat-tests")
            build_list+=("-A" "$name.qlat-env")
        done
        if time nix-build "$src"/q-pkgs.nix "${build_list[@]}" --no-out-link --argstr version "$nix_version" "$@" ; then
            echo "Build successful."
        else
            echo "Build failed."
            exit 1
        fi
        for name in "${name_list[@]}" ; do
            echo
            echo "Building $nix_version $name"
            echo
            time nix-build "$src"/q-pkgs.nix -A "$name".qlat-tests -A "$name".qlat-env -o result-"$nix_version-$name" --argstr version "$nix_version" "$@" &
        done
        time wait
    done

)
