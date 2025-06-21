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

time (

for nix_version in "25.05" "24.11" ; do
    time nix-build "$src"/q-pkgs.nix -A qlat-name-list-file -o result-qlat-name-"$nix_version" --argstr version "$nix_version" "$@"
    echo
    echo result-qlat-name-"$nix_version"
    echo
    cat result-qlat-name-"$nix_version"
    yes q-pkgs | head -n "$(cat result-qlat-name-"$nix_version" | wc -l)" | paste -d '' - result-qlat-name-"$nix_version" >q-pkgs-list-"$nix_version".txt
    echo
    for name in $(cat q-pkgs-list-"$nix_version".txt) ; do
        echo
        echo "Building $nix_version $name"
        echo
        time nix-build "$src"/q-pkgs.nix -A "$name".qlat-tests -A "$name".qlat-env -o result-"$nix_version-$name" --argstr version "$nix_version" "$@"
    done
done

)
