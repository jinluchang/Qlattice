#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

src="$script_path"
dst="$HOME/qlat-build"
mkdir -p "$dst"
cd "$dst"

time (

for nix_version in "24.11" "24.05" ; do
    time nix-build "$src"/q-pkgs.nix -A qlat-name-list-file -o result-qlat-name-"$nix_version" --arg nixpkgs "import (fetchTarball \"https://channels.nixos.org/nixos-$nix_version/nixexprs.tar.xz\")" "$@"
    echo
    echo result-qlat-name-"$nix_version"
    echo
    cat result-qlat-name-"$nix_version"
    echo
    for name in $(cat result-qlat-name-"$nix_version") ; do
        echo
        echo "Building $nix_version $name"
        echo
        time nix-build "$src"/q-pkgs.nix -A q-pkgs"$name" -o result-"$nix_version$name" --arg nixpkgs "import (fetchTarball \"https://channels.nixos.org/nixos-$nix_version/nixexprs.tar.xz\")" "$@"
    done
done

)
