#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/qlat-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file qlat-pkgs.nix and other files."
    exit 1
fi

src="$script_path"
dst="$HOME/qlat-build"
mkdir -p "$dst"
cd "$dst"
for name in
    ""
    "-ucxless"
    "-std-clang"
    "-std-cubaquadless"
    "-cuda"
    "-cuda-ucxless"
    "-pypi"
    ; do
    time nix-build "$src"/qlat-pkgs.nix -A qlat-jhub-tests"$name" -o result"$name" "$@"
done
