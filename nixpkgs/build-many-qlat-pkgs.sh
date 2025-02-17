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
for name in \
    "" \
    "-cu" \
    "-cuda" \
    "-cudasupport" \
    "-std" \
    "-gridless" \
    "-cpsless" \
    "-gridless-cubaquadless" \
    "-std-clang" \
    "-ucxless" \
    "-pypi" \
    "-std-ucxless" \
    "-std-clang-ucxless" \
    "-cu-ucxless" \
    "-cuda-ucxless" \
    "-std-cu" \
    "-std-cuda" \
    "-std-cu-ucxless" \
    "-std-cuda-ucxless" \
    "-pypi" \
    ; do
    time nix-build "$src"/q-pkgs.nix -A qlat-jhub-tests"$name" -o result-24-11"$name" --arg nixpkgs 'import (fetchTarball "https://channels.nixos.org/nixos-24.11/nixexprs.tar.xz")' "$@"
done
for name in \
    "" \
    "-cu" \
    "-cuda" \
    "-cudasupport" \
    "-std" \
    "-gridless" \
    "-cpsless" \
    "-gridless-cubaquadless" \
    "-std-clang" \
    "-ucxless" \
    "-pypi" \
    "-std-ucxless" \
    "-std-clang-ucxless" \
    "-cu-ucxless" \
    "-cuda-ucxless" \
    "-std-cu" \
    "-std-cuda" \
    "-std-cu-ucxless" \
    "-std-cuda-ucxless" \
    "-pypi" \
    ; do
    time nix-build "$src"/q-pkgs.nix -A qlat-jhub-tests"$name" -o result-24-05"$name" --arg nixpkgs 'import (fetchTarball "https://channels.nixos.org/nixos-24.05/nixexprs.tar.xz")' "$@"
done
)
