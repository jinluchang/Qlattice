#!/usr/bin/env bash
#
# Build cuda packages (skip cuda tests).
# Output: ~/qlat-build/nix/cuda/result
# Names: ["", "-clang", "-cu", "-cudasupport", "-ucxless", "-pypi"]
# Extra args passed to nix-build via "$@".

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/many-qlat-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file many-qlat-pkgs.nix and other files."
    exit 1
fi

set -ex

src="$script_path"
dst="$HOME/qlat-build/nix"
mkdir -p "$dst"
cd "$dst"

time nix-build \
    "$src"/many-qlat-pkgs.nix \
    -o "$dst/cuda/result" \
    --arg qlat-name-list '["" "-clang" "-cu" "-cudasupport" "-ucxless" "-pypi"]' \
    --argstr qlat-cuda-tests none \
    --log-format internal-json -v \
    -j 6 --cores 15 \
    "$@" \
    |& nom --json
