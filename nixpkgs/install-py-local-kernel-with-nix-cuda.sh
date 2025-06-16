#!/usr/bin/env bash

script_path="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to place the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

time nix-build "$script_path"/q-pkgs.nix \
  -A pkgs.qlat-jhub-env \
  -A pkgs-cu.qlat-jhub-env \
  -A pkgs-cuda.qlat-jhub-env \
  -A pkgs-cudasupport.qlat-jhub-env \
  --no-out-link "$@"

time (
  time name="" "$script_path"/install-py-local-kernel-with-nix.sh "$@"
  time name="-cu" "$script_path"/install-py-local-kernel-with-nix.sh "$@"
  time name="-cuda" "$script_path"/install-py-local-kernel-with-nix.sh "$@"
  time name="-cudasupport" "$script_path"/install-py-local-kernel-with-nix.sh "$@"
)
