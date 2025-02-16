#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

name=py-local

src="$script_path"
dst="$HOME/.local/share/jupyter/kernels/nix-build-$name"
mkdir -p "$dst"
cd "$dst"
time nix-build "$src"/q-pkgs.nix -A qlat-jhub-env "$@"
ls -l
./result/bin/python3 -m ipykernel \
    install --user \
    --env "PATH" "$dst/result/bin" \
    --env "PYTHONPATH" "" \
    --env "LD_LIBRARY_PATH" "$dst/result/lib" \
    --env "LD_RUN_PATH" "$dst/result/lib" \
    --env "CPLUS_INCLUDE_PATH" "$dst/result/include" \
    --env "C_INCLUDE_PATH" "$dst/result/include" \
    --env "LIBRARY_PATH" "$dst/result/lib" \
    --env "PKG_CONFIG_PATH" "$dst/result/lib/pkgconfig" \
    --env "CUBACORES" "0" \
    --env "OMP_NUM_THREADS" "2" \
    --name=$name
