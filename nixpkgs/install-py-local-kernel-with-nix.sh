#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/default-jhub.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file default-jhub.nix and other files."
    exit 1
fi

src="$script_path"
dst="$HOME/.local/share/jupyter/kernels"
mkdir -p "$dst"
cd "$dst"
time nix-build "$src"/default-jhub.nix "$@"
ls -l
./result/bin/python3 -m ipykernel \
    install --user \
    --env "PATH" "$dst/result/bin" \
    --env "PYTHONPATH" "" \
    --env "LD_LIBRARY_PATH" "" \
    --env "LD_RUN_PATH" "" \
    --env "CPLUS_INCLUDE_PATH" "" \
    --env "C_INCLUDE_PATH" "" \
    --env "LIBRARY_PATH" "" \
    --env "PKG_CONFIG_PATH" "" \
    --env "CUBACORES" "0" \
    --env "OMP_NUM_THREADS" "2" \
    --name=py-local
