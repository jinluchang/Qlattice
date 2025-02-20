#!/usr/bin/env bash

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to run the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

if [ -z "$name" ] ; then
    name=""
    echo "You can export name=-cuda to enable CUDA support."
else
    echo "name='$name'"
    echo "Trying to build 'pkgs$name.qlat-jhub-env'."
fi

py_kernel_name=py-local$name

src="$script_path"
dst="$HOME/.local/share/jupyter/kernels/nix-build-$py_kernel_name"
mkdir -p "$dst"
time nix-build "$src"/q-pkgs.nix -A pkgs$name.qlat-jhub-env -o "$dst/result" "$@"
if [ ! -e "$dst"/result/bin/python3 ] ; then
    rmdir "$dst"
fi
ls -l "$dst"
"$dst"/result/bin/python3 -m ipykernel \
    install --user \
    --env "PATH" "$dst/result/bin:/run/current-system/sw/bin" \
    --env "PYTHONPATH" "" \
    --env "LD_LIBRARY_PATH" "$dst/result/lib" \
    --env "LIBRARY_PATH" "$dst/result/lib" \
    --env "CPATH" "$dst/result/include" \
    --env "PKG_CONFIG_PATH" "$dst/result/lib/pkgconfig" \
    --env "CUBACORES" "0" \
    --env "OMP_NUM_THREADS" "2" \
    --name=$py_kernel_name
