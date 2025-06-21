#!/usr/bin/env bash

cat <<EOF
Usage:

    name='' ./install-py-local-kernel-with-nix.sh
    # Install a local kernel with Qlattice and Grid/GPT.

    name='-cu' ./install-py-local-kernel-with-nix.sh
    # Install a local kernel with some CUDA utilities. Qlattice and Grid/GPT are NOT compiled with CUDA.

    name='-cuda' ./install-py-local-kernel-with-nix.sh
    # Install a local kernel with some CUDA utilities. Qlattice and Grid/GPT are compiled with CUDA.

    name='-cudasupport' ./install-py-local-kernel-with-nix.sh
    # Install a local kernel with full CUDA support when possible.

EOF

script_path="$( cd -- "$(dirname -- "$0")" >/dev/null 2>&1 ; pwd -P )"

if [ -f "$script_path"/q-pkgs.nix ]; then
    :
else
    echo "Need to place the script inside its original directory, which also have the file q-pkgs.nix and other files."
    exit 1
fi

set -ex

if [ -z "${name+x}" ] ; then
    echo "name is not set."
    echo "name='' (use default)."
else
    echo "name='$name'."
fi
echo "Trying to build 'pkgs$name.qlat-jhub-env'."

py_kernel_name=py-local$name

src="$script_path"
dst="$HOME/.local/share/jupyter/kernels/nix-build-$py_kernel_name"
mkdir -p "$dst"
time nix-build "$src"/q-pkgs.nix -A pkgs$name.qlat-jhub-env -o "$dst/result" "$@"
if [ ! -e "$dst"/result/bin/python3 ] ; then
    rmdir "$dst"
    exit 1
fi
ls -l "$dst"
"$dst"/result/bin/python3 -m ipykernel \
    install --user \
    --env "SHELL" "$dst/result/bin/bash" \
    --env "LOCALE_ARCHIVE" "/run/current-system/sw/lib/locale/locale-archive" \
    --env "PATH" "$dst/result/bin:/run/current-system/sw/bin" \
    --env "PKG_CONFIG_PATH" "$dst/result/lib/pkgconfig" \
    --env "LD_LIBRARY_PATH" "/run/opengl-driver/lib:$dst/result/lib" \
    --env "LIBRARY_PATH" "$dst/result/lib" \
    --env "CPATH" "$dst/result/include" \
    --env "PYTHONPATH" "" \
    --env "CUBACORES" "0" \
    --env "OMP_NUM_THREADS" "2" \
    --env "JAX_ENABLE_X64" "True" \
    --name=$py_kernel_name
