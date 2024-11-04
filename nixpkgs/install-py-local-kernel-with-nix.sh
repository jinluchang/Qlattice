# Need to run from this directory

src="$PWD"
dst="$HOME/.local/share/jupyter/kernels"
mkdir -p "$dst"
cd "$dst"
time nix-build "$src"/default-jhub.nix
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
