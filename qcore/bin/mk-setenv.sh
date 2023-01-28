#!/bin/bash

# Usage:
# "$wd"/qcore/bin/mk-setenv.sh
# "$wd"/qcore/bin/mk-setenv.sh --keep

if ! [ -f "$wd/qcore/conf.sh" -a -f "$wd/qcore/set-prefix.sh" ] ; then
    echo "Need to set 'wd' to be the root directory of the repository."
    echo "Currently, wd='$wd'"
    return 1
fi

if "$wd"/qcore/bin/python-check-version.py >/dev/null 2>&1 && "$wd"/qcore/bin/mk-setenv.py "$@" ; then
    # successfully finished.
    exit 0
fi

echo "Python mk-setenv.py scripts failed. Fallback to bash version."

cat >"$prefix"/setenv-new.sh <<EOF
#!/bin/bash

func() {

local setenv_prefix
local v

setenv_prefix="$prefix"

EOF

if [ "--keep" = "$1" ] ; then
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
    cat "$prefix"/setenv.sh >>"$prefix"/setenv-new.sh
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
fi

cat >>"$prefix"/setenv-new.sh <<EOF

for v in \
    "\$setenv_prefix/bin" \
    ; do
    if [ -d "\$v" ] ; then
        export PATH="\$v":"\$PATH"
    fi
done
for v in \
    "\$setenv_prefix/lib" \
    "\$setenv_prefix/lib64" \
    ; do
    if [ -d "\$v" ] ; then
        export LD_RUN_PATH="\$v":"\$LD_RUN_PATH"
        export LD_LIBRARY_PATH="\$v":"\$LD_LIBRARY_PATH"
        export LIBRARY_PATH="\$v":"\$LIBRARY_PATH"
    fi
done
for v in \
    "\$setenv_prefix/include" \
    ; do
    if [ -d "\$v" ] ; then
        export C_INCLUDE_PATH="\$v":"\$C_INCLUDE_PATH"
        export CPLUS_INCLUDE_PATH="\$v":"\$CPLUS_INCLUDE_PATH"
    fi
done
for v in \
    "\$setenv_prefix/lib/pkgconfig" \
    "\$setenv_prefix/lib64/pkgconfig" \
    ; do
    if [ -d "\$v" ] ; then
        export PKG_CONFIG_PATH="\$v":"\$PKG_CONFIG_PATH"
    fi
done
for v in \
    "\$setenv_prefix/lib/python3" \
    ; do
    if [ -d "\$v" ] ; then
        export PYTHONPATH="\$v":"\$PYTHONPATH"
    fi
done

}

func

if python-check-version.py >/dev/null 2>&1 && which organize-env-path.py >/dev/null 2>&1 ; then
    source <(organize-env-path.py)
fi
EOF

mv -v "$prefix"/setenv-new.sh "$prefix"/setenv.sh
