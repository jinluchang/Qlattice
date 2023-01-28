#!/bin/bash

# Usage:
# "$wd"/qcore/bin/mk-setenv-dir.sh
# "$wd"/qcore/bin/mk-setenv-dir.sh --keep

if ! [ -f "$wd/qcore/conf.sh" -a -f "$wd/qcore/set-prefix.sh" ] ; then
    echo "Need to set 'wd' to be the root directory of the repository."
    echo "Currently, wd='$wd'"
    return 1
fi

if "$wd"/qcore/bin/python-check-version.py >/dev/null 2>&1 && "$wd"/qcore/bin/mk-setenv-dir.py "$@" ; then
    # successfully finished.
    exit 0
fi

echo "Python mk-setenv-dir.py scripts failed. Fallback to bash version."

cat >"$prefix"/setenv-new.sh <<EOF
#!/bin/bash

setenv_prefix="$prefix"

EOF

if [ "--keep" = "$1" ] ; then
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
    cat "$prefix"/setenv.sh >>"$prefix"/setenv-new.sh
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
fi

cat >"$prefix"/setenv-new.sh <<EOF

for v in "\$setenv_prefix"/*/setenv.sh ; do
    if [ -f "\$v" ] ; then
        echo "Loading:" "\$v"
        source "\$v"
        echo "Loaded: " "\$v"
    fi
done
unset v

unset setenv_prefix

if which organize-env-path.sh >/dev/null 2>&1 ; then
    source <(organize-env-path.sh)
fi
EOF

mv -v "$prefix"/setenv-new.sh "$prefix"/setenv.sh
