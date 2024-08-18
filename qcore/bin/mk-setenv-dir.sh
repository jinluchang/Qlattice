#!/usr/bin/env bash

# Usage:
# "$wd"/qcore/bin/mk-setenv-dir.sh
# "$wd"/qcore/bin/mk-setenv-dir.sh --keep

if [ -z "$wd" ] ; then
    export wd="$PWD"
fi

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
#!/usr/bin/env bash

func() {

local setenv_prefix
local v

# setenv_prefix="$prefix"
setenv_prefix="\$(builtin cd -- "\$(dirname -- "\${BASH_SOURCE[0]}")" &> /dev/null && builtin pwd)"

setenv_path="\$(basename "\$setenv_prefix")"
export SETENV_PATH="\$setenv_path:\$SETENV_PATH"

EOF

if [ "--keep" = "$1" ] ; then
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
    cat "$prefix"/setenv.sh >>"$prefix"/setenv-new.sh
    echo "# ------------------------------------------------" >>"$prefix"/setenv-new.sh
fi

cat >>"$prefix"/setenv-new.sh <<EOF

if ! [ "--nr" = "\$1" ] ; then
    for v in "\$setenv_prefix"/*/setenv.sh ; do
        if [ -f "\$v" ] ; then
            echo "Loading:" "\$v"
            source "\$v"
            echo "Loaded: " "\$v"
        fi
    done
fi

}

func "\$@"

if python-check-version.py >/dev/null 2>&1 && which organize-env-path.py >/dev/null 2>&1 ; then
    eval "\$(organize-env-path.py)"
elif which organize-env-path.sh >/dev/null 2>&1 ; then
    eval "\$(organize-env-path.sh)"
fi
EOF

mv -v "$prefix"/setenv-new.sh "$prefix"/setenv.sh
