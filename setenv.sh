if [ -z "$prefix" ] ; then
    prefix=$HOME/qlat-build/default
fi

if [ "$(uname -m)" = x86_64 ] ; then
    arch=amd64
else
    arch=unknown
fi

add-to-colon-list () {
    local name="$1"
    local new_value="$2"
    local value="${!name}"
    if [ -z "$value" ] ; then
        export "$name"="$new_value"
    else
        IFS=':' read -a vs <<< "$value"
        local value=''
        for v in "${vs[@]}" ; do
            if [ "$new_value" != "$v" ] ; then
                value+=:"$v"
            fi
        done
        export "$name"="$new_value""$value"
    fi
}

add-to-colon-list PYTHONPATH "$prefix/pylib"

