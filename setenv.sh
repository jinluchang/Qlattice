echo
echo "prefix=$prefix"
echo
echo "num_proc=$num_proc"

export QLAT_PREFIX="$prefix"

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

organize-colon-list() {
    local name="$1"
    local value="${!name}"
    if [ -n "$value" ] ; then
        IFS=':' read -a vs <<< "$value"
        value=''
        for v in "${vs[@]}" ; do
            value="$v":"$value"
        done
        value="${value%:}"
        IFS=':' read -a vs <<< "$value"
        export "$name"=""
        for v in "${vs[@]}" ; do
            add-to-colon-list "$name" "$v"
        done
    fi
}

add-to-colon-list PATH "$prefix/bin"
add-to-colon-list PYTHONPATH "$prefix/pylib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib/cgpt/build"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib"

organize-colon-list PATH
organize-colon-list PYTHONPATH
organize-colon-list LD_LIBRARY_PATH

echo
for v in PATH PYTHONPATH LD_LIBRARY_PATH QLAT_PREFIX QLAT_CXX QLAT_CXXLD QLAT_FLAGS ; do
export | grep --color=never " $v="'"' || true
done
echo
