echo
echo "prefix=$prefix"
echo
echo "num_proc=$num_proc"

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

add-to-colon-list PATH "$prefix/bin"
add-to-colon-list PYTHONPATH "$prefix/pylib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib"
add-to-colon-list PYTHONPATH "$prefix/gpt/lib/cgpt/build"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib"

echo
echo "PATH=$PATH"
echo
echo "PYTHONPATH=$PYTHONPATH"
echo
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
