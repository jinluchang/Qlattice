#!/usr/bin/env bash

add-to-colon-list() {
    local name="$1"
    local new_value="$2"
    local value="${!name}"
    local v
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
    local v
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

func() {
    local name
    local value
    for name in SETENV_PATH PATH PYTHONPATH LD_RUN_PATH LD_LIBRARY_PATH LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH PKG_CONFIG_PATH ; do
        organize-colon-list "$name"
        value="${!name}"
        echo "export $name=\"$value\""
    done
}

func
