if [ -z "$prefix" ] ; then
    prefix=$HOME/qlat-build/1.0
fi

if [ "$(uname -m)" = ppc64 ] ; then
    arch=bgq
elif [ "$(uname -m)" = x86_64 ] ; then
    arch=amd64
elif [ "$(uname -m)" = alpha ] ; then
    arch=taihu
else
    arch=unknown
fi

export HOME=~
export LANG=en_US.UTF-8

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

if [ $arch = bgq ] ; then
    add-to-colon-list PATH "/bgsys/drivers/ppcfloor/gnu-linux/bin"
    add-to-colon-list PATH "/bgsys/drivers/ppcfloor/comm/gcc/bin"
fi

add-to-colon-list PATH "$prefix/bin"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib"
add-to-colon-list LD_LIBRARY_PATH "$prefix/lib64"
add-to-colon-list LIBRARY_PATH "$prefix/lib"
add-to-colon-list LIBRARY_PATH "$prefix/lib64"
add-to-colon-list C_INCLUDE_PATH "$prefix/include/ncurses"
add-to-colon-list C_INCLUDE_PATH "$prefix/include"
add-to-colon-list CPLUS_INCLUDE_PATH "$prefix/include/ncurses"
add-to-colon-list CPLUS_INCLUDE_PATH "$prefix/include"
add-to-colon-list PKG_CONFIG_PATH "$prefix/lib/pkgconfig"

export MANPATH="$(man --path)"
add-to-colon-list MANPATH "$prefix/share/man"

if [ -z "$TEXINPUTS" ] ; then
    export TEXINPUTS=".:"
fi
add-to-colon-list TEXINPUTS "$prefix/share/texmf/tex/latex/local//"
add-to-colon-list TEXINPUTS "."
export TEXINPUTS="$TEXINPUTS":

