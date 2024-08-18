#!/usr/bin/env bash

set -e

{
    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/ncurses.sh
    ./scripts/libevent.sh

    ./scripts/tmux.sh

    date
} 2>&1 | tee $prefix/log.build.txt
