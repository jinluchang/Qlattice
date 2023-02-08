#!/bin/bash

export
log="$PJM_STDOUT_PATH"
{
    export PS_FORMAT=user,pid,ppid,pgid,pmem,pcpu,etime,cputime,comm,cmd
    while : ; do
        if [ "$(find "$log" -mmin +5)" ] ; then
            echo Program Hang. Killing everything
            ps -u $(whoami)
            free -g
            uptime
        else
            {
                echo
                echo monitor
                uptime
                ps -u $(whoami)
                free -g
            } >> "$log"-monitor
            sleep 10
            continue
        fi
        sync
        sleep 5
        mv "$log" "$log"-failed
        kill -9 -1
    done
} &
