#!/bin/bash

{
    while : ; do
        if [ "$(find "$log" -mmin +15)" ] ; then
            echo Program Hang. Killing everything
        else
            sleep 180
            continue
        fi
        sync
        sleep 5
        mv "$log" failed-"$log"
        kill -9 -1
    done
} &
