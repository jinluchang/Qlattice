#!/bin/bash

source qcore/conf.sh

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        rm -v "$log".new
    fi
done
