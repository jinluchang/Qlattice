#!/bin/bash

source qcore/conf.sh

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        mv -v "$log".new "$log"
    fi
done
