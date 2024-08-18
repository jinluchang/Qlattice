#!/usr/bin/env bash

source qcore/conf.sh

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        mv -v "$log".new "$log"
    fi
    if [ -f "$log".json.new ] ; then
        mv -v "$log".json.new "$log".json
    fi
    if [ -f "$log".full.txt ] ; then
        rm -v "$log".full.txt
    fi
done
