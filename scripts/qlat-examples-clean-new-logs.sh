#!/usr/bin/env bash

source qcore/conf.sh

for log in examples-py*/*.log examples-cpp*/*/log ; do
    if [ -f "$log".new ] ; then
        rm -v "$log".new
    fi
    if [ -f "$log".json.new ] ; then
        rm -v "$log".json.new
    fi
    if [ -f "$log".full.txt ] ; then
        rm -v "$log".full.txt
    fi
done
