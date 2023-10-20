#!/bin/bash

source qcore/conf.sh

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" || true
        if diff "$log".new "$log" | grep 'CHECK: ' ; then
            cat "$log".full.txt || true
        fi
    fi
done

echo
echo "Summary:"
echo

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" | grep 'CHECK: ' || true
    fi
done
