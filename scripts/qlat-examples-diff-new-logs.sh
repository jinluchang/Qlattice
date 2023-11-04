#!/bin/bash

source qcore/set-prefix.sh

source qcore/conf.sh

prefix_main="$prefix"

for prefix in "$prefix_main/qlat-examples-cpp" "$prefix_main/qlat-examples-cpp-grid" ; do
    for log in examples-cpp/*/log ; do
        echo diff "$prefix/$log" "$log"
        diff "$prefix/$log" "$log" | grep 'CHECK: ' && ( echo "$log" ; cat "$prefix/$log" || true )
        if diff "$prefix/$log" "$log" >/dev/null 2>&1 ; then
            :
        else
            cp -rpv "$prefix/$log" "$log".new || true
            cp -rpv "$prefix/$log.full" "$log".full.txt || true
        fi
    done
done

for prefix in "$prefix_main/qlat-examples-py" "$prefix_main/qlat-examples-py-gpt" "$prefix_main/qlat-examples-py-cps" ; do
    for log in examples-py/*.log ; do
        echo diff "$prefix/$log" "$log"
        diff "$prefix/$log" "$log" | grep 'CHECK: ' && ( echo "$log" ; cat "$prefix/$log" || true )
        if diff "$prefix/$log" "$log" >/dev/null 2>&1 ; then
            :
        else
            cp -rpv "$prefix/$log" "$log".new || true
            cp -rpv "$prefix/${log%.log}.py.p/log.full.txt" "$log".full.txt || true
        fi
    done
done

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" || true
        if diff "$log".new "$log" | grep 'CHECK: ' ; then
            cat "$log".full.txt || true
            echo diff "$log".new "$log"
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
