#!/usr/bin/env bash

source scripts/qlat-examples-clean-new-logs.sh

source qcore/set-prefix.sh

prefix_main="$prefix"

for prefix in "$prefix_main/qlat-examples-cpp" "$prefix_main/qlat-examples-cpp-grid" ; do
    for log in examples-cpp/*/log ; do
        # echo diff "$log" "$prefix/$log"
        if diff "$prefix/$log" "$log" >/dev/null 2>&1 ; then
            :
        else
            if [ -f "$prefix/$log" ] ; then
                cp -rpv "$prefix/$log" "$log".new || true
            fi
            if [ -f "$prefix/$log".full ] ; then
                cp -rpv "$prefix/$log".full "$log".full.txt || true
            fi
        fi
    done
done

for prefix in "$prefix_main/qlat-examples-py" "$prefix_main/qlat-examples-py-gpt" "$prefix_main/qlat-examples-py-cps" ; do
    for log in examples-py/*.log ; do
        # echo diff "$log" "$prefix/$log"
        if diff "$prefix/$log" "$log" >/dev/null 2>&1 ; then
            :
        else
            if [ -f "$prefix/$log" ] ; then
                cp -rpv "$prefix/$log" "$log".new || true
            fi
            if [ -f "$prefix/$log".json ] ; then
                cp -rpv "$prefix/$log".json "$log".json.new || true
            fi
            if [ -f "$prefix/${log%.log}.py.p/log.full.txt" ] ; then
                cp -rpv "$prefix/${log%.log}.py.p/log.full.txt" "$log".full.txt || true
            fi
        fi
    done
done

echo
echo "Summary:"
echo

for log in examples-py*/*.log examples-cpp*/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" | grep 'CHECK: ' || true
    fi
done

echo
echo "Summary for all diff:"
echo

for log in examples-py*/*.log examples-cpp*/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" || true
    fi
done

echo
echo "Full log:"
echo

for log in examples-py*/*.log examples-cpp*/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" || true
        if diff "$log".new "$log" | grep 'CHECK: ' ; then
            echo diff log start "$log".new "$log"
            cat "$log".full.txt || true
            echo diff log end "$log".new "$log"
        fi
    fi
done

echo
echo "Final Summary:"
echo

for log in examples-py*/*.log examples-cpp*/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log".new "$log"
        diff "$log".new "$log" | grep 'CHECK: ' || true
    fi
done
