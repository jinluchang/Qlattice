#!/usr/bin/env bash

source scripts/qlat-examples-clean-new-logs.sh

source qcore/set-prefix.sh

prefix_main="$prefix"

for name in "-cpp" "-cpp-grid" ; do
    prefix="$prefix_main/qlat-examples$name"
    for log in examples$name/*/log ; do
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

for name in "-py" "-py-gpt" "-py-cps" ; do
    prefix="$prefix_main/qlat-examples$name"
    for log in examples$name/*.log ; do
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
