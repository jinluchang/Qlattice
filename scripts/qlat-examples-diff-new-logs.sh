#!/bin/bash

./scripts/qlat-examples-clean-new-logs.sh

source qcore/set-prefix.sh

source qcore/conf.sh

prefix_main="$prefix"

for prefix in "$prefix_main/qlat-examples-cpp" "$prefix_main/qlat-examples-cpp-grid" ; do
    for log in examples-cpp/*/log ; do
        echo diff "$log" "$prefix/$log"
        if diff "$log" "$prefix/$log" >/dev/null 2>&1 ; then
            :
        else
            echo "$log differ"
            diff "$log" "$prefix/$log" || true
            cp -rpv "$prefix/$log" "$log".new || true
            cp -rpv "$prefix/$log.full" "$log".full.txt || true
        fi
    done
done

for prefix in "$prefix_main/qlat-examples-py" "$prefix_main/qlat-examples-py-gpt" "$prefix_main/qlat-examples-py-cps" ; do
    for log in examples-py/*.log ; do
        echo diff "$log" "$prefix/$log"
        if diff "$log" "$prefix/$log" >/dev/null 2>&1 ; then
            :
        else
            echo "$log differ"
            diff "$log" "$prefix/$log" || true
            cp -rpv "$prefix/$log" "$log".new || true
            cp -rpv "$prefix/$log".json "$log".json.new || true
            cp -rpv "$prefix/${log%.log}.py.p/log.full.txt" "$log".full.txt || true
        fi
    done
done

for log in examples-py/*.log examples-cpp/*/log ; do
    if [ -f "$log".new ] ; then
        echo diff "$log" "$log".new
        diff "$log" "$log".new || true
        if diff "$log" "$log".new | grep 'CHECK: ' ; then
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
