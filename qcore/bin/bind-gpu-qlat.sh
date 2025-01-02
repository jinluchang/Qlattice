#!/usr/bin/env bash

if [ -z "$NGPU" ] ; then
    if which nvidia-smi >/dev/null 2>&1 ; then
        NGPU="$(nvidia-smi -L | wc -l)"
    else
        NGPU=1
    fi
fi

if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ] ; then
    export CUDA_VISIBLE_DEVICES="$((OMPI_COMM_WORLD_LOCAL_RANK % NGPU))"
fi

"$@"
