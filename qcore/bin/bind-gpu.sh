#!/usr/bin/env bash

if which nvidia-smi >/dev/null 2>&1 ; then
    NGPU="$(nvidia-smi -L | wc -l)"
fi

if [ -z "$NGPU" ] ; then
    NGPU=1
fi

if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ] ; then
    export CUDA_VISIBLE_DEVICES="$((OMPI_COMM_WORLD_LOCAL_RANK % NGPU))"
fi

"$@"
