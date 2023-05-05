#!/bin/bash

NGPU="$(nvidia-smi -L | wc -l)"

if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ] ; then
    export CUDA_VISIBLE_DEVICES="$((OMPI_COMM_WORLD_LOCAL_RANK % NGPU))"
fi

"$@"
