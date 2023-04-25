#!/bin/bash

if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ] ; then
    export CUDA_VISIBLE_DEVICES="$OMPI_COMM_WORLD_LOCAL_RANK"
fi

"$@"
