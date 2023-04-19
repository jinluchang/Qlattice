#!/bin/bash

time {
./scripts/qlat-utils.sh
./scripts/qlat.sh
./scripts/qlat-grid.sh
./scripts/qlat-cps.sh
}
