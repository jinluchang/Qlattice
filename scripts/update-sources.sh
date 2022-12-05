#!/bin/bash

wd=$(pwd)

( cd "$wd/qlat-utils/pylib/cqlat_utils" ; bash update.sh )

( cd "$wd/qlat/pylib/cqlat" ; bash update.sh )

( cd "$wd/qlat-grid/pylib/cqlat_grid" ; bash update.sh )
