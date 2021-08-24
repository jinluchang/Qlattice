#!/bin/bash

. conf.sh

mkdir -p $distfiles

cd $distfiles

sha256sum *.tar.* > sha256sums.txt

echo >> sha256sums.txt

echo -n "Grid: " >> sha256sums.txt
(
cd Grid
git rev-parse HEAD >> ../sha256sums.txt
)

echo -n "gpt: " >> sha256sums.txt
(
cd gpt
git rev-parse HEAD >> ../sha256sums.txt
)

cat sha256sums.txt
