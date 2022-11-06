#!/bin/bash

. scripts/conf.sh

mkdir -p $distfiles

cd $distfiles

dget() {
    name="$1"
    url="$2"
    if [ -f "$name" ] ; then
        echo "$name is downloaded"
    else
        wget --no-check-certificate -O "$name" -c "$url"
    fi
}

dget "c-lime.tar.gz" "https://github.com/usqcd-software/c-lime/tarball/master"

dget "hdf5-1.10.7.tar.bz2" "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.bz2"

if [ -d Grid-lehner ] ; then
    ( cd Grid-lehner ; git pull https://github.com/jinluchang/Grid.git )
else
    git clone https://github.com/jinluchang/Grid.git Grid-lehner
fi

if [ -e Grid-lehner/configure ] ; then
    echo "Grid-lehner bootstrapped."
else
    ( cd Grid-lehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d gpt ] ; then
    ( cd gpt ; git pull https://github.com/jinluchang/gpt.git )
else
    git clone https://github.com/jinluchang/gpt.git
fi

if [ -d Grid-tblum ] ; then
    ( cd Grid-tblum ; git pull )
else
    git clone https://github.com/tblum2/Grid.git Grid-tblum
fi

if [ -e Grid-tblum/configure ] ; then
    echo "Grid-tblum bootstrapped."
else
    ( cd Grid-tblum ; git clean -f ; ./bootstrap.sh )
fi

if [ -d Hadrons-tblum ] ; then
    ( cd Hadrons-tblum ; git pull )
else
    git clone https://github.com/tblum2/Hadrons.git Hadrons-tblum
fi

if [ -e Hadrons-tblum/configure ] ; then
    echo "Hadrons-tblum bootstrapped."
else
    ( cd Hadrons-tblum ; git clean -f ; ./bootstrap.sh )
fi

cd $wd
