#!/bin/bash

source qcore/conf.sh

(

cd "$wd/qlat/cqlat"

bash update.sh

) || echo "Run update.sh not successful."

(

version="$(cat VERSION)"
version=${version#v}

echo "Current version from file VERSION is: '$version'."

sed -i "s/^release = '.*'$/release = '$version'/" docs/conf.py

sed -i "s/^  version: '.*',$/  version: '$version',/" qlat*/meson.build

sed -i "s/^    version=\".*\" # default version$/    version=\"v$version-current\" # default version/" qlat-utils/qlat_utils/lib/version_gen.sh

echo "Version info updated."

) || echo "Update version not successful."

(

mkdir -p "$distfiles"

cd "$distfiles"

sha256sum *.tar.* | sort > sha256sums.txt

echo >> sha256sums.txt

sha256sum python-packages/*.* | sort >> sha256sums.txt

echo >> sha256sums.txt

for fn in * ; do
    if [ -e "$fn"/.git ] ; then
        echo -n "$fn: "
        ( cd "$fn" ; git rev-parse HEAD )
    fi
done | sort >> sha256sums.txt

cat sha256sums.txt

) || echo "Update sha256sums.txt not successful."
