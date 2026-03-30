#!/usr/bin/env bash

source qcore/conf.sh

(

echo "Update cqlat export."

cd "$wd/qlat/cqlat"

bash update.sh

) || echo "Update cqlat export not successful."

(

echo "Update version."

version="$(cat VERSION)"
version=${version#v}

echo "Current version from file VERSION is: '$version'."

sed -i "s/^release = '.*'$/release = '$version'/" docs/conf.py

sed -i "s/^  version: '.*',$/  version: '$version',/" qlat*/meson.build

sed -i "s/^    version=\".*\" # default version$/    version=\"v$version-current\" # default version/" qlat-utils/qlat_utils/lib/version_gen.sh

echo "Version info updated."

) || exit 1

(

echo "Update 'depend-qlat/meson.build'."

for dir in $(find qlat* examples-* -type d -name depend-qlat) ; do
    fn="$dir/meson.build"
    if [ -f "$fn" ] ; then
        cp -pv qlat-cps/depend-qlat/meson.build "$fn"
    fi
done

) || exit 1

(

echo "Update 'depend-grid/meson.build'."

for dir in $(find qlat* examples-* -type d -name depend-grid) ; do
    fn="$dir/meson.build"
    if [ -f "$fn" ] ; then
        cp -pv qlat-grid/depend-grid/meson.build "$fn"
    fi
done

) || exit 1

(

echo "Create '.clangd'."

fn=.clangd

echo "CompileFlags:" >"$fn"
echo "    Add: [" >>"$fn"
echo "        '-DQLAT_USE_GRID_EIGEN'," >>"$fn"
echo "        '-D__QLAT_NO_FLOAT128__'," >>"$fn"
echo "        '-DNO_CPS'," >>"$fn"
for pkg in zlib ompi python fftw3 fftw3f gsl ; do
    if pkg-config --cflags $pkg >/dev/null 2>&1 ; then
        echo "        '$(pkg-config --cflags $pkg)'," >>"$fn"
    else
        echo "pkg-config '$pkg' is not available."
    fi
done
if grid-config --cxxflags >/dev/null 2>&1 ; then
    for arg in $(grid-config --cxxflags) ; do
        if [ "${arg#-I}" == "${arg}" ] ; then
            echo "skip '$arg'"
        else
            echo "        '$arg'," >>"$fn"
        fi
    done
else
    echo "grid-config is not available."
fi
echo "        '-I$PWD/qlat-utils/qlat_utils/include'," >>"$fn"
echo "        '-I$PWD/qlat/qlat/include'," >>"$fn"
echo "        '-I$PWD/qlat-grid/qlat_grid/include'," >>"$fn"
echo "        '-I$PWD/qlat-cps/qlat_cps/include'" >>"$fn"
echo "    ]" >>"$fn"

) || echo "Create '.clangd' not successful."

(

find . -name "__pycache__" -exec rm -rfv '{}' \;

) || echo "Clean '__pycache__' not successful."

(

echo "Update 'sha256sums.txt'."

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

) || echo "Update 'sha256sums.txt' not successful."
