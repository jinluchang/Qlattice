SHELL=/bin/bash

all: run

compile:
	meson build
	ninja -C build -j3

run: compile
	time ninja -C build run
	-if diff build/log log | grep 'CHECK: ' ; then cat build/log.full ; else echo passed ; fi
	-cp build/log log
	-cp build/log.full log.full

clean:
	-rm -rf build log.full
