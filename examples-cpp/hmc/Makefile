SHELL=/bin/bash

all: run

compile:
	meson build
	ninja -C build -j3

run: compile
	time ninja -C build run
	-cp build/log log
	-cp build/log.full log.full
	-diff log log.ref | grep 'CHECK: ' && cat log.full

clean:
	-rm -rf build

clean-logs: clean
	-mv log log.ref
