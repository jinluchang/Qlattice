SHELL=/bin/bash

all: run

compile:
	meson build
	ninja -C build

run: compile
	ninja -C build run
	cp build/log log

clean:
	-rm -rf build

clean-logs: clean
	-rm -rf log
