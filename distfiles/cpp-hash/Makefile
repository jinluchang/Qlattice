SHELL = /bin/bash

prefix = $(shell \
		 if [ -z "$$prefix" ] ; then \
		 echo local ; \
		 else echo "$$prefix" ; \
		 fi)

CC = g++ -Wall -O2

objects = crc32.o \
		  keccak.o \
		  md5.o \
		  sha1.o \
		  sha256.o \
		  sha3.o \

all: header lib

header:
	-mkdir -p "$(prefix)/include/cpp-hash"
	cp *.h "$(prefix)/include/cpp-hash"

%.o: %.cpp *.h
	$(CC) -c -o $@ $<

lib: $(objects)
	ar rcs cpp-hash.a *.o
	-mkdir -p "$(prefix)/lib"
	cp cpp-hash.a "$(prefix)/lib"

clean:
	-rm *.o cpp-hash.a
	-rm -rf local
