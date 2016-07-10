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
	-mkdir -p "$(prefix)/include/hash-cpp"
	cp *.h "$(prefix)/include/hash-cpp"

%.o: %.cpp *.h
	$(CC) -c -o $@ $<

lib: $(objects)
	ar rcs libhash-cpp.a *.o
	-mkdir -p "$(prefix)/lib"
	cp libhash-cpp.a "$(prefix)/lib"

clean:
	-rm *.o libhash-cpp.a
	-rm -rf local
