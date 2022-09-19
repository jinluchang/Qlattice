CXXFLAGS += -fPIC
LDFLAGS += --shared

include $(QLAT_PREFIX)/include/qlat/Makefile.inc

CXXFLAGS += $$(python3-config --includes)
# CXXFLAGS += $$(python3 -c "import numpy; print(numpy.get_include())")

NAME = cqlat.so

OBJS := $(patsubst %.C, %.o, $(wildcard *.C))

all: build

build: $(OBJS) lib.o
	$(MPICXX) $(LDFLAGS) -o $(NAME) *.o base/*.o $(LIBS)
	ls -lh $(NAME)
	mv $(NAME) ..

exports.h: *.C base/exports.h
	./base/update.sh base/exports.h

base/exports.h:
	rsync -av --delete $(QLAT_PREFIX)/pylib/cqlat/ base/

lib.C: base/exports.h
	rm -v base/lib.o
	cp -v base/lib.C .

%.o : %.C exports.h
	$(MPICXX) $(CXXFLAGS) -Ibase -c -o $@ $<
	@echo "$@" "$<"

clean:
	-rm -rfv *.o lib.C exports.h base
