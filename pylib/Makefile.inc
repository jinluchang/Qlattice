include $(QLAT_PREFIX)/include/qlat/Makefile.inc

CXXFLAGS+= $$(python3-config --includes)
# CXXFLAGS+= $$(python3 -c "import numpy; print(numpy.get_include())") 

NAME=cqlat.so

OBJS:=$(patsubst %.C, %.o, $(wildcard *.C))

all: build

build: $(OBJS) lib.o
	$(CXXLD) $(FLAGS) -o $(NAME) *.o base/*.o $(LDFLAGS)
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
	$(CXX) $(FLAGS) -c -o $@ $< $(CXXFLAGS) -Ibase
	@echo "$@" "$<"

clean:
	-rm -rfv *.o lib.C exports.h base
