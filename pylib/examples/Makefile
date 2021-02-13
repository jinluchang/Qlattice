all: run

run: py-examples cplat-example

tests = make-sample-gauge-field.log \
		hmc-pure-gauge.log

py-examples : $(tests)

%.log : %.py
	mpirun --np 8 ./"$<" | tee "$@".full
	-rm -rfv results
	-grep -v "^Timer" "$@".full > "$@"

cplat-example:
	make -C cplat -j 2
	python3 -c "import cplat ; cplat.hello_world()" > cplat.log

clean:
	-rm -v *.log.full
	-make -C cplat clean
	-rm -v cplat.so

clean-logs:
	-rm -v *.log
