SHELL=/bin/bash

all: run

run: cqlat-example py-examples

tests = make-sample-gauge-field.log \
		field-utils.log \
		qcd-utils.log \
		hmc-pure-gauge.log \
		field-selection.log \
		selected-field.log \
		lat-io.log \
		fields-io.log \
		free-invert.log \
		gpt-qlat-convert.log \
		gpt-qlat-smear.log \
		gpt-qlat-sample-gauge-field.log \
		gpt-qlat-lanc.log \
		gpt-qlat-free-invert.log \
		gpt-qlat-mdwf.log \
		gpt-qlat-madwf.log \
		gpt-qlat-gfix-mdwf.log \
		gpt-qlat-rbc-ukqcd-invert.log

py-examples: $(tests)

%.log: %.py
	mpirun --np 2 ./"$<" --mpi 1.1.1.2 --mpi_split 1.1.1.1 --mpi 1.1.2 --mpi_split 1.1.1 | tee "$@".full
	-rm -rfv results
	-grep -v "^Timer\|^GPT\|^Grid\|timing: \|Timing: \|local time slice" "$@".full > "$@"

cqlat-example: cqlat.log

cqlat.log:
	-make -C cqlat -j 2
	python3 -c "import cqlat ; cqlat.hello_world()" > cqlat.log

clean:
	-rm -v *.log.full
	-make -C cqlat clean
	-rm -v cqlat.so

clean-logs:
	-rm -v *.log
