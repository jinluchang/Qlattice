all: run

run: py-examples cplat-example

tests = make-sample-gauge-field.log \
		field-utils.log \
		hmc-pure-gauge.log \
		field-selection.log \
		selected-field.log \
		lat-io.log \
		fields-io.log \
		free-invert.log \
		gpt-qlat-lanc.log \
		gpt-qlat-free-invert.log \
		gpt-qlat-mdwf.log \
		gpt-qlat-madwf.log \
		gpt-qlat-gfix-mdwf.log \
		gpt-qlat-rbc-ukqcd-invert.log

py-examples : $(tests)

%.log : %.py
	mpirun --np 2 ./"$<" --mpi 1.1.1.2 --mpi_split 1.1.1.1 | tee "$@".full
	-rm -rfv results
	-grep -v "^Timer\|^GPT\|^Grid\|timing: \|Timing: " "$@".full > "$@"

cplat-example:
	make -C cplat -j 2
	python3 -c "import cplat ; cplat.hello_world()" > cplat.log

clean:
	-rm -v *.log.full
	-make -C cplat clean
	-rm -v cplat.so

clean-logs:
	-rm -v *.log
