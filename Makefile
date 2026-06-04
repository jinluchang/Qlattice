all:

clean:
	make -C examples-py clean
	make -C examples-py-gpt clean
	make -C examples-py-cps clean
	make -C examples-cpp clean
	make -C examples-cpp-grid clean
