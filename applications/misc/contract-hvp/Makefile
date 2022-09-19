run:
	mpirun --np 2 ./run.py --mpi 1.1.1.2
	find results -type f | sort | xargs -d '\n' md5sum > md5sum.txt

clean:
	-rm -rf locks results
