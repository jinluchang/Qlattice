export OMP_NUM_THREADS=2

export q_verbose=2
export q_time_limit=$((12 * 60 * 60))
export q_mem_cache_max_size=512 # MB
export q_alloc_mem_max_size=$((256 * 1024)) # MB
export q_malloc_mmap_threshold=8192

grep -a "^CHECK: " log.txt >log.check.txt

timeout -s KILL 60m time mpiexec --np 4 $mpi_options python3 -m mpi4py ./run.py --test -qmp-geom 1 1 1 4 --mpi 1.1.1.4 >log.full.txt 2>&1

name="$(basename "$PWD")"

grep -a "^CHECK: \|^INFO: " log.full.txt >log.txt
grep -a "^CHECK: " log.txt >log.check.txt.new

if diff log.check.txt log.check.txt.new ; then
    echo passed
else
    tail -n 100 log.full.txt
    echo failed
fi
