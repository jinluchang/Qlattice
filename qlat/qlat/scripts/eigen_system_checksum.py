# Author: Luchang Jin 2022

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: eigen-system-checksum path1 path2 ...")
    exit(1)

q.begin_with_mpi()

for path in sys.argv[1:]:
    q.check_compressed_eigen_vectors(path)

q.timer_display()

q.end_with_mpi()

exit()
