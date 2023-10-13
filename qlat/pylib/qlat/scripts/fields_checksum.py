# Author: Luchang Jin 2022

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: fields-checksum path1 path2 ...")
    exit()

q.begin_with_mpi()

for path in sys.argv[1:]:
    q.check_fields(path)

q.timer_display()

q.end_with_mpi()

exit()
