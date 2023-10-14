# Author: Luchang Jin 2022

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: fields-list path1 path2 ...")
    exit()

q.begin_with_mpi()

for path in sys.argv[1:]:
    tags = q.list_fields(path)
    for tag in tags:
        q.displayln_info(tag)

q.timer_display()

q.end_with_mpi()

exit()
