# Author: Luchang Jin 2022

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: fields-list path1 path2 ...")
    exit()

q.begin_with_mpi()

for path in sys.argv[1:]:
    q.displayln_info(0, path)
    sfr = q.open_fields(path, "r")
    tags = sfr.list()
    has_dup = sfr.has_duplicates()
    sfr.close()
    for tag in tags:
        q.displayln_info(tag)
    if has_dup:
        q.displayln_info(0, f"INFO: '{path}' has_duplicates.")

q.timer_display()

q.end_with_mpi()

exit()
