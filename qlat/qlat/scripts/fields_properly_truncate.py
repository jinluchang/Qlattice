# Author: Luchang Jin 2023

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: fields-properly-truncate [--check-all] [--only-check] path1 path2 ...")
    exit()

q.begin_with_mpi()

argv = sys.argv[1:]

is_check_all = q.get_option("--check-all", argv=argv, is_removing_from_argv=True)
is_only_check = q.get_option("--only-check", argv=argv, is_removing_from_argv=True)

path_list = argv

for path in path_list:
    tags = q.properly_truncate_fields(
            path, is_check_all = is_check_all, is_only_check = is_only_check)
    for tag in tags:
        q.displayln_info(tag)

q.timer_display()

q.end_with_mpi()

exit()
