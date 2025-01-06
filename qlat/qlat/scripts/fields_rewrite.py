# Author: Luchang Jin 2023

import qlat as q
import sys

if len(sys.argv) < 2:
    q.displayln_info("Usage: fields-rewrite [--force] path1 path2 ...")
    q.displayln_info("Only rewrite if there is duplicates unless --force is given.")
    exit()

q.begin_with_mpi()

argv = sys.argv[1:]

is_force = q.get_option("--force", argv=argv, is_removing_from_argv=True)

path_list = argv

for path in path_list:
    path = q.remove_trailing_slashes(path)
    q.displayln_info(-1, f"INFO: Consider '{path}'.")
    sfr = q.open_fields(path, "r")
    has_duplicates = sfr.has_duplicates()
    if has_duplicates:
        q.displayln_info(-1, f"INFO: '{path}' has_duplicates.")
    do_rewrite = has_duplicates or is_force
    if do_rewrite:
        new_path = path + ".rewrite.acc"
        q.displayln_info(-1, f"INFO: '{path}' start to rewrite to '{new_path}'.")
        sfw = q.open_fields(new_path, "w", sfr.new_size_node())
        tags = sfr.list()
        for tag in tags:
            if tag in sfw:
                q.displayln_info(-1, f"INFO: Skip duplicated '{tag}' of '{sfr.path()}'.")
                continue
            q.displayln_info(-1, f"INFO: Read '{tag}' of '{sfr.path()}'.")
            obj = sfr.read_as_char(tag)
            q.displayln_info(-1, f"INFO: Write '{tag}' of '{sfw.path()}'.")
            obj.save_direct(sfw, tag)
        sfw.close()
    sfr.close()
    if do_rewrite:
        bak_path = path + ".bak"
        q.displayln_info(-1, f"INFO: Rename '{path}' to '{bak_path}'.")
        q.qrename_info(path, bak_path)
        q.displayln_info(-1, f"INFO: Rename '{new_path}' to '{path}'.")
        q.qrename_info(new_path, path)
    q.displayln_info(-1, f"INFO: Done '{path}'.")

q.timer_display()

q.end_with_mpi()

exit()
