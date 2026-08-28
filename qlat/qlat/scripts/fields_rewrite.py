# Author: Luchang Jin 2023

import argparse

import qlat as q

def parse_args():
    parser = argparse.ArgumentParser(description="Rewrite fields to remove duplicates.")
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force rewrite even if no duplicates",
    )
    parser.add_argument(
        "path_list",
        nargs="+",
        help="Paths to fields to rewrite",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()

    q.begin_with_mpi()

    for path in sys_args.path_list:
        path = q.remove_trailing_slashes(path)
        q.displayln_info(-1, f"INFO: Consider '{path}'.")
        sfr = q.open_fields(path, "r")
        has_duplicates = sfr.has_duplicates()
        if has_duplicates:
            q.displayln_info(-1, f"INFO: '{path}' has_duplicates.")
        do_rewrite = has_duplicates or sys_args.force
        if do_rewrite:
            new_path = path + ".rewrite.acc"
            q.displayln_info(-1, f"INFO: '{path}' start to rewrite to '{new_path}'.")
            sfw = q.open_fields(new_path, "w", sfr.new_size_node())
            tags = sfr.list()
            for tag in tags:
                if tag in sfw:
                    q.displayln_info(
                        -1, f"INFO: Skip duplicated '{tag}' of '{sfr.path()}'."
                    )
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
