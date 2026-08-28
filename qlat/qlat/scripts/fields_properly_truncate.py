# Author: Luchang Jin 2023

import argparse

import qlat as q

def parse_args():
    parser = argparse.ArgumentParser(
        description="Properly truncate fields and remove duplicates."
    )
    parser.add_argument(
        "--check_all",
        action="store_true",
        default=False,
        help="Check all fields",
    )
    parser.add_argument(
        "--only_check",
        action="store_true",
        default=False,
        help="Only check, do not truncate",
    )
    parser.add_argument(
        "path_list",
        nargs="+",
        help="Paths to fields to truncate",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()

    q.begin_with_mpi()

    for path in sys_args.path_list:
        tags = q.properly_truncate_fields(
            path, is_check_all=sys_args.check_all, is_only_check=sys_args.only_check
        )
        for tag in tags:
            q.displayln_info(tag)

    q.timer_display()

    q.end_with_mpi()
