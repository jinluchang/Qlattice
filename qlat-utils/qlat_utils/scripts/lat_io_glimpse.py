import argparse
import sys

import qlat_utils as q

def parse_args():
    parser = argparse.ArgumentParser(description="Display contents of lat data files.")
    parser.add_argument(
        "filenames",
        nargs="+",
        help="Lat data files to display",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    ld = q.LatData()
    for fn in sys_args.filenames:
        ld.load(fn)
        sys.stdout.write(f"# '{fn}'\n")
        sys.stdout.write(ld.show())
    sys.exit()
