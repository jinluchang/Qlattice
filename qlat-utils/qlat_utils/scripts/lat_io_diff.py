import argparse
import sys

import qlat_utils as q

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two lat data files and show differences."
    )
    parser.add_argument(
        "filenames",
        nargs=2,
        help="Two lat data files to compare",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    filenames = sys_args.filenames

    ld1 = q.load_lat_data(filenames[0])
    ld2 = q.load_lat_data(filenames[1])

    print(f"{q.qnorm(ld1 - ld2)}")
    print(f"{q.qnorm(ld1)}")
    print(f"{q.qnorm(ld2)}")

    sys.exit()
