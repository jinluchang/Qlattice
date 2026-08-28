import argparse
import pprint
import sys

import qlat_utils as q

def parse_args():
    parser = argparse.ArgumentParser(description="Display contents of pickle files.")
    parser.add_argument(
        "filenames",
        nargs="+",
        help="Pickle files to display",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    for fn in sys_args.filenames:
        obj = q.load_pickle_obj(fn)
        sys.stdout.write(f"# '{fn}'\n")
        sys.stdout.write(pprint.pformat(obj))
        sys.stdout.write("\n")
    sys.exit()
