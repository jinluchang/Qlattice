import argparse
import sys
import zlib

def crc32(fileName):
    with open(fileName, "rb") as fh:
        hash = 0
        while True:
            s = fh.read(65536)
            if not s:
                break
            hash = zlib.crc32(s, hash)
        return "%08x" % (hash & 0xFFFFFFFF)

def parse_args():
    parser = argparse.ArgumentParser(description="Compute CRC32 checksums for files.")
    parser.add_argument(
        "filenames",
        nargs="+",
        help="Files to compute checksums for",
    )
    args, _ = parser.parse_known_args()
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    for v in sys_args.filenames:
        print(f"{crc32(v)} '{v}'")
    sys.exit()
