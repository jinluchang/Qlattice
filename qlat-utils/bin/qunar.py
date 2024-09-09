#!/usr/bin/env python3
#
# Standalone unarchiver for https://jinluchang.github.io/Qlattice/contents/qar-format.html
# Author: Christoph Lehner 2024
#
import sys, os

if len(sys.argv) == 4:
    src = sys.argv[1]
    tag = sys.argv[2]
    dst = sys.argv[3]
    os.makedirs(dst)

elif len(sys.argv) == 2:
    src = sys.argv[1]
    tag = None
    dst = None

else:
    print("qunar qar-file [dataset output-directory]")
    sys.exit(1)

f = open(src + "/index.qar", "rb")
header = f.readline().decode("utf-8")
header_expected = "#!/usr/bin/env qar-glimpse\n"

if header != header_expected:
    sys.stderr.write("Invalid qar file\n")
    sys.exit(2)

def expect_new_lines(n):
    for i in range(n):
        x = f.readline().decode("utf-8")
        if x != "\n":
            sys.stderr.write(f"Expected new line, got {x}")
            sys.exit(3)

def collapse_whitespace(x):
    return ' '.join(x.split())

expect_new_lines(1)

while True:
    ln = f.readline().decode("utf-8")
    if ln == "":
        break

    tags = collapse_whitespace(ln).split(" ")
    assert tags[0] == "QAR-FILE"
    sizes = [int(x) for x in tags[1:]]
    assert len(sizes) == 3
    name = f.read(sizes[0]).decode("utf-8")
    expect_new_lines(1)
    info = f.read(sizes[1])
    expect_new_lines(1)
    data = f.read(sizes[2])

    metadata = data.decode("utf-8").split("\n")[0:-1]
    ttag = metadata[1]

    if tag is None:
        print(f"{name}: {ttag}")
    elif tag == name:
        start = [int(x) for x in metadata[-2].split(" ")]
        end = [int(x) for x in metadata[-1].split(" ")]
        assert len(start) == len(end)
        ntotal = len(start)
        for rank in range(ntotal):
            dirs = 32
            nperdir = ntotal // dirs
            if nperdir < 1:
                nperdir = 1
            dirrank = rank // nperdir
            directory = "%2.2d" % (dirrank)
            filename = "%s/%10.10d" % (directory, rank)

            g = open(f"{src}/{filename}", "rb")
            g.seek(start[rank])
            rank_data = g.read(end[rank] - start[rank])
            g.close()

            os.makedirs(f"{dst}/{directory}", exist_ok=True)
            g = open(f"{dst}/{filename}", "wb")
            g.write(rank_data)
            g.close()

            print(f"Wrote {dst}/{filename}")
        sys.exit(0)

    expect_new_lines(2)




