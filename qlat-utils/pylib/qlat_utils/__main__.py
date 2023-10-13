import sys

usage_message = """
Usage:
    python3 -m qlat_utils crc32 ...
    python3 -m qlat_utils lat-io-glimpse ...
    python3 -m qlat_utils lat-io-diff ...
    python3 -m qlat_utils pickle-glimpse ...
    python3 -m qlat_utils qar-glimpse ...
    python3 -m qlat_utils qar ...
    python3 -m qlat_utils qlat-utils-config ...
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

action = sys.argv[1]
sys.argv = sys.argv[1:]

if action == "qlat-utils-config":
    from .scripts import qlat_utils_config
elif action == "crc32":
    from .scripts import crc32
elif action == "lat-io-glimpse":
    from .scripts import lat_io_glimpse
elif action == "lat-io-diff":
    from .scripts import lat_io_diff
elif action == "pickle-glimpse":
    from .scripts import pickle_glimpse
elif action == "qar-glimpse":
    from .scripts import qar_glimpse
elif action == "qar":
    from .scripts import qar
else:
    print(usage_message)
    exit()
