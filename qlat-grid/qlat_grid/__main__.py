import sys

usage_message = """
Usage:
    python3 -m qlat_grid qlat-grid-config ...
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

action = sys.argv[1]
sys.argv = sys.argv[1:]

if action == "qlat-grid-config":
    from .scripts import qlat_grid_config
else:
    print(usage_message)
    exit()
