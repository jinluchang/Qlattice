import sys

usage_message = """
Usage:
    python3 -m qlat_cps qlat-cps-config ...
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

action = sys.argv[1]
sys.argv = sys.argv[1:]

if action == "qlat-cps-config":
    import qlat_cps_config.__main__
else:
    print(usage_message)
    exit()
