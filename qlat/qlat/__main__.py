import sys

usage_message = """
Usage:
    python3 -m qlat qlat-config ...
    python3 -m qlat eigen-system-checksum ...
    python3 -m qlat eigen-system-repartition ...
    python3 -m qlat fields-checksum ...
    python3 -m qlat fields-list ...
    python3 -m qlat fields-build-index ...
    python3 -m qlat fields-rewrite ...
    python3 -m qlat fields-properly-truncate ...
    python3 -m qlat gauge-fix-coulomb ...
    python3 -m qlat topo-measure ...
    python3 -m qlat topo-measure-wilson-flow ...
    python3 -m qlat topo-measure-gpt ...
""".strip()

if len(sys.argv) < 2:
    print(usage_message)
    exit()

action = sys.argv[1]
sys.argv = sys.argv[1:]

if action == "qlat-config":
    import qlat_config.__main__
elif action == "eigen-system-checksum":
    from .scripts import eigen_system_checksum
elif action == "eigen-system-repartition":
    from .scripts import eigen_system_repartition
elif action == "fields-checksum":
    from .scripts import fields_checksum
elif action == "fields-list":
    from .scripts import fields_list
elif action == "fields-build-index":
    from .scripts import fields_build_index
elif action == "fields-rewrite":
    from .scripts import fields_rewrite
elif action == "fields-properly-truncate":
    from .scripts import fields_properly_truncate
elif action == "gauge-fix-coulomb":
    from .scripts import gauge_fix_coulomb
elif action == "topo-measure":
    from .scripts import topo_measure
elif action == "topo-measure-wilson-flow":
    from .scripts import topo_measure_wilson_flow
elif action == "topo-measure-gpt":
    from .scripts import topo_measure_gpt
else:
    print(usage_message)
    exit()
