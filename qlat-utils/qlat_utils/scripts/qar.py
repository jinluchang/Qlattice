# Author: Luchang Jin 2024

import argparse
import sys

import qlat_utils as q

def remove_trailing_slashes(path):
    while True:
        if path == "":
            return ""
        if path[-1] == "/":
            path = path[:-1]
        else:
            break
    return path

def show_list_qar(path_qar, idx=0, is_recursive=True, drop_prefix=""):
    assert path_qar[-4:] == ".qar"
    fns = q.list_qar(path_qar)
    for i, fn in enumerate(fns):
        pfn_full = path_qar[:-4] + "/" + fn
        assert pfn_full.startswith(drop_prefix)
        pfn = pfn_full[len(drop_prefix) :]
        q.displayln_info(f"{idx:8} '{path_qar}' {i:8} '{pfn}'")
        idx += 1
        if is_recursive and fn[-4:] == ".qar":
            idx = show_list_qar(
                pfn_full, idx, is_recursive=is_recursive, drop_prefix=drop_prefix
            )
    return idx

def build_index_qar(path_qar, is_recursive=True):
    assert path_qar[-4:] == ".qar"
    q.qar_build_index_info(path_qar)
    if not is_recursive:
        return
    fns = q.list_qar(path_qar)
    for fn in fns:
        if fn[-4:] == ".qar":
            pfn_full = path_qar[:-4] + "/" + fn
            build_index_qar(pfn_full, is_recursive=is_recursive)

def parse_args():
    parser = argparse.ArgumentParser(description="Manage qar archive files.")
    subparsers = parser.add_subparsers(dest="action", help="Action to perform")
    #
    # list: exactly 1 path, non-recursive
    subparsers.add_parser(
        "list", help="List contents of a qar archive (non-recursive, single path)"
    )
    # l: multiple paths, non-recursive (cumulative index)
    subparsers.add_parser(
        "l", help="List contents of qar archives (non-recursive, multiple paths)"
    )
    # lr: multiple paths, recursive (cumulative index)
    subparsers.add_parser("lr", help="List contents of qar archives (recursive)")
    # build-idx: exactly 1 path, non-recursive
    subparsers.add_parser(
        "build-idx", help="Build index for a qar archive (non-recursive, single path)"
    )
    # b: multiple paths, non-recursive
    subparsers.add_parser(
        "b", help="Build index for qar archives (non-recursive, multiple paths)"
    )
    # br: multiple paths, recursive
    subparsers.add_parser("br", help="Build index for qar archives (recursive)")
    # create: 2 explicit args (qar_path, dir_path)
    subparsers.add_parser(
        "create", help="Create a qar archive from a directory (explicit qar path + dir)"
    )
    # c: multiple dirs, auto-appends .qar
    subparsers.add_parser(
        "c", help="Create qar archives from directories (auto-appends .qar)"
    )
    # cr: multiple dirs, auto-appends .qar, remove dirs after
    subparsers.add_parser(
        "cr", help="Create qar archives and remove source directories"
    )
    # extract: 2 explicit args (qar_path, dest_path)
    subparsers.add_parser(
        "extract", help="Extract a qar archive to a destination (explicit qar + dest)"
    )
    # x: multiple qars, auto-derives dest by stripping .qar
    subparsers.add_parser(
        "x", help="Extract qar archives (auto-derives destination from .qar name)"
    )
    # xr: multiple qars, auto-derives dest, remove qars after
    subparsers.add_parser("xr", help="Extract qar archives and remove archive files")
    # cp: 2 args (src, dst)
    subparsers.add_parser("cp", help="Copy a qar archive file")
    # cat: multiple qars, concatenate
    subparsers.add_parser("cat", help="Concatenate and display qar archive contents")
    #
    args, remaining = parser.parse_known_args()
    #
    if args.action is None:
        parser.print_help()
        sys.exit(1)
    #
    # Store remaining args for manual processing
    args.paths = remaining
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    action = sys_args.action
    paths = sys_args.paths

    if action == "list":
        assert len(paths) == 1, f"list requires exactly 1 path, got {len(paths)}"
        show_list_qar(paths[0], 0, is_recursive=False)
    elif action == "l":
        idx = 0
        for path_qar in paths:
            idx = show_list_qar(path_qar, idx, is_recursive=False)
    elif action == "lr":
        idx = 0
        for path_qar in paths:
            idx = show_list_qar(path_qar, idx, is_recursive=True)
    elif action == "build-idx":
        assert len(paths) == 1, f"build-idx requires exactly 1 path, got {len(paths)}"
        build_index_qar(paths[0], is_recursive=False)
    elif action == "b":
        for path_qar in paths:
            build_index_qar(path_qar, is_recursive=False)
    elif action == "br":
        for path_qar in paths:
            build_index_qar(path_qar, is_recursive=True)
    elif action == "create":
        assert len(paths) == 2, (
            f"create requires exactly 2 paths (qar_path, dir_path), got {len(paths)}"
        )
        path_qar = paths[0]
        path = paths[1]
        assert not q.does_file_exist(path_qar)
        assert q.is_directory(path)
        q.qar_create_info(path_qar, path)
    elif action == "c":
        path_list = paths
        for path in path_list:
            path = remove_trailing_slashes(path)
            path_qar = path + ".qar"
            assert not q.does_file_exist(path_qar)
            assert q.is_directory(path)
        for path in path_list:
            path = remove_trailing_slashes(path)
            path_qar = path + ".qar"
            q.qar_create_info(path_qar, path)
    elif action == "cr":
        path_list = paths
        for path in path_list:
            path = remove_trailing_slashes(path)
            path_qar = path + ".qar"
            assert not q.does_file_exist(path_qar)
            assert q.is_directory(path)
        for path in path_list:
            path = remove_trailing_slashes(path)
            path_qar = path + ".qar"
            q.qar_create_info(path_qar, path, is_remove_folder_after=True)
    elif action == "extract":
        assert len(paths) == 2, (
            f"extract requires exactly 2 paths (qar_path, dest_path), got {len(paths)}"
        )
        path_qar = paths[0]
        path = paths[1]
        assert q.does_file_exist_qar(path_qar)
        assert not q.does_file_exist(path)
        q.qar_extract_info(path_qar, path)
    elif action == "x":
        path_qar_list = paths
        for path_qar in path_qar_list:
            assert path_qar[-4:] == ".qar"
            path = path_qar[:-4]
            assert path != ""
            assert q.does_file_exist_qar(path_qar)
            assert not q.does_file_exist(path)
        for path_qar in path_qar_list:
            path = path_qar[:-4]
            q.qar_extract_info(path_qar, path)
    elif action == "xr":
        path_qar_list = paths
        for path_qar in path_qar_list:
            assert path_qar[-4:] == ".qar"
            path = path_qar[:-4]
            assert path != ""
            assert q.does_file_exist_qar(path_qar)
            assert not q.does_file_exist(path)
        for path_qar in path_qar_list:
            path = path_qar[:-4]
            q.qar_extract_info(path_qar, path, is_remove_qar_after=True)
    elif action == "cp":
        assert len(paths) == 2, (
            f"cp requires exactly 2 paths (src, dst), got {len(paths)}"
        )
        path_src = paths[0]
        path_dst = paths[1]
        assert q.does_file_exist_qar(path_src)
        assert not q.does_file_exist(path_dst)
        q.qcopy_file_info(path_src, path_dst)
    elif action == "cat":
        for path in paths:
            assert q.does_file_exist_qar(path)
            content = q.qcat_bytes(path)
            sys.stdout.buffer.write(content)
    else:
        assert False, f"Unknown action: {action}"

    q.clear_all_caches()

    sys.exit()
