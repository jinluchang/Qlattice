#!/usr/bin/env python3

"""
Lint tool to enforce a coding style where empty lines inside function bodies
are replaced with indented '#' comment lines.\n
Design goals:
- Remove all empty lines within function/method bodies to enforce compact code.
- Use indented '#' comment lines as visual separators when logical grouping is needed.
- Preserve the indentation of the surrounding context for each replaced empty line.
- Leave empty lines outside function bodies untouched (module level, between functions).
- Remove empty lines inside multi-line string literals (e.g., docstrings) while preserving content by appending "\\n" to the preceding line.
- Handle nested functions correctly by scoping replacements to each function's body.
- Support both .py and .pyx (Cython) files.\n
Usage:
    python lint.py <file_or_directory> [--check] [--diff]\n
Options:
    --check   Check if files need changes (exit 1 if so)
    --diff    Show diff of proposed changes
"""

import ast
import sys
import os
import difflib
import argparse
import tokenize
import io
import subprocess
import tempfile
import re


def get_function_body_lines(source_lines):
    """Return set of line numbers that are inside a function/method body.\n
    Args:
        source_lines (list[str]): Source code lines.\n
    Returns:
        set[int]: 1-indexed line numbers inside function bodies.
    """
    source = "".join(source_lines)
    tree = ast.parse(source)
    body_lines = set()
    #
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.body:
                start = node.body[0].lineno
                end = node.end_lineno
                for line_no in range(start, end + 1):
                    body_lines.add(line_no)
    #
    return body_lines


def get_function_body_lines_cython(source_lines):
    """Return set of line numbers that are inside a function/method body for Cython files.\n
    Uses regex-based detection since ast.parse cannot handle .pyx syntax.\n
    Args:
        source_lines (list[str]): Source code lines.\n
    Returns:
        set[int]: 1-indexed line numbers inside function bodies.
    """
    body_lines = set()
    func_pattern = re.compile(
        r"^(\s*)(def|cdef|cpdef)\s+(?!class\b)(?:\w+\s+)*\w+\s*\("
    )
    i = 0
    while i < len(source_lines):
        m = func_pattern.match(source_lines[i])
        if m:
            base_indent = len(m.group(1))
            j = i + 1
            while j < len(source_lines) and source_lines[j].strip() == "":
                j += 1
            if j < len(source_lines):
                body_indent = len(get_indent(source_lines[j]))
                if body_indent > base_indent:
                    body_start = j
                    body_end = j
                    for k in range(j + 1, len(source_lines)):
                        line = source_lines[k]
                        if line.strip() == "":
                            continue
                        curr_indent = len(get_indent(line))
                        if curr_indent > base_indent:
                            body_end = k
                        else:
                            break
                    for ln in range(body_start, body_end + 1):
                        body_lines.add(ln + 1)
        i += 1
    return body_lines


def get_indent(line):
    """Return leading whitespace of a line.\n
    Args:
        line (str): Source code line.\n
    Returns:
        str: Leading whitespace string.
    """
    return line[: len(line) - len(line.lstrip())]


def get_string_lines(source):
    """Return set of line numbers that are inside multi-line string literals.\n
    Args:
        source (str): Python source code.\n
    Returns:
        set[int]: 1-indexed line numbers inside multi-line strings.
    """
    string_lines = set()
    tokens = tokenize.generate_tokens(io.StringIO(source).readline)
    for tok in tokens:
        if tok.type == tokenize.STRING:
            start, end = tok.start, tok.end
            if end[0] > start[0]:
                for ln in range(start[0], end[0]):
                    string_lines.add(ln)
    return string_lines


def collapse_consecutive_comments(lines):
    """Collapse multiple consecutive '#' comment lines into one.\n
    Args:
        lines (list[str]): Source code lines.\n
    Returns:
        list[str]: Lines with consecutive '#' comments collapsed.
    """
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped == "#":
            j = i + 1
            while j < len(lines) and lines[j].strip() == "#":
                j += 1
            if j > i + 1:
                result.append(line)
                i = j
                continue
        result.append(line)
        i += 1
    return result


def transform_file(source, is_cython=False):
    """Replace empty lines inside function bodies with properly indented '#' comment lines.\n
    Args:
        source (str): Python source code.
        is_cython (bool): If True, use regex-based function detection for .pyx files.\n
    Returns:
        str: Transformed source code.
    """
    lines = source.splitlines(keepends=True)
    if is_cython:
        body_lines = get_function_body_lines_cython(lines)
    else:
        body_lines = get_function_body_lines(lines)
    try:
        string_lines = get_string_lines(source)
    except tokenize.TokenError:
        string_lines = set()
    #
    # Build indent lookup: for each empty line in a function body,
    # use the indent of the nearest non-empty line (prefer next, then prev).
    # Exclude lines inside string literals to avoid modifying string content.
    indent_map = {}
    for i, line in enumerate(lines):
        line_no = i + 1
        if line_no in body_lines and line_no not in string_lines and line.strip() == "":
            indent = None
            # Search forward for next non-empty line
            for j in range(i + 1, len(lines)):
                if lines[j].strip() != "":
                    indent = get_indent(lines[j])
                    break
            # Fall back to previous non-empty line
            if indent is None:
                for j in range(i - 1, -1, -1):
                    if lines[j].strip() != "":
                        indent = get_indent(lines[j])
                        break
            if indent is not None:
                indent_map[line_no] = indent
    #
    # Identify empty lines in strings to remove (content preserved via trailing \n).
    remove_lines = set()
    for i, line in enumerate(lines):
        line_no = i + 1
        if line.strip() == "" and line_no in string_lines:
            remove_lines.add(line_no)
    #
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]
        line_no = i + 1
        if line_no in remove_lines:
            count = 1
            i += 1
            while i < len(lines) and (i + 1) in remove_lines:
                count += 1
                i += 1
            if result:
                prev = result[-1]
                if prev.endswith("\n"):
                    result[-1] = prev[:-1] + "\\n" * count + "\n"
                else:
                    result[-1] = prev + "\\n" * count
        elif line.strip() == "" and line_no in indent_map:
            indent = indent_map[line_no]
            if line.endswith("\n"):
                result.append(f"{indent}#\n")
            else:
                result.append(f"{indent}#")
            i += 1
        else:
            result.append(line)
            i += 1
    #
    result = collapse_consecutive_comments(result)
    return "".join(result)


def run_ruff(source):
    """Run ruff format and ruff check --fix on source code.\n
    Args:
        source (str): Python source code.\n
    Returns:
        str: Formatted source code after ruff processing.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        f.write(source)
        tmp = f.name
    try:
        subprocess.run(["ruff", "format", tmp], capture_output=True)
        subprocess.run(
            ["ruff", "check", "--fix", "--unsafe-fixes", tmp],
            capture_output=True,
        )
        with open(tmp, "r") as f:
            return f.read()
    finally:
        os.unlink(tmp)


def process_file(filepath, check=False, diff=False):
    """Process a single Python or Cython file.\n
    Args:
        filepath (str): Path to the file.
        check (bool): If True, only check if changes are needed.
        diff (bool): If True, show diff of changes.\n
    Returns:
        bool: True if file was changed (or would be changed in check mode).
    """
    is_cython = filepath.endswith(".pyx")
    with open(filepath, "r") as f:
        original = f.read()
    #
    if is_cython:
        transformed = transform_file(original, is_cython=True)
    else:
        ruffed = run_ruff(original)
        transformed = transform_file(ruffed)
    #
    if original == transformed:
        return False
    #
    if diff:
        d = difflib.unified_diff(
            original.splitlines(keepends=True),
            transformed.splitlines(keepends=True),
            fromfile=f"a/{filepath}",
            tofile=f"b/{filepath}",
        )
        sys.stdout.writelines(d)
    #
    if not check:
        with open(filepath, "w") as f:
            f.write(transformed)
    #
    return True


def find_source_files(path):
    """Find all Python and Cython files at path (file or directory).\n
    Args:
        path (str): File or directory path.\n
    Returns:
        list[str]: List of source file paths.
    """
    if os.path.isfile(path) and (path.endswith(".py") or path.endswith(".pyx")):
        return [path]
    if os.path.isdir(path):
        files = []
        for root, _, filenames in os.walk(path):
            for fn in sorted(filenames):
                if fn.endswith(".py") or fn.endswith(".pyx"):
                    files.append(os.path.join(root, fn))
        return files
    return []


def main():
    parser = argparse.ArgumentParser(
        description="Replace empty lines in function bodies with '#' comment lines."
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help="Files or directories to process",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check mode: exit 1 if any file needs changes",
    )
    parser.add_argument(
        "--diff",
        action="store_true",
        help="Show unified diff of proposed changes",
    )
    args = parser.parse_args()
    #
    changed_files = []
    for path in args.paths:
        for filepath in find_source_files(path):
            changed = process_file(filepath, check=args.check, diff=args.diff)
            if changed:
                changed_files.append(filepath)
    #
    if changed_files:
        action = "would modify" if args.check else "modified"
        for f in changed_files:
            print(f"{action}: {f}")
    else:
        print("No changes needed.")
    #
    if args.check and changed_files:
        sys.exit(1)


if __name__ == "__main__":
    main()
