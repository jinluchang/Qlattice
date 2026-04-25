#!/usr/bin/env python3

"""
Lint tool: replace empty lines inside function bodies with indented comment lines.\n
Supported file types:
    .py          Python      — formatted with ruff, AST-based function detection
    .pyx         Cython      — regex-based function detection
    .c .cc .cpp .cxx .h .hpp  C/C++ — formatted with clang-format, brace-based function detection\n
Comment markers:
    Python/Cython: #
    C/C++:         //\n
Design goals:
    - Remove all empty lines within function/method bodies to enforce compact code.
    - Use indented comment lines as visual separators when logical grouping is needed.
    - Preserve the indentation of the surrounding context for each replaced empty line.
    - Leave empty lines outside function bodies untouched (module level, between functions).
    - Remove empty lines inside multi-line string literals (e.g., docstrings) while
      preserving content by appending "\\n" to the preceding line.
    - Handle nested functions correctly by scoping replacements to each function's body.\n
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
    """Return 1-indexed line numbers inside function/method bodies for Python files.\n
    Uses the ast module to walk the parse tree and collect line ranges
    of FunctionDef and AsyncFunctionDef nodes (including nested ones).
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
    """Return 1-indexed line numbers inside function/method bodies for Cython files.\n
    Uses regex-based detection (def/cdef/cpdef) since ast.parse cannot handle .pyx syntax.
    Scans for lines indented deeper than the function signature as the body extent.
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


CPP_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".h", ".hpp"}


def get_function_body_lines_cpp(source_lines):
    """Return 1-indexed line numbers inside function/method bodies for C/C++ files.\n
    Uses brace-depth tracking to find function bodies. Detects functions by looking
    for a ')' followed by optional qualifiers and then '{', while excluding control
    structures (if, for, while, switch, catch, else). Skips preprocessor directives
    and comments.
    """
    body_lines = set()
    i = 0
    while i < len(source_lines):
        line = source_lines[i]
        stripped = line.strip()
        #
        # Skip preprocessor directives
        if stripped.startswith("#"):
            i += 1
            continue
        #
        # Skip single-line comments
        if stripped.startswith("//"):
            i += 1
            continue
        #
        # Skip multi-line comments (advance i past them)
        if "/*" in stripped:
            while i < len(source_lines) and "*/" not in source_lines[i]:
                i += 1
            i += 1
            continue
        #
        if "{" in stripped:
            # Collect lines from '{' backwards to find function signature
            sig_lines = []
            found_paren = False
            for k in range(i, max(-1, i - 15), -1):
                sig_lines.insert(0, source_lines[k].rstrip())
                if "(" in source_lines[k]:
                    found_paren = True
                    break
            if not found_paren:
                i += 1
                continue
            sig = " ".join(sig_lines)
            #
            # Remove comments and strings from signature for analysis
            clean = re.sub(r"//.*", "", sig)
            clean = re.sub(r"/\*.*?\*/", "", clean)
            clean = re.sub(r'"(?:[^"\\]|\\.)*"', '""', clean)
            clean = re.sub(r"'(?:[^'\\]|\\.)*'", "''", clean)
            #
            # Detect '{' position in cleaned signature
            brace_pos = clean.find("{")
            before_brace = clean[:brace_pos].rstrip()
            #
            # Check for function call pattern: '...' before '{' ends with ')'
            # (after removing trailing qualifiers like 'const', 'override', 'noexcept')
            qualifiers = r"(?:\b(?:const|override|noexcept|final|volatile|&|&&|\=\s*0|=\s*default|=\s*delete)\s*)*"
            func_pat = re.compile(r"\)\s*" + qualifiers + r"$")
            is_func = bool(func_pat.search(before_brace))
            #
            # Exclude control structures: if, for, while, switch, catch, else
            if is_func:
                control_pat = re.compile(r"\b(?:if|for|while|switch|catch|else)\s*\(")
                if control_pat.search(before_brace):
                    is_func = False
            #
            if not is_func:
                i += 1
                continue
            #
            # Find opening brace position in original source
            brace_line = i
            brace_col = source_lines[i].find("{")
            depth = 1
            j = i
            col = brace_col + 1
            while j < len(source_lines) and depth > 0:
                line_text = source_lines[j]
                while col < len(line_text):
                    ch = line_text[col]
                    if ch == "{":
                        depth += 1
                    elif ch == "}":
                        depth -= 1
                        if depth == 0:
                            break
                    elif ch == '"':
                        col += 1
                        while col < len(line_text) and line_text[col] != '"':
                            if line_text[col] == "\\":
                                col += 1
                            col += 1
                    elif ch == "'":
                        col += 1
                        while col < len(line_text) and line_text[col] != "'":
                            if line_text[col] == "\\":
                                col += 1
                            col += 1
                    elif ch == "/" and col + 1 < len(line_text):
                        if line_text[col + 1] == "/":
                            break
                        elif line_text[col + 1] == "*":
                            col += 2
                            while col < len(line_text):
                                if (
                                    line_text[col] == "*"
                                    and col + 1 < len(line_text)
                                    and line_text[col + 1] == "/"
                                ):
                                    col += 2
                                    break
                                col += 1
                            continue
                    col += 1
                if depth == 0:
                    break
                j += 1
                col = 0
            #
            if depth == 0:
                for ln in range(brace_line + 1, j + 2):
                    body_lines.add(ln)
            i = j + 1
            continue
        #
        i += 1
    return body_lines


def get_indent(line):
    """Return the leading whitespace string of a line."""
    return line[: len(line) - len(line.lstrip())]


def get_string_lines(source):
    """Return 1-indexed line numbers inside multi-line string literals.\n
    Uses the tokenize module to find STRING tokens that span multiple lines.
    Only meaningful for Python/Cython sources.
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
    """Collapse runs of identical bare comment lines ('#' or '//') into a single line."""
    result = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped in ("#", "//"):
            j = i + 1
            while j < len(lines) and lines[j].strip() == stripped:
                j += 1
            if j > i + 1:
                result.append(line)
                i = j
                continue
        result.append(line)
        i += 1
    return result


def transform_file(source, is_cython=False, is_cpp=False):
    """Replace empty lines inside function bodies with indented comment lines.\n
    The comment marker is '#' for Python/Cython and '//' for C/C++.
    Exactly one of is_cython/is_cpp should be True for non-Python files;
    when both are False, Python AST-based detection is used.
    """
    lines = source.splitlines(keepends=True)
    if is_cpp:
        body_lines = get_function_body_lines_cpp(lines)
        string_lines = set()
        comment_marker = "//"
    elif is_cython:
        body_lines = get_function_body_lines_cython(lines)
        string_lines = set()
        comment_marker = "#"
    else:
        body_lines = get_function_body_lines(lines)
        comment_marker = "#"
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
                result.append(f"{indent}{comment_marker}\n")
            else:
                result.append(f"{indent}{comment_marker}")
            i += 1
        else:
            result.append(line)
            i += 1
    #
    result = collapse_consecutive_comments(result)
    return "".join(result)


def run_ruff(source):
    """Run ruff format and ruff check --fix on Python source code."""
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


def run_clang_format(source, filepath):
    """Run clang-format on C/C++ source code.\n
    Uses --assume-filename so clang-format picks up the correct language and
    .clang-format config from the file's original directory.
    """
    suffix = os.path.splitext(filepath)[1]
    with tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False) as f:
        f.write(source)
        tmp = f.name
    try:
        subprocess.run(
            ["clang-format", "-i", f"--assume-filename={filepath}", tmp],
            capture_output=True,
        )
        with open(tmp, "r") as f:
            return f.read()
    finally:
        os.unlink(tmp)


def process_file(filepath, check=False, diff=False):
    """Process a single source file through its language-specific formatter and transform.\n
    Pipeline by extension:
        .py   — ruff format + ruff check --fix, then transform
        .pyx  — transform only (no formatter)
        C/C++ — clang-format, then transform\n
    Returns True if the file was (or would be) modified.
    """
    ext = os.path.splitext(filepath)[1]
    with open(filepath, "r") as f:
        original = f.read()
    #
    if ext in CPP_EXTENSIONS:
        formatted = run_clang_format(original, filepath)
        transformed = transform_file(formatted, is_cpp=True)
    elif ext == ".pyx":
        transformed = transform_file(original, is_cython=True)
    elif ext == ".py":
        ruffed = run_ruff(original)
        transformed = transform_file(ruffed)
    else:
        return False
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
    """Find all supported source files (.py, .pyx, .c, .cc, .cpp, .cxx, .h, .hpp) at path.\n
    If path is a single file, returns it (if supported). If a directory, walks
    recursively and returns all matching files.
    """
    all_extensions = {".py", ".pyx"} | CPP_EXTENSIONS
    if os.path.isfile(path) and os.path.splitext(path)[1] in all_extensions:
        return [path]
    if os.path.isdir(path):
        files = []
        for root, _, filenames in os.walk(path):
            for fn in sorted(filenames):
                if os.path.splitext(fn)[1] in all_extensions:
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
