#!/usr/bin/env python3
"""Number the display equations of one or more Markdown files and keep `Eq. (N)` references in sync.\n
Every display-math block (`$$ ... $$`) receives a `\\tag{N}`, N counted in document order. Display
blocks must open and close `$$` on their own line; a single-line `$$ ... $$` block is converted to
that canonical form. Any pre-existing `\\tag{...}` is stripped first, and the old tag value is used to
map old numbers to new ones, so the command is IDEMPOTENT and it RENUMBERS the prose references.\n
References rewritten (outside fenced code and outside inline code spans):
    Eq. (N)   Eqs. (N)   Eqs. (N) and (M)   Eqs. (N)-(M)   Eq. (N)-(M)
A reference whose number has no mapping (the equation was deleted, or was never numbered) is
reported as DANGLING with its line number and makes the command exit non-zero — the tags are still
written, so the author can fix the prose.\n
Usage:
    python3 nixpkgs/number-docs-equations.py notes/foo.md [notes/bar.md ...]
    python3 nixpkgs/number-docs-equations.py --check notes/foo.md     # exit 1 on stale tags, stale refs, or dangling
"""

from __future__ import annotations

import argparse
import os
import re
import sys

TAG_RE = re.compile(r"[ \t]*\\tag\{([^}]*)\}")
FENCE_RE = re.compile(r"^[ \t]*(```|~~~)")
CODE_SPAN_RE = re.compile(r"`[^`]*`")
REF_RE = re.compile(r"Eqs?\.\s*(?:\(\d+\)(?:\s*(?:,|and|–|—|-|to)\s*)?)+")
NUM_RE = re.compile(r"\((\d+)\)")
RANGE_SEP = re.compile(r"[ \t]*(?:–|—|-|to)[ \t]*")

def code_fence_mask(lines: list[str]) -> list[bool]:
    """True for the fence lines themselves and for every line inside a fenced block."""
    inside = False
    mask = []
    for ln in lines:
        if FENCE_RE.match(ln):
            inside = not inside
            mask.append(True)
            continue
        mask.append(inside)
    return mask

def scan_blocks(lines: list[str], code: list[bool], warn) -> list[dict]:
    """Locate display blocks; warn about an unclosed `$$` instead of skipping it silently."""
    blocks: list[dict] = []
    i = 0
    while i < len(lines):
        if code[i]:
            i += 1
            continue
        s = lines[i].strip()
        if s.startswith("$$"):
            if "$$" in s[2:]:
                blocks.append({"start": i, "end": i, "single": True})
                i += 1
                continue
            j = i + 1
            while j < len(lines) and not lines[j].strip().endswith("$$"):
                j += 1
            if j >= len(lines):
                warn(
                    f"unclosed display block: `$$` at line {i + 1} has no closing `$$`"
                )
                i += 1
                continue
            blocks.append({"start": i, "end": j, "single": False})
            i = j + 1
            continue
        if "$$" in CODE_SPAN_RE.sub("", s):
            warn(
                f"stray `$$` at line {i + 1}: a display fence must sit on its own line"
            )
        i += 1
    return blocks

def old_tag(lines: list[str], block: dict) -> int | None:
    for k in range(block["start"], block["end"] + 1):
        m = TAG_RE.search(lines[k])
        if m and m.group(1).strip().isdigit():
            return int(m.group(1))
    return None

def block_lines(lines: list[str], block: dict, number: int) -> list[str]:
    """Canonical form: opening `$$`, stripped body, `\\tag{N}`, closing `$$`."""
    first, last = block["start"], block["end"]
    indent = lines[first][: len(lines[first]) - len(lines[first].lstrip())]
    tag = f"{indent}\\tag{{{number}}}"
    if block["single"]:
        body = TAG_RE.sub("", lines[first]).strip()
        inner = body[2:-2].strip() if body.endswith("$$") else body[2:].strip()
        return [f"{indent}$$", f"{indent}{inner}", tag, f"{indent}$$"]
    opening = TAG_RE.sub("", lines[first]).rstrip()
    body = []
    for k in range(first + 1, last):
        piece = TAG_RE.sub("", lines[k]).rstrip()
        if piece.strip():
            body.append(piece)
    closing = TAG_RE.sub("", lines[last]).rstrip()
    if closing.strip() != "$$":
        extra = closing[:-2].rstrip()
        if extra.strip():
            body.append(extra)
        closing = f"{indent}$$"
    return [opening, *body, tag, closing]

def rewrite_refs(line: str, mapping: dict[int, int], lineno: int, stats: dict) -> str:
    """Rewrite every `Eq.(s) (N)...` reference in a line that is neither fenced nor inline code."""
    #
    def process(segment: str) -> str:
        def repl(match: re.Match) -> str:
            text = match.group(0)
            nums = [int(x) for x in NUM_RE.findall(text)]
            new_nums = []
            for n in nums:
                if n in mapping:
                    new_nums.append(mapping[n])
                else:
                    stats["dangling"].append((lineno, n))
                    new_nums.append(n)
            if new_nums != nums:
                stats["changed"].append((lineno, tuple(nums), tuple(new_nums)))
                if len(nums) > 1 and RANGE_SEP.search(text):
                    stats["ranges"].append((lineno, tuple(nums), tuple(new_nums)))
            it = iter(new_nums)
            return re.sub(r"\(\d+\)", lambda _m: f"({next(it)})", text)
        #
        return REF_RE.sub(repl, segment)
    #
    out, pos = [], 0
    for m in CODE_SPAN_RE.finditer(line):
        out.append(process(line[pos : m.start()]))
        out.append(m.group(0))  # inline code is never a reference
        pos = m.end()
    out.append(process(line[pos:]))
    return "".join(out)

def number_file(path: str, warn=print) -> dict:
    with open(path, encoding="utf-8") as f:
        original = f.read()
    trailing_newline = original.endswith("\n")
    lines = original.split("\n")
    if trailing_newline:
        lines = lines[:-1]
    code = code_fence_mask(lines)
    blocks = scan_blocks(lines, code, warn)
    #
    mapping: dict[int, int] = {}
    for k, block in enumerate(blocks, start=1):
        t = old_tag(lines, block)
        if t is not None:
            mapping[t] = k
    #
    stats = {
        "path": path,
        "blocks": len(blocks),
        "changed": [],
        "ranges": [],
        "dangling": [],
    }
    out: list[str] = []
    inset: set[int] = set()
    cursor = 0
    for k, block in enumerate(blocks, start=1):
        out.extend(lines[cursor : block["start"]])
        new_lines = block_lines(lines, block, k)
        inset.update(range(len(out), len(out) + len(new_lines)))
        out.extend(new_lines)
        cursor = block["end"] + 1
    out.extend(lines[cursor:])
    #
    code_out = code_fence_mask(out)
    final: list[str] = []
    for idx, ln in enumerate(out):
        if idx in inset or code_out[idx]:
            final.append(ln)
            continue
        final.append(rewrite_refs(ln, mapping, idx + 1, stats))
    text = "\n".join(final) + ("\n" if trailing_newline else "")
    stats["text"] = text
    stats["old"] = original
    stats["mapping"] = mapping
    return stats

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="number-docs-equations.py",
        description=__doc__,
        epilog="Exit codes: 0 success; 1 stale tags / stale references / dangling reference / "
        "missing or unreadable file; 2 command-line misuse.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "files",
        nargs="+",
        metavar="FILE",
        help="Markdown file(s) to number and to check references in",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate only: exit 1 when tags or references are stale, or a "
        "reference dangles",
    )
    return parser

def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    files, check = args.files, args.check
    #
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        for f in missing:
            print(f"error: no such file: {f}", file=sys.stderr)
        return 1
    #
    warnings: list[str] = []
    stale_tags: list[str] = []
    stale_refs: list[dict] = []
    dangling: list[tuple[str, int, int]] = []
    total = 0
    for path in files:
        st = number_file(path, warn=lambda msg, p=path: warnings.append(f"{p}: {msg}"))
        total += st["blocks"]
        old_wo_tags = TAG_RE.sub("", st["old"])
        new_wo_tags = TAG_RE.sub("", st["text"])
        tags_changed = old_wo_tags != new_wo_tags
        refs_changed = [c for c in st["changed"]]
        dangling += [(path, ln, n) for ln, n in st["dangling"]]
        if check:
            if tags_changed:
                stale_tags.append(path)
            if refs_changed:
                stale_refs.append({"path": path, "changes": refs_changed})
            continue
        if st["text"] != st["old"]:
            with open(path, "w", encoding="utf-8") as f:
                f.write(st["text"])
        print(f"{path}: {st['blocks']} display equation(s) numbered 1..{st['blocks']}")
        if st["mapping"]:
            shifted = {o: n for o, n in st["mapping"].items() if o != n}
            if shifted:
                sample = ", ".join(f"{o}->{n}" for o, n in sorted(shifted.items())[:6])
                print(f"  renumbered {len(shifted)} existing equation(s): {sample}")
        if refs_changed:
            for ln, old, new in refs_changed:
                print(f"  line {ln}: Eq(s). {old} -> {new}")
        for ln, old, new in st["ranges"]:
            contiguous = (
                len(old) == 2 and len(new) == 2 and new[1] - new[0] == old[1] - old[0]
            )
            print(
                f"  line {ln}: range {old} -> {new}"
                + ("" if contiguous else "  (NOT contiguous: check the range by hand)")
            )
    for w in warnings:
        print(f"warning: {w}", file=sys.stderr)
    for path, ln, n in dangling:
        print(
            f"DANGLING reference: {path} line {ln}: Eq. ({n}) has no numbered equation",
            file=sys.stderr,
        )
    #
    if check:
        reasons = []
        if stale_tags:
            reasons.append("stale tag numbering in " + ", ".join(stale_tags))
        if stale_refs:
            reasons.append(
                "stale references in " + ", ".join(d["path"] for d in stale_refs)
            )
        if dangling:
            reasons.append(f"{len(dangling)} dangling reference(s)")
        for d in stale_refs:
            for ln, old, new in d["changes"]:
                reasons.append(f"  {d['path']} line {ln}: Eq(s). {old} -> {new}")
        if reasons:
            for r in reasons:
                print(r, file=sys.stderr)
            return 1
        print(f"{len(files)} file(s) checked, {total} display equation(s)")
        return 0
    if dangling:
        return 1
    print(f"{len(files)} file(s), {total} display equation(s)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
