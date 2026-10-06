#!/usr/bin/env python3

"""
Check that the math of a built Qlattice documentation tree renders in a
browser.\n
Sphinx only reports markup problems, so a LaTeX construct that MathJax rejects
is silently rendered as an error box in the generated html.  The common case
is a display environment (align, eqnarray, ...) nested in the "split" wrapper
that Sphinx adds around display math containing a line break, which MathJax
rejects with "Erroneous nesting of equation structures".\n
This script serves the built html, opens every page that contains math in a
headless browser, and reports every expression whose rendered output is a
MathJax error, together with its LaTeX source.\n
Chromium is required: it can dump the rendered DOM, which is how the error
boxes are detected.  Firefox cannot, so it can only be used for visual spot
checks.\n
Usage:
  ./nixpkgs/check-qlat-docs-math.py [html-dir] [options]\n
The default html-dir is tmp/result-docs/share/doc/qlat/html, i.e. the doc tree
produced by:
  nix-build nixpkgs/q-pkgs.nix -A pkgs.qlat_docs -o tmp/result-docs\n
Options:
  --browser NAME             Browser to use (default: chromium)
  --virtual-time-budget MS   Virtual time given to MathJax per page (default: 30000)
  --verbose                  List every page, not only the failing ones
  --help, -h                 Show this help message\n
Examples:
  ./nixpkgs/check-qlat-docs-math.py
  ./nixpkgs/check-qlat-docs-math.py tmp/result-docs/share/doc/qlat/html --verbose
"""

import argparse
import functools
import glob
import html
import http.server
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading

def show_help():
    print(__doc__.strip())
    sys.exit(0)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("html_dir", nargs="?", default="")
    parser.add_argument("--browser", default="chromium")
    parser.add_argument("--virtual-time-budget", type=int, default=30000)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--help", "-h", action="store_true")
    return parser.parse_args()

def get_project_root():
    script_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(script_path)

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

def resolve_html_dir(args):
    html_dir = args.html_dir
    if not html_dir:
        project_root = get_project_root()
        html_dir = os.path.join(
            project_root, "tmp", "result-docs", "share", "doc", "qlat", "html"
        )
    html_dir = os.path.abspath(html_dir)
    if not os.path.isfile(os.path.join(html_dir, "index.html")):
        print(f"Error: no built html found in: {html_dir}", file=sys.stderr)
        print("Build the documentation first, e.g.:", file=sys.stderr)
        print(
            "  nix-build nixpkgs/q-pkgs.nix -A pkgs.qlat_docs -o tmp/result-docs",
            file=sys.stderr,
        )
        sys.exit(1)
    return html_dir

def resolve_browser(name):
    if os.path.sep in name:
        if os.path.isfile(name) and os.access(name, os.X_OK):
            return name
    else:
        found = shutil.which(name)
        if found:
            return found
    print(f"Error: browser not found: {name}", file=sys.stderr)
    print(
        "Chromium is required: dumping the rendered DOM is how the MathJax "
        "error boxes are detected.",
        file=sys.stderr,
    )
    sys.exit(1)

def start_server(html_dir):
    handler = functools.partial(QuietHandler, directory=html_dir)
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"

def find_math_pages(html_dir):
    pages = []
    pattern = os.path.join(html_dir, "**", "*.html")
    for path in sorted(glob.glob(pattern, recursive=True)):
        with open(path, encoding="utf-8", errors="replace") as fp:
            if 'class="math' in fp.read():
                pages.append(os.path.relpath(path, html_dir))
    return pages

def clean_node(inner):
    inner = re.sub(r'<span class="eqno">.*?</span>', " ", inner, flags=re.S)
    text = html.unescape(re.sub(r"<[^>]+>", "", inner))
    for env in ("split", "align", "aligned"):
        text = text.replace(f"\\begin{{{env}}}", " ").replace(f"\\end{{{env}}}", " ")
    text = text.strip()
    for prefix, suffix in (("\\(", "\\)"), ("\\[", "\\]")):
        if text.startswith(prefix) and text.endswith(suffix):
            text = text[len(prefix) : len(text) - len(suffix)]
            break
    return re.sub(r"\s+", " ", text).strip()

def math_sources(path):
    with open(path, encoding="utf-8", errors="replace") as fp:
        text = fp.read()
    pattern = re.compile(r'<(span|div) class="[^"]*math[^"]*"[^>]*>(.*?)</\1>', re.S)
    return [
        re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", inner))).strip()
        for _, inner in pattern.findall(text)
    ]

def dump_dom(browser, url, virtual_time_budget):
    profile = tempfile.mkdtemp(prefix="qlat-docs-check-")
    try:
        cmd = [
            browser,
            "--headless",
            "--no-sandbox",
            "--disable-gpu",
            f"--user-data-dir={profile}",
            f"--virtual-time-budget={virtual_time_budget}",
            "--dump-dom",
            url,
        ]
        timeout = max(120, virtual_time_budget // 250)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    finally:
        shutil.rmtree(profile, ignore_errors=True)
    return result.stdout

def check_page(browser, base_url, html_dir, rel, virtual_time_budget):
    sources = math_sources(os.path.join(html_dir, rel))
    dom = dump_dom(browser, f"{base_url}/{rel}", virtual_time_budget)
    containers = re.findall(r"<mjx-container.*?</mjx-container>", dom, re.S)
    problems = []
    for index, source in enumerate(sources):
        if index >= len(containers):
            problems.append((index, "no rendered output", source))
            continue
        match = re.search(r'data-mjx-error="([^"]*)"', containers[index])
        if match:
            problems.append((index, html.unescape(match.group(1)), source))
    return sources, problems

def main():
    args = parse_args()
    #
    if args.help:
        show_help()
    #
    html_dir = resolve_html_dir(args)
    browser = resolve_browser(args.browser)
    #
    pages = find_math_pages(html_dir)
    if not pages:
        print(f"Error: no page with math found in: {html_dir}", file=sys.stderr)
        sys.exit(1)
    #
    print(f"Checking {len(pages)} page(s) with math in {html_dir}")
    print(f"Browser: {browser}")
    print()
    #
    server, base_url = start_server(html_dir)
    #
    n_math = 0
    n_problems = 0
    n_bad_pages = 0
    try:
        for rel in pages:
            sources, problems = check_page(
                browser, base_url, html_dir, rel, args.virtual_time_budget
            )
            n_math += len(sources)
            n_problems += len(problems)
            if problems:
                n_bad_pages += 1
                print(f"FAIL {rel}: {len(problems)} of {len(sources)} expression(s)")
                for index, message, source in problems:
                    print(f"     [{index}] {message}")
                    print(f"         {source[:150]}")
            elif args.verbose:
                print(f"OK   {rel}: {len(sources)} expression(s)")
    finally:
        server.shutdown()
        server.server_close()
    #
    print()
    print(
        f"Summary: {n_math} expression(s) checked, "
        f"{n_problems} error(s) in {n_bad_pages} of {len(pages)} page(s)"
    )
    if n_problems:
        sys.exit(1)

main()
