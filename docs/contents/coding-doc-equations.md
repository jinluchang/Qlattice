# Coding doc equations

The rules for equations in the Markdown prose pages under `docs/`: every
formula is LaTeX, inline math is one token, and everything larger is a numbered
display. Both the tags and the `Eq. (N)` references come from
`nixpkgs/number-docs-equations.py`, never from a hand edit.

## Core rules

1. **LaTeX only.** Every formula is LaTeX, never plain text.

2. **Inline math is one token.** Use `$…$` only for a symbol (`$n$`,
   `$\alpha$`), a subscripted or superscripted symbol (`$f_n$`, `$g_{n+1}$`,
   `$N^2$`), a plain number (`$0.1$`), or a ratio of two such tokens (`$t/6$`,
   `$1/6$`). Nothing else.

3. **Everything else is a numbered display.** Any relation or operator (`=`,
   `<`, `>`, `+`, `\times`, `\equiv`, `\sim`, `\to`, `\approx`), any `\sum`,
   `\prod`, `\int`, a multi-term `\frac` or `\binom`, and any multi-token
   expression goes on its own lines inside `$$ … $$`. The fences sit on their
   own lines too — never inline with text.

4. **Numbers only from the script.** Run
   `python3 nixpkgs/number-docs-equations.py <file.md>`; never hand-edit
   `\tag`; run `--check` before hand-off.

5. **Reference as `Eq. (N)`.** Run the script first to number a new display,
   then cite the number it printed; never hand-pick one. The script also
   rewrites existing references when the numbering shifts.

6. **One relation per display.** Split chained derivations into separate
   numbered lines. A chained equality, an initial-value pair, or a relation
   with its domain of validity counts as one relation.

7. **Fenced code blocks are exempt.** They hold commands and captured output
   verbatim.

8. **Rephrase, do not smuggle.** State the value in prose: with `$\alpha$`
   equal to `$0.1$`, never with `$\alpha = 0.1$`.

## Allowed and disallowed inline math

| allowed inline | not allowed inline |
|---|---|
| `$n$`, `$f_n$`, `$N^2$`, `$0.1$`, `$t/6$`, `$1/6$` | `$t = 0.1$`, `$\alpha = 0.1$`, `$\|t\| < 1/2$` |
| `$\alpha$`, `$g_{n+1}$`, `$12$`, `$\Phi$` | `$(1-2t)^{-1/2}$`, `$t^n/n!$`, `$\sum_n f_n$` |

The right column only becomes a display. Every entry is a code span so this
page shows the syntax instead of rendering it — and so the numbering tool,
which skips inline code, leaves it alone.

## Worked example

Write the display, then let the script number it:

```
Eq. (4) is
$$
S(t) = (1-2t)^{-1/2}
$$

after `python3 nixpkgs/number-docs-equations.py notes/file.md` the same block is

$$
S(t) = (1-2t)^{-1/2}
\tag{4}
$$
```

## The numbering tool

`nixpkgs/number-docs-equations.py` numbers every `$$` display in document
order: it strips any pre-existing `\tag{...}` and writes `\tag{N}` as the last
body line before the closing fence. A single-line `$$ … $$` is expanded to the
canonical four-line form (open fence, body, `\tag{N}`, close fence). The old tag
values build an old→new map, so the command is idempotent and renumbering
rewrites the prose references `Eq. (N)`, `Eqs. (N)`, `Eqs. (N) and (M)`,
`Eqs. (N)-(M)` and `Eq. (N)-(M)` outside fenced code blocks and inline code
spans.

- Markdown only (`.md`); reStructuredText and source files are never touched.
- All-or-nothing per file: the new text is built in full and written in one
  pass, only when it differs.
- A reference with no numbered equation is reported as `DANGLING` with its line
  number; the tags are still written so the prose can be fixed.
- `--check` validates instead of rewriting. Exit `0` is clean; exit `1` means
  stale tags, a stale reference, a dangling reference, or a missing/unreadable
  file; exit `2` is command-line misuse.
- An unclosed `$$` and a `$$` fence sharing a line with text are warned about on
  stderr; the warning does not by itself fail the run.

This tool only reads the Markdown source; it cannot tell whether MathJax can
render the result. That is `./nixpkgs/check-qlat-docs-math.py`, which builds the
docs, opens every page containing math in headless chromium and reports the
expressions MathJax could not render. Sphinx does not report MathJax failures,
so both checks matter: one keeps the `\tag` numbers honest, the other keeps the
page from showing error boxes.

Numbering is only needed on a page that cites its displays: if no `Eq. (N)`
reference is used, do not tag the page and do not run the tool on it.

## Scope

- Applies to new or edited prose pages that reference equations.
- The existing `\begin{aligned}` groups in `docs/contents/contraction.md`,
  `docs/contents/hmc.md`, `docs/contents/shmc.md` and
  `docs/contents/jackknife-random.md` are grandfathered: 54 aligned blocks, no
  tags, no references. Do not retrofit tags into them.
- Never run the tool on a page whose tagging is not intended: it inserts a
  `\tag{N}` into every display and rewrites every reference it can match.

## Verification

Check one page before hand-off, and build and render the site after any
documentation change that touches equations:

```bash
python3 nixpkgs/number-docs-equations.py --check <page>
nix-build nixpkgs/q-pkgs.nix -A pkgs.qlat_docs -o tmp/result-docs
./nixpkgs/check-qlat-docs-math.py
```

`<page>` is the Markdown file that changed, e.g.
`docs/contents/coding-doc-equations.md`; this page itself has no display
equations, so `--check` reports zero. `AGENTS.md` `## Documentation` is the
authority for the build-and-render step; follow it together with the `--check`
gate above.
