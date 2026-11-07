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

   A numbered display must not contain a LaTeX line break (`\\`): Sphinx wraps
   such a display in `\begin{split}`, where MathJax rejects `\tag`, and the
   page renders an error box. Write the relation on one line, split
   independent relations into separate displays, or use `\cr` inside a
   `\begin{array}` for a piecewise definition. A numbered display never uses
   `\begin{aligned}`, `\begin{align}` or `\nonumber` — numbering is
   unconditional.

4. **Numbers only from the script.** Run
   `python3 nixpkgs/number-docs-equations.py <file.md>`; never hand-edit
   `\tag`; run `--check` before hand-off.

5. **Reference a display on this page as `Eq. (N)`.** Run the script first to
   number a new display, then cite the number it printed; never hand-pick one.
   The script also rewrites existing references when the numbering shifts.
   Citations of another paper's equation are not this form — see
   `## Carve-outs`.

6. **One relation per display.** A chained equality (`A = B = C`), an
   initial-value pair, or a relation with its domain of validity is one
   relation: it stays one display, written on a single line and without a
   LaTeX line break. Genuinely independent relations get separate displays —
   different left-hand sides, a `\gets` update followed by an equality, a
   distribution followed by its Boltzmann weight.

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

- Markdown only (`.md`): run it on Markdown files alone. It has no extension
  guard, so pointing it at a `.rst` or a source file would rewrite that file
  the same way.
- All-or-nothing per file: the new text is built in full and written in one
  pass, only when it differs.
- A reference with no numbered equation is reported as `DANGLING` with its line
  number; the tags are still written so the prose can be fixed.
- `--check` validates instead of rewriting. Exit `0` is clean; exit `1` means
  stale tags, a stale reference, a dangling reference, a render-unsafe display,
  or a missing/unreadable file; exit `2` is command-line misuse.
- An unclosed `$$` and a `$$` fence sharing a line with text are warned about on
  stderr; the warning does not by itself fail the run.

This tool only reads the Markdown source; it cannot tell whether MathJax can
render the result. That is `./nixpkgs/check-qlat-docs-math.py`, which builds the
docs, opens every page containing math in headless chromium and reports the
expressions MathJax could not render. Sphinx does not report MathJax failures,
so both checks matter: one keeps the `\tag` numbers honest, the other keeps the
page from showing error boxes.

Numbering is unconditional: every display gets a `\tag{N}`, even on a page
whose prose never writes `Eq. (N)`, so a reference can be added later without
renumbering the page.

## Scope

- Applies to every page under `docs/` that contains math, including the module
  reference pages and the `docs/contents/` guides. Nothing is grandfathered:
  every display on every such page is numbered.
- The four pages that used chained `\begin{aligned}` displays
  (`docs/contents/contraction.md`, `docs/contents/hmc.md`,
  `docs/contents/shmc.md` and `docs/contents/jackknife-random.md`) are being
  retrofitted under these rules — one relation per display, no line break
  inside a numbered display, the `## Carve-outs` respected — so a chained or
  multi-relation display is not allowed anywhere.
- The tool numbers the displays in file order; never hand-edit a `\tag`, and
  re-run the tool after inserting or deleting a display so the tags and the
  `Eq. (N)` references stay in sync.
- Table cells, size labels, code and external-literature references cannot
  become displays: they are exempt, and `## Carve-outs` says exactly what to
  do there.

## Carve-outs

Places that cannot follow the display rule. Leave them as they are unless they
are plainly wrong.

- **Markdown table cells.** A `$$` display cannot exist inside a table cell,
  so cell math stays concise inline LaTeX (or a code span) and is exempt from
  the "everything else is a display" rule.
- **Dimension and size labels.** `12×12`, `3×3`, `1 x 1 x 1 x 1` and
  `288×288×288×576` are labels, not formulas: leave them as plain text or a
  simple inline/code token, and never promote them to a display.
- **Fenced code blocks and inline code.** They stay exempt from every rule
  here; they hold commands, captured output and literal syntax examples.
- **External-literature citations.** For another paper's equation write
  `equation (4)`, `Eq. 4`, or a code span — never `Eq. (4)` with digits: the
  numbering tool reads that as a reference to a display on this page and
  reports it as `DANGLING`.

## Verification

Check each page you touch before hand-off, and build and render the site after
any documentation change that touches equations:

```bash
python3 nixpkgs/number-docs-equations.py --check <page>
nix-build nixpkgs/q-pkgs.nix -A pkgs.qlat_docs -o tmp/result-docs
./nixpkgs/check-qlat-docs-math.py
```

`<page>` is the Markdown file that changed, e.g.
`docs/contents/coding-doc-equations.md`; this page itself has no display
equations, so `--check` reports zero. `--check` also rejects a render-unsafe
display: one containing a LaTeX line break (`\\`) or an `aligned`/`align`
environment, which MathJax cannot tag. That keeps the source honest, but it
cannot prove the rendered output — the chromium math check above is the final
gate. `AGENTS.md` `## Documentation` is the authority for the
build-and-render step; follow it together with the `--check` gate above.
