# Plan: Fix TODOs in nixpkgs/

There are 7 TODOs across two files: `options.nix` (4) and `q-pkgs.nix` (2),
plus 1 new TODO (define `import-nixgl`).

---

## 1. Move `nixpkgs-release` to `options.nix`

**File:** `q-pkgs.nix:22`

```
# TODO: move to options.nix
nixpkgs-release = lib.trivial.release;
```

### Procedure

**Step 1.1 — Add definition to `options.nix`**

In `options.nix`, insert after line 27 (`lib = o-pkgs.lib;`):

```nix
  nixpkgs-release = lib.trivial.release;
```

**Step 1.2 — Export from `options.nix`**

In `options.nix` return set (line 188-195), add `nixpkgs-release` to the
existing `inherit` lines. For example, change line 191 from:

```nix
  inherit nixgl-src is-linux lib;
```

to:

```nix
  inherit nixgl-src is-linux lib nixpkgs-release;
```

**Step 1.3 — Update `q-pkgs.nix` to inherit instead of define**

In `q-pkgs.nix`, replace lines 22-23:

```nix
  # TODO: move to options.nix
  nixpkgs-release = lib.trivial.release;
```

with nothing (delete them). Then add `nixpkgs-release` to the `inherit` block
at line 17. Change:

```nix
  inherit (opts) nixgl-src is-linux lib o-pkgs;
```

to:

```nix
  inherit (opts) nixgl-src is-linux lib o-pkgs nixpkgs-release;
```

**Step 1.4 — Verify `nixpkgs-release` is still used in `q-pkgs.nix`**

`nixpkgs-release` appears in the `builtins.deepSeq` list at line 732. No
change needed there — it will resolve via `inherit` instead of local def.

---

## 2. Move `options-list` section to `options.nix`

**File:** `q-pkgs.nix:657`

```
# TODO: move this section to options.nix
```

The section spans lines 659-698 (`options-list`, `qlat-name-list`,
`qlat-name-list-file-from-str`, `qlat-name-list-file`).

### Procedure

**Step 2.1 — Add `options-list` to `options.nix`**

In `options.nix`, insert after the `mk-qlat-name` definition (after line 186,
before the `in` block). Copy the entire block from `q-pkgs.nix` lines 659-686:

```nix
  options-list = [
    {}
    { use-cps = false; use-grid-gpt = false; }
    { use-cps = false; }
    { use-grid-gpt = false; }
    { use-cuda-software = true; }
    { use-cuda = true; }
    { use-cudasupport = true; }
    { use-ucx = false; }
    { use-clang = true; }
    { use-pypi = version-pypi; }
    #
    { use-clang = true; use-ucx = false; }
    { use-cuda = true; use-ucx = false; }
    { use-grid-gpt = false; use-cubaquad = false; }
    { use-grid-gpt = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-ucx = false; }
    { use-cps = false; use-grid-gpt = false; use-clang = true; use-ucx = false; }
    #
    { use-cps = false; use-ucx = false; }
    { use-cps = false; use-clang = true; use-ucx = false; }
    { use-cps = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-cuda-software = true; }
    { use-cps = false; use-grid-gpt = false; use-cuda = true; }
    { use-cps = false; use-grid-gpt = false; use-cudasupport = true; }
  ];
```

**Step 2.2 — Add `version-pypi` to `options.nix`**

`options-list` references `version-pypi` at line 669. This value is defined
in `q-pkgs.nix` line 20 as `"1.5"`. Add it to `options.nix` (e.g., after
`mk-qlat-name`):

```nix
  version-pypi = "1.5";
```

**Step 2.3 — Add `qlat-name-list` to `options.nix`**

Copy from `q-pkgs.nix` line 688:

```nix
  qlat-name-list = let r = lib.lists.unique (builtins.map mk-qlat-name (builtins.map mk-options options-list)); in builtins.deepSeq r r;
```

**Step 2.4 — Add `qlat-name-list-file-from-str` to `options.nix`**

Copy from `q-pkgs.nix` lines 690-691:

```nix
  qlat-name-list-file-from-str = builtins.toFile "qlat-name-list"
  (builtins.foldl' (s: v: s + "${v}\n") "" qlat-name-list);
```

**Step 2.5 — Add `qlat-name-list-file` to `options.nix`**

Copy from `q-pkgs.nix` lines 693-698. Replace `o-pkgs.runCommand` with
`o-pkgs.runCommand` (or `pkgs.runCommand` if TODO 4 is already done):

```nix
  qlat-name-list-file = o-pkgs.runCommand
  "qlat-name-list"
  {}
  ''
    cp -v ${qlat-name-list-file-from-str} $out
  '';
```

**Step 2.6 — Export all four from `options.nix`**

Add to the return set (line 188-195):

```nix
  inherit options-list qlat-name-list qlat-name-list-file-from-str;
  inherit qlat-name-list-file version-pypi;
```

**Step 2.7 — Remove definitions from `q-pkgs.nix`**

Delete `q-pkgs.nix` lines 655-698 (the `# TODO` comment, `options-list`,
`qlat-name-list`, `qlat-name-list-file-from-str`, `qlat-name-list-file`).

Also remove `version-pypi = "1.5";` from line 20 if it is no longer used
elsewhere in `q-pkgs.nix`. (It is referenced at line 669 inside
`options-list` which is now in `options.nix`. Verify it is not used in
`mk-q-pkgs` or other functions in `q-pkgs.nix`.)

**Step 2.8 — Inherit from `opts` in `q-pkgs.nix`**

Add to the `inherit` block at line 18:

```nix
  inherit (opts) options-list qlat-name-list qlat-name-list-file-from-str;
  inherit (opts) qlat-name-list-file version-pypi;
```

**Step 2.9 — Update `o-pkgs` usage in `q-pkgs.nix`**

After this move, `o-pkgs` is no longer used in `q-pkgs.nix` (its only usage
was `o-pkgs.runCommand` at line 693 which is now in `options.nix`). Verify
all `o-pkgs` references in `q-pkgs.nix` are gone before proceeding to TODO 4.

---

## 3. Copy `force` to `q-pkgs.nix`

**File:** `options.nix:21`

```
# TODO: copy to q-pkgs.nix (and use it there when needed)
force = x: builtins.deepSeq x x;
```

### Procedure

**Step 3.1 — Add `force` to `q-pkgs.nix`**

In `q-pkgs.nix`, insert after line 18 (`inherit (opts) options-default ...`):

```nix
  force = x: builtins.deepSeq x x;
```

**Step 3.2 — Keep `force` in `options.nix`**

Do NOT remove `force` from `options.nix` — it is used at line 133:
`in force opts;`. Leave it as-is. The TODO only says "copy", not "move".

**Step 3.3 — Remove the TODO comment from `options.nix`**

Change lines 21-22 from:

```nix
  # TODO: copy to q-pkgs.nix (and use it there when needed)
  force = x: builtins.deepSeq x x;
```

to:

```nix
  force = x: builtins.deepSeq x x;
```

---

## 4. Rename `o-pkgs` to `pkgs` in `options.nix`

**File:** `options.nix:24`

```
# TODO: rename to pkgs (and do not return it)
o-pkgs = import-nixpkgs-wd {};
```

### Procedure

**Step 4.1 — Rename in `options.nix`**

Change line 25 from:

```nix
  o-pkgs = import-nixpkgs-wd {};
```

to:

```nix
  pkgs = import-nixpkgs-wd {};
```

Update all references inside `options.nix`:

| Line | Before | After |
|------|--------|-------|
| 24 | `# TODO: rename to pkgs (and do not return it)` | *(delete line)* |
| 27 | `lib = o-pkgs.lib;` | `lib = pkgs.lib;` |

Note: line 58 has `pkgs = pkgs;` inside the `options-default` `let` block.
This is a different scope — the inner `pkgs` shadows the outer one. This is
already correct and does not conflict because `o-pkgs` was only used at lines
25, 27, and in the return set.

**Step 4.2 — Remove `o-pkgs` from `options.nix` return set**

Change line 192 from:

```nix
  inherit o-pkgs;
```

to (delete the line entirely, or if other things are on it, just remove `o-pkgs`).

**Step 4.3 — Update `q-pkgs.nix`**

Remove `o-pkgs` from the `inherit` at line 17. Change:

```nix
  inherit (opts) nixgl-src is-linux lib o-pkgs;
```

to:

```nix
  inherit (opts) nixgl-src is-linux lib;
```

**Step 4.4 — Replace `o-pkgs` references in `q-pkgs.nix`**

Search for all `o-pkgs` usages in `q-pkgs.nix`:

- **Line 693:** `qlat-name-list-file = o-pkgs.runCommand` — If TODO 2 is
  already done, this line is gone. If not yet done, replace with
  `pkgs.runCommand` and add a local `pkgs = import-nixpkgs-wd {};` definition.
- **Line 728:** `use-gitee = use-gitee-wd;` — This references `use-gitee-wd`,
  not `o-pkgs`. No change needed.
- **Line 737 area:** No direct `o-pkgs` usage.

If `o-pkgs` is still needed in `q-pkgs.nix` (e.g., TODO 2 not yet done),
define a local alias:

```nix
  pkgs = import-nixpkgs-wd {};
```

and replace `o-pkgs` with `pkgs` everywhere in `q-pkgs.nix`.

---

## 5. Move `is-linux` to `q-pkgs.nix`

**File:** `options.nix:35`

```
# TODO: move to q-pkgs.nix
is-linux = (lib.lists.elem builtins.currentSystem lib.platforms.linux);
```

### Procedure

**Step 5.1 — Remove from `options.nix`**

Delete lines 35-36 from `options.nix`:

```nix
  # TODO: move to q-pkgs.nix
  is-linux = (lib.lists.elem builtins.currentSystem lib.platforms.linux);
```

**Step 5.2 — Remove from `options.nix` return set**

Change line 191 from:

```nix
  inherit nixgl-src is-linux lib;
```

to:

```nix
  inherit nixgl-src lib;
```

**Step 5.3 — Verify `is-linux` is not used inside `options.nix`**

Grep `options.nix` for `is-linux` — it only appears at line 36 (definition)
and line 191 (export). No internal usage. Safe to remove.

**Step 5.4 — Add to `q-pkgs.nix`**

In `q-pkgs.nix`, insert after line 18 (the `inherit` block):

```nix
  is-linux = (lib.lists.elem builtins.currentSystem lib.platforms.linux);
```

**Step 5.5 — Remove `is-linux` from `q-pkgs.nix` inherit**

Change line 17 from:

```nix
  inherit (opts) nixgl-src is-linux lib o-pkgs;
```

to:

```nix
  inherit (opts) nixgl-src lib o-pkgs;
```

(Adjust if TODO 4 already removed `o-pkgs`.)

**Step 5.6 — Verify downstream usage**

`is-linux` is used in `q-pkgs.nix` line 733 (`builtins.deepSeq` list) and
inside `mk-overlay` (search for `is-linux` in the file). All usages are in
`q-pkgs.nix`, so the move is safe.

---

## 6. Return processed values with original names from `options.nix`

**File:** `options.nix:189`

```
# TODO: For all args, return processed value with its original name (i.e., without -wd suffix)
```

Currently exported: `version-wd`, `use-gitee-wd`, `nixpkgs-wd`, `import-nixpkgs-wd`.

### Procedure

**Step 6.1 — Update `options.nix` return set**

Change lines 189-190 from:

```nix
  # TODO: For all args, return processed value with its original name (i.e., without -wd suffix)
  inherit version-wd use-gitee-wd nixpkgs-wd import-nixpkgs-wd;
```

to:

```nix
  version = version-wd;
  use-gitee = use-gitee-wd;
  nixpkgs = nixpkgs-wd;
  import-nixpkgs = import-nixpkgs-wd;
```

Keep the internal `-wd` variable definitions (lines 12-19) unchanged — they
are used throughout `options.nix` internally.

**Step 6.2 — Update `q-pkgs.nix` inherit**

Change line 16 from:

```nix
  inherit (opts) version-wd use-gitee-wd nixpkgs-wd import-nixpkgs-wd;
```

to:

```nix
  inherit (opts) version use-gitee nixpkgs import-nixpkgs;
```

**Step 6.3 — Update all references in `q-pkgs.nix`**

Search and replace in `q-pkgs.nix`:

| Line | Before | After |
|------|--------|-------|
| 20 | `version-pypi = "1.5";` | *(no change — not related)* |
| 35 | `use-gitee = use-gitee-wd;` | `use-gitee = use-gitee-wd;` *(inside mk-overlay, different scope — keep as-is)* |
| 723 | `version = version-wd;` | `version = version-wd;` *(inside q-pkgs attrset — keep as-is)* |
| 724 | `nixpkgs = nixpkgs-wd;` | `nixpkgs = nixpkgs-wd;` *(inside q-pkgs attrset — keep as-is)* |
| 728 | `use-gitee = use-gitee-wd;` | `use-gitee = use-gitee-wd;` *(inside q-pkgs attrset — keep as-is)* |

Wait — lines 723, 724, 728 are inside the `q-pkgs` return attrset and use
`version-wd`, `nixpkgs-wd`, `use-gitee-wd` as *values*. After TODO 6, these
variables are no longer in scope in `q-pkgs.nix` (they were inherited from
`opts`). They need to be updated:

| Line | Before | After |
|------|--------|-------|
| 723 | `version = version-wd;` | `version = version;` |
| 724 | `nixpkgs = nixpkgs-wd;` | `nixpkgs = nixpkgs;` |
| 728 | `use-gitee = use-gitee-wd;` | `use-gitee = use-gitee;` |

**Step 6.4 — Handle `mk-overlay` internal references**

Inside `mk-overlay` (line 25+), `use-gitee-wd` is used at line 35. After
TODO 6, the inherited name is `use-gitee`, so change line 35 from:

```nix
    use-gitee = use-gitee-wd;
```

to:

```nix
    use-gitee = use-gitee;
```

But this is a `let` binding that shadows the outer `use-gitee`. This is fine
in Nix (the inner `use-gitee` takes the value of the outer one). Verify no
other `-wd` references exist in `mk-overlay` by grepping for `-wd` in
`q-pkgs.nix`.

---

## 7. Define `import-nixgl` in `options.nix` and use it in `q-pkgs.nix`

Both files contain the same `import nixgl-src { pkgs = pkgs; }` expression.
Define a reusable `import-nixgl` function in `options.nix` and use it in both
files to eliminate duplication.

### Procedure

**Step 7.1 — Define `import-nixgl` in `options.nix`**

In `options.nix`, insert after `nixgl-src` definition (after line 33):

```nix
  import-nixgl = import nixgl-src;
```

**Step 7.2 — Use `import-nixgl` in `options.nix`**

Change line 58 from:

```nix
    nixgl = (import nixgl-src { pkgs = pkgs; }).auto.nixGLDefault;
```

to:

```nix
    nixgl = (import-nixgl { pkgs = pkgs; }).auto.nixGLDefault;
```

**Step 7.3 — Export `import-nixgl` from `options.nix`**

Add `import-nixgl` to the return set. Change line 191 from:

```nix
  inherit nixgl-src is-linux lib;
```

to:

```nix
  inherit nixgl-src import-nixgl is-linux lib;
```

**Step 7.4 — Inherit `import-nixgl` in `q-pkgs.nix`**

Update line 17 to include `import-nixgl`:

```nix
  inherit (opts) nixgl-src import-nixgl is-linux lib o-pkgs;
```

(Adjust if TODO 4 already removed `o-pkgs` or TODO 5 already removed
`is-linux`.)

**Step 7.5 — Use `import-nixgl` in `mk-overlay` in `q-pkgs.nix`**

Change line 178 from:

```nix
    nixgl = (import nixgl-src { pkgs = pkgs; }).auto.nixGLDefault;
```

to:

```nix
    nixgl = (import-nixgl { pkgs = pkgs; }).auto.nixGLDefault;
```

**Step 7.6 — Remove unused `nixgl-src` if applicable**

After this change, grep both files for `nixgl-src`. If it is no longer
referenced directly (only used via `import-nixgl`), remove it from
`options.nix` definition (lines 29-33) and from the return set. If it is
still used elsewhere, keep it.

---

## Suggested implementation order

1. **TODO 5** (move `is-linux`) — simplest, one variable, no dependencies
2. **TODO 3** (copy `force`) — simple copy, no breaking changes
3. **TODO 1** (move `nixpkgs-release`) — small, self-contained
4. **TODO 7** (define `import-nixgl`) — small, self-contained, no breaking changes
5. **TODO 6** (rename `-wd` exports) — touches both files, moderate
6. **TODO 4** (rename `o-pkgs` -> `pkgs`) — broader impact, needs care
7. **TODO 2** (move `options-list` section) — largest move, do last

After each TODO, run the parse check (Verification Step 1) to catch errors early.

---

## Verification

Run each step after completing all TODOs. If a step fails, bisect by reverting
one TODO at a time (in reverse order) to find the culprit.

### Step 1: Parse check (instantiation)

Confirms both files parse and evaluate without errors:
```bash
nix-instantiate --eval nixpkgs/options.nix -A options-default --strict >/dev/null
nix-instantiate nixpkgs/q-pkgs.nix --dry-run
```

### Step 2: Evaluated attribute spot-checks

Verify key attributes are accessible and have expected types:
```bash
# version string
nix-instantiate --eval nixpkgs/default.nix -A version

# qlat-name-list is a non-empty list
nix-instantiate --eval nixpkgs/q-pkgs.nix -A qlat-name-list --json | jq 'type == "array" and length > 0'

# options-default is an attrset with expected keys
nix-instantiate --eval nixpkgs/q-pkgs.nix -A options-default --json | jq 'has("use-grid-gpt", "use-cuda", "use-clang", "ngpu")'

# is-linux is a bool
nix-instantiate --eval nixpkgs/q-pkgs.nix -A options-default.ngpu --json
```

### Step 3: Name list comparison (regression guard)

The `qlat-name-list` must be identical before and after the refactor.
Capture it before making changes, then compare after:
```bash
# Before changes (run once, save baseline)
nix-instantiate --eval nixpkgs/q-pkgs.nix -A qlat-name-list --json > /tmp/qlat-names-before.json

# After changes (compare)
nix-instantiate --eval nixpkgs/q-pkgs.nix -A qlat-name-list --json > /tmp/qlat-names-after.json
diff /tmp/qlat-names-before.json /tmp/qlat-names-after.json
```

### Step 4: Build a minimal package variant (no CUDA)

Smoke-test that a full derivation graph evaluates and builds:
```bash
nix-build nixpkgs/q-pkgs.nix -A pkgs-std-ucxless.qlat-env -j 1 --cores 3
nix-build nixpkgs/q-pkgs.nix -A pkgs-std-ucxless.qlat-tests -j 1 --cores 3
```

### Step 5: Build a second variant (clang, no UCX)

Covers the `use-clang = true` code path:
```bash
nix-build nixpkgs/q-pkgs.nix -A pkgs-std-clang-ucxless.qlat-env -j 1 --cores 3
nix-build nixpkgs/q-pkgs.nix -A pkgs-std-clang-ucxless.qlat-tests -j 1 --cores 3
```

### Step 6: Build a pypi variant

Covers the `use-pypi` code path:
```bash
nix-build nixpkgs/q-pkgs.nix -A pkgs-pypi.qlat-env -j 1 --cores 3
```

### Step 7: Full CI-equivalent build (optional, slow)

Matches `.github/workflows/qlat.yml` — run only if time permits:
```bash
nix-build nixpkgs/q-pkgs.nix -A pkgs-ucxless.qlat-env -j 1 --cores 3
nix-build nixpkgs/q-pkgs.nix -A pkgs-ucxless.qlat-tests -j 1 --cores 3
```

### Step 8: Build-many script smoke test (optional)

Runs the same logic the team uses for batch builds:
```bash
cd nixpkgs && bash build-many-qlat-pkgs-core.sh --dry-run
```

### Step 9: CUDA kernel install and example run

End-to-end test: installs the CUDA Jupyter kernel and runs a real example
with cudasupport, covering the full CUDA code path:
```bash
make clean ; ./nixpkgs/install-py-local-kernel-with-nix-cuda.sh -j 4 --cores 15 ; ./nixpkgs/run-one-example-py.sh auto-contract-01 --cudasupport
```

### Step 10: Nix-based example run (`run-one-example-py.nix`)

Builds and runs a single example purely through Nix (no shell wrapper).
Tests both CPU and CUDA evaluation paths of `q-pkgs.nix`:
```bash
# CPU build
nix-build nixpkgs/run-one-example-py.nix --argstr testName auto-contract-01 -j 4 --cores 15

# CUDA build (requires GPU runner)
nix-build nixpkgs/run-one-example-py.nix --argstr testName auto-contract-01 --arg cudaSupport true -j 4 --cores 15
```
