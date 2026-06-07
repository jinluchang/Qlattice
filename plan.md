# Plan: Add `qlat-cuda-tests` argument to `many-qlat-pkgs.nix`

## Problem

`nixpkgs/many-qlat-pkgs.nix:9` has a TODO to add a `qlat-cuda-tests` argument that controls which CUDA-related test variants are built.

## Current behavior

`many-qlat-pkgs.nix` iterates over all `version-list` entries and all `qlat-name-list` entries from `q-pkgs`, building `qlat-tests` and `qlat-env` for every combination. There is no way to filter CUDA variants.

## CUDA variant identification

From `options.nix:197-199`, name suffixes are built incrementally:
- `use-cuda-software` → `-cu`
- `use-cuda` → appends `da` → `-cuda`
- `use-cudasupport` → appends `support` → `-cudasupport`

Only `-cuda` and `-cudasupport` should be filtered (not `-cu` alone, which is just `use-cuda-software`).

A variant is "CUDA test" if its name contains `-cuda` (which matches both `-cuda` and `-cudasupport` but not `-cu` alone):
```nix
builtins.match ".*-cuda.*" name != null
```

## Changes

### File: `nixpkgs/many-qlat-pkgs.nix`

1. **Add the new argument** (line 9): Replace the TODO comment with a proper argument declaration:
   ```nix
   qlat-cuda-tests ? null,
   ```

2. **Add helper functions** to determine if tests/env should be included:
   ```nix
   is-cuda-name = name: builtins.match ".*-cuda.*" name != null;
   include-tests = if qlat-cuda-tests == null then (name: true)
     else if qlat-cuda-tests == "none" then (name: ! is-cuda-name name)
     else if qlat-cuda-tests == "only" then (name: is-cuda-name name)
     else builtins.throw "qlat-cuda-tests must be null, \"none\", or \"only\", got: ${builtins.toString qlat-cuda-tests}";
   include-env = qlat-cuda-tests != "only";
   ```

3. **Filter only `qlat-tests`** in `mk-version-entries` (around line 39-49): One function, one list, conditionally include tests:
   ```nix
   name-list = builtins.map (n: "q-pkgs${n}") (
     if qlat-name-list != null
     then builtins.filter (n: builtins.elem n qlat-name-list) q-pkgs.qlat-name-list
     else q-pkgs.qlat-name-list
   );
   ver-suffix = if version == "" then "" else "-${builtins.replaceStrings ["."] ["-"] version}";
   mk-name-entries = name: (if include-env then {
     "${name}${ver-suffix}-qlat-env" = q-pkgs.${name}.qlat-env;
   } else {}) // (if include-tests name then {
     "${name}${ver-suffix}-qlat-tests" = q-pkgs.${name}.qlat-tests;
   } else {});
   in builtins.foldl' (s: v: s // v) {} (builtins.map mk-name-entries name-list);
   ```

4. **Assert validation** at the top of the `let` block for safety:
   ```nix
   assert builtins.elem qlat-cuda-tests [ null "none" "only" ];
   ```

## Resulting behavior

| `qlat-cuda-tests` | `qlat-tests` | `qlat-env` |
|---|---|---|
| `null` (default) | All variants | All variants |
| `"none"` | Exclude CUDA variants | All variants |
| `"only"` | Only CUDA variants | None |

## Verification

Syntax check:
```bash
nix-instantiate --parse nixpkgs/many-qlat-pkgs.nix
```

Evaluate without error for each option:
```bash
# Default (all)
nix-instantiate nixpkgs/many-qlat-pkgs.nix
# No CUDA
nix-instantiate nixpkgs/many-qlat-pkgs.nix --arg qlat-cuda-tests '"none"'
# Only CUDA
nix-instantiate nixpkgs/many-qlat-pkgs.nix --arg qlat-cuda-tests '"only"'
```

Spot-check specific attributes to verify filtering works:
```bash
# Default: both exist
nix-instantiate nixpkgs/many-qlat-pkgs.nix -A 'q-pkgs-cuda-qlat-tests'
nix-instantiate nixpkgs/many-qlat-pkgs.nix -A 'q-pkgs-cuda-qlat-env'
# "none": tests absent, env present
nix-instantiate nixpkgs/many-qlat-pkgs.nix --arg qlat-cuda-tests '"none"' -A 'q-pkgs-cuda-qlat-env'
# "only": tests present, env absent
nix-instantiate nixpkgs/many-qlat-pkgs.nix --arg qlat-cuda-tests '"only"' -A 'q-pkgs-cuda-qlat-tests'
```
