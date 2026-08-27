# Command-Line Argument Handling

All Qlattice Python scripts use **standard `argparse`** for command-line argument
parsing. The previous `q.get_arg()` / `q.get_option()` helpers are deprecated in
favor of the standard library.

## Core Rules

1. **Use underscores in option names.**
   Option names must use underscores (`_`), not hyphens (`-`), so that the
   `argparse`-generated attribute name matches the option name directly:

   ```python
   parser.add_argument("--no_inversion", ...)   # ✓  attribute: no_inversion
   parser.add_argument("--no-inversion", ...)   # ✗  attribute: no_inversion (ambiguous)
   parser.add_argument("--random_permute_job_tag_traj_list", ...)
   ```

2. **Do not use `dest`.**
   Since option names use underscores, `argparse`'s default attribute name
   (`--foo_bar` → `foo_bar`) always matches the option. The `dest` parameter
   should never be needed.

3. **Use `parse_known_args()` for MPI safety.**
   Scripts launched via `mpiexec` may receive extra arguments from the MPI
   runtime. Always use `parse_known_args()` instead of `parse_args()` to
   ignore unknown arguments gracefully:

   ```python
   args, _ = parser.parse_known_args()
   ```

4. **Use the raw `no_*` flags directly in logic.**
   Do **not** convert negative flags (e.g. `--no_inversion`) into positive
   booleans (e.g. `is_performing_inversion`). Invert the logic inline where
   needed:

   ```python
   if not sys_args.no_inversion:     # ✓
   if sys_args.is_performing_inversion:  # ✗
   ```

5. **Post-process list/composite arguments in `parse_args()`.**
   Simple transformations (e.g. splitting a comma-separated string into a list)
   belong inside the `parse_args()` function:

   ```python
   args, _ = parser.parse_known_args()
   args.job_tag_list = args.job_tag_list.split(",")
   return args
   ```

6. **Call `parse_args()` inside `if __name__ == "__main__":`.**
   The `parse_args()` function is defined at module level but called only when
   the script is executed directly, not when imported:

   ```python
   if __name__ == "__main__":
       sys_args = parse_args()
       job_tag_list = sys_args.job_tag_list
   ```

## Reference Example

```python
import argparse
import sys

# ... imports and function definitions ...

def parse_args():
    parser = argparse.ArgumentParser(
        description="Script description."
    )
    parser.add_argument(
        "--job_tag_list",
        type=str,
        default="test-4nt8-checker",
        help="Comma-separated list of job tags (default: 'test-4nt8-checker')",
    )
    parser.add_argument(
        "--no_inversion",
        action="store_true",
        default=False,
        help="Skip the inversion step",
    )
    parser.add_argument(
        "--no_contract",
        action="store_true",
        default=False,
        help="Skip the contraction step",
    )
    parser.add_argument(
        "--random_permute_job_tag_traj_list",
        action="store_true",
        default=False,
        help="Randomly permute the job_tag/traj list",
    )
    args, _ = parser.parse_known_args()
    args.job_tag_list = args.job_tag_list.split(",")
    return args

if __name__ == "__main__":
    sys_args = parse_args()
    # use sys_args.no_inversion, sys_args.no_contract, etc.
    if not sys_args.no_inversion:
        run_inversion()
    if not sys_args.no_contract:
        run_contraction()
```

## Migration from `q.get_arg()`

| Old pattern | New pattern |
|---|---|
| `q.get_arg("--job_tag_list", default="tag1").split(",")` | `sys_args.job_tag_list` (split in `parse_args()`) |
| `q.get_arg("--no-inversion", default=None) is None` | `not sys_args.no_inversion` |
| `q.get_option("--random-permute-job_tag-traj-list")` | `sys_args.random_permute_job_tag_traj_list` |
