# Adding a new test script

## 1. Create the test file

Create `<name>.py` with this skeleton:

```python
#!/usr/bin/env python3

import qlat as q

size_node_list = [
    [1, 1, 1, 1],
]

q.begin_with_mpi(size_node_list)

q.json_results_append("test description")
# ... test logic using q.json_results_append(...) to record results ...

q.check_log_json(__file__, check_eps=1e-14)
q.timer_display()
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
```

- Use `q.json_results_append(str)` to record test output for log comparison.
- End with `q.displayln_info("CHECK: finished successfully.")` — the trailing dot is required.

## 2. Register in the Makefile

Add `<name>.log` to the `tests` list in `Makefile` (keep alphabetical order):

```makefile
tests = \
    ...
    <name>.log
```

## 3. Recompile (if source code changed)

If the test exercises new or modified library code (`qlat-utils`, `qlat`, etc.),
recompile before running:

```bash
bash scripts/qlat-all.sh
```

If meson errors occur during configuration, do a clean build first:

```bash
bash scripts/qlat-clean-build.sh
bash scripts/qlat-all.sh
```

## 4. Generate the initial reference logs

```bash
cd examples-py
rm -rf <name>.py.p <name>.log <name>.log.json
make <name>.log
```

The first run will show "failed" (no reference yet), but it generates both files:
- `<name>.log` — `CHECK:` lines for final validation
- `<name>.log.json` — detailed test output for `q.check_log_json()`

## 5. Verify the test passes

```bash
rm -rf <name>.py.p && touch <name>.py && make <name>.log
```

Should print `passed`. Do not delete `<name>.log` or `<name>.log.json` — they
are the committed reference files used for cross-check.

## How it works

The Makefile rule `%.log: %.py`:
1. Copies `%.py` and `%.log.json` into a temp `.p/` directory.
2. Runs the script via `mpiexec -n 2` with MPI geometry `1.1.1.2`.
3. Greps `CHECK:` / `INFO:` / `WARNING` lines from the full output into `log.txt`.
4. Extracts `CHECK:` lines into `log.check.txt.new` and diffs against the
   committed `%.log` reference file.
5. Copies `log.txt` → `%.log` and `%.log.json.new` → `%.log.json`.

## Running all tests

```bash
cd examples-py
make run
```
