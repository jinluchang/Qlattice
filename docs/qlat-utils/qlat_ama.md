# `qlat_utils.ama` — All-mode Averaging (AMA) for Multi-accuracy Measurements

Source: `qlat-utils/qlat_utils/ama.py`

> **Note:** Update this document when updating the source file.

## Outline

- `AmaVal` — container for sloppy value plus accuracy corrections
- `mk_ama_val` — construct an AMA value from a list of measurements
- `ama_extract` / `ama_extract_ama_val` — produce the variance-reduced final estimate
- `ama_apply1`, `ama_apply2`, `ama_apply` — propagate operations through AMA values
- `ama_list` — collect multiple AMA values into a list
- `ama_msc_mult`, `ama_msc_add` — arithmetic helpers used by `AmaVal.__mul__` and `__add__`
- `merge_description_dict` — merge source-specification dictionaries

## Overview

All-mode Averaging (AMA) is a variance-reduction technique used in lattice
QCD.  Measurements are taken at multiple accuracy levels (e.g. sloppy,
fine, ultra-fine) with known sampling probabilities.  The `ama` module
provides a container, `AmaVal`, that stores the sloppy result together
with all correction pairs `(value, description_dict)` and propagates
arbitrary binary operations through all valid correction combinations.
`ama_extract` then combines them into a single unbiased estimate.

The **description dictionary** maps each `source_specification` key to a
tuple `(accuracy_level, probability)` and ensures that corrections from
different sources are only combined when their specifications are
compatible.

## Detailed Sections

### `AmaVal`

```python
class AmaVal:
    val           # sloppy result (may be None after operations)
    corrections   # [(value, {source_spec: (level, prob), ...}), ...]
```

Supports `+`, `*` with scalars or other `AmaVal` instances.  The
operators delegate to `ama_msc_add` and `ama_msc_mult` which in turn
call `ama_apply2`.

### Constructing AMA Values

`mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list)`:

- `val` — the sloppy result (must equal `val_list[0]`)
- `source_specification` — hashable key uniquely identifying the propagator
  source, e.g. `("point", (12, 2, 3, 4))`
- `val_list` — list of measurements at increasing accuracy; `None` entries
  are skipped
- `rel_acc_list` — integer accuracy levels (0 = sloppy, 1 = fine, ...)
- `prob_list` — probability of each accuracy level being sampled

If only one non-`None` value is provided the function returns the plain
value (no wrapping).

### Applying Operations

| Function | Description |
|----------|-------------|
| `ama_apply1(f, x)` | Apply unary `f` to every correction value |
| `ama_apply2(f, x, y)` | Apply binary `f` to all compatible correction pairs |
| `ama_apply(f, *args)` | N-ary version; collects args via `ama_list` then applies `f` |

When both operands are `AmaVal`, the correction dictionaries are merged
with `merge_description_dict`: dictionaries are compatible iff all shared
keys map to the same `(level, prob)` pair.

### Extracting the Final Result

`ama_extract(x, *, is_sloppy=False)`:

- If `x` is not an `AmaVal`, returns `x` unchanged.
- If `is_sloppy=True`, returns only the sloppy component.
- Otherwise, recursively computes the correction sum:

```
result = sloppy
for each source key k:
    for accuracy level i >= 1:
        diff = val(level_i) - val(level_{i-1})
        result += diff / prob_i
```

This yields an unbiased estimate with reduced variance compared to the
sloppy measurement alone.

### Merging Description Dictionaries

`merge_description_dict(d1, d2)` returns a merged dictionary if
compatible, or `None` if any shared key has conflicting values.
Compatibility ensures that two correction terms can be combined
(e.g. when multiplying AMA values from different sources).

### Utility

`ama_counts(x)` returns the number of correction terms (1 for plain
values).  Useful for estimating computational cost.

## Examples

### Basic AMA construction and extraction

```python
import qlat_utils as q

# Sloppy measurement = 1.0, fine = 1.01, ultra-fine = 1.011
v = q.mk_ama_val(
    1.0,
    ("point", (0, 0, 0, 0)),
    [1.0, 1.01, 1.011],
    [0, 1, 2],
    [1.0, 0.1, 0.02],
)
# v is an AmaVal; extract the corrected result
result = q.ama_extract(v)
```

### Arithmetic on AMA values

```python
import qlat_utils as q

v1 = q.mk_ama_val(
    1.0, "src_a",
    [1.0, 1.01], [0, 1], [1.0, 0.1],
)
v2 = q.mk_ama_val(
    2.0, "src_b",
    [2.0, 2.01], [0, 1], [1.0, 0.1],
)

# Operations propagate through corrections automatically
v_sum = v1 + v2
v_prod = v1 * v2

print(q.ama_extract(v_sum))   # corrected sum
print(q.ama_extract(v_prod))  # corrected product
```

### Mixed AMA and plain values

```python
import qlat_utils as q

v = q.mk_ama_val(
    1.0, "src",
    [1.0, 1.01], [0, 1], [1.0, 0.1],
)

# Multiplying by a plain scalar works transparently
scaled = v * 3.0
print(q.ama_extract(scaled))
```

### Extracting the sloppy-only result

```python
import qlat_utils as q

v = q.mk_ama_val(
    1.0, "src",
    [1.0, 1.01, 1.011], [0, 1, 2], [1.0, 0.1, 0.02],
)
sloppy = q.ama_extract(v, is_sloppy=True)
print("Sloppy only:", sloppy)   # 1.0
```

### Counting correction terms

```python
import qlat_utils as q

plain = 3.14
print(q.ama_counts(plain))  # 1

v = q.mk_ama_val(1.0, "src", [1.0, 1.01, 1.011], [0, 1, 2], [1.0, 0.1, 0.02])
print(q.ama_counts(v))      # 3
```
