# field_trunc_inv_acc_list Parameter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use compose:subagent (recommended) or compose:execute to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `inv_acc` configurable in `run_prop_wsrc_truncated_save()` via a parameter `measurement.field_trunc_inv_acc_list` mapping `inv_type` index to the desired accuracy.

**Architecture:** Add a 2-element list parameter `field_trunc_inv_acc_list` where index 0 = light, index 1 = strange. Read it in `run_prop_wsrc_truncated_save()` instead of hardcoding `inv_acc = 0`.

**Tech Stack:** Python, qlat params API (`set_param`/`get_param`)

## Global Constraints

- Default to `[0, 0]` for all ensembles to preserve existing behavior
- Parameter goes under `measurement` namespace, consistent with `field_trunc_half_width_list` and `field_trunc_t_size_divisor`

---

### Task 1: Add param and use it in run_prop_wsrc_truncated_save

**Files:**
- Modify: `examples-py-gpt/gpt-qlat-data-gen-truncated-wsrc-prop.py`

**Interfaces:**
- Consumes: `get_param(job_tag, "measurement")` — existing params
- Produces: `inv_acc` integer read from param

- [ ] **Step 1: Replace hardcoded inv_acc with param lookup**

In `run_prop_wsrc_truncated_save` (line 302), replace:
```python
    inv_acc = 0
```
with:
```python
    inv_acc = get_param(job_tag, "measurement", "field_trunc_inv_acc_list")[inv_type]
```

- [ ] **Step 2: Set default param for all ensembles**

After each existing `set_param(job_tag, "measurement", "field_trunc_t_size_divisor")` line, add:
```python
set_param(job_tag, "measurement", "field_trunc_inv_acc_list")([0, 0])
```
This applies to: 24D (line 421), 32Dfine (line 426), 64I (line 431), test-4nt16-checker (line 499).

- [ ] **Step 3: Run test to verify**

```bash
./nixpkgs/run-one-example-py-gpt.py gpt-qlat-data-gen-truncated-wsrc-prop
```

Expected: existing CHECK lines match, behavior unchanged (defaults to `[0, 0]`).
