# `qlat_utils.json` — JSON Serialization for NumPy and Numeric Types

## Outline

1. [Overview](#overview)
2. [Supported Types](#supported-types)
3. [API](#api)
4. [Encoding Format](#encoding-format)
5. [Examples](#examples)

---

## Overview

The `qlat_utils.json` module provides `json_dumps` and `json_loads` functions
that extend Python's standard `json` module to handle NumPy and extended numeric
types commonly used in lattice QCD computations.

Standard `json.dumps` raises `TypeError` on `numpy.float32`, `numpy.ndarray`,
`complex`, and similar types. This module encodes them transparently and
decodes them back to their original types, ensuring faithful round-tripping.

```python
import qlat_utils as q
s = q.json_dumps({"value": np.float32(0.5)})
obj = q.json_loads(s)   # recovers np.float32(0.5)
```

---

## Supported Types

| Python / NumPy Type | Round-Trips As |
|---|---|
| `complex` | `complex` |
| `numpy.complex64` | `numpy.complex64` |
| `numpy.complex128` | `numpy.complex128` |
| `numpy.float32` | `numpy.float32` |
| `numpy.int32` | `numpy.int32` |
| `numpy.int64` | `numpy.int64` |
| `numpy.ndarray` | `numpy.ndarray` (via `tolist()`) |
| `range` | `range` |

All standard JSON types (`int`, `float`, `str`, `list`, `dict`, `None`, `bool`)
are handled by the standard encoder and pass through unchanged.

---

## API

### `json_dumps(obj, *, indent=2) -> str`

Serialize `obj` to a JSON string. Extended types are encoded with an
`__extended_json_type__` marker (see [Encoding Format](#encoding-format)).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `obj` | any | — | Object to serialize |
| `indent` | `int \| str \| None` | `2` | JSON indentation; `None` for compact output |

```python
s = q.json_dumps({"mass": np.float32(0.5), "corr": np.array([1, 2, 3])})
compact = q.json_dumps(data, indent=None)
```

### `json_loads(s) -> any`

Deserialize a JSON string back into Python/NumPy objects. Extended types
(encoded with `__extended_json_type__`) are restored to their original types.

| Parameter | Type | Description |
|---|---|---|
| `s` | `str` | JSON string to deserialize |

```python
obj = q.json_loads(s)
```

---

## Encoding Format

Each extended type is stored as a JSON object with an `__extended_json_type__`
key identifying the original type. The decoder uses this key to restore the
correct Python/NumPy type automatically.

| Type | Encoded Representation |
|---|---|
| `complex` | `{"real": r, "imag": i, "__extended_json_type__": "complex"}` |
| `complex64` | `{"real": r, "imag": i, "__extended_json_type__": "complex64"}` |
| `complex128` | `{"real": r, "imag": i, "__extended_json_type__": "complex128"}` |
| `float32` | `{"value": v, "__extended_json_type__": "float32"}` |
| `int32` | `{"value": v, "__extended_json_type__": "int32"}` |
| `int64` | `{"value": v, "__extended_json_type__": "int64"}` |
| `ndarray` | `{"value": [...], "__extended_json_type__": "ndarray"}` |
| `range` | `{"start": s, "stop": e, "step": t, "__extended_json_type__": "range"}` |

Plain Python types (`int`, `float`, `str`, `list`, `dict`, `None`, `bool`) are
encoded by the standard JSON encoder and contain no `__extended_json_type__`
marker.

---

## Examples

### Basic Round-Trip

```python
import qlat_utils as q
import numpy as np

data = {
    "mass": np.float32(0.5),
    "energy": complex(1.2, -0.3),
    "propagator": np.array([[1+0j, 2+0j], [3+0j, 4+0j]]),
    "sweeps": np.int64(1000),
    "step_range": range(0, 100, 5),
}

s = q.json_dumps(data)
print(s)

restored = q.json_loads(s)
print(type(restored["mass"]))        # <class 'numpy.float32'>
print(type(restored["energy"]))      # <class 'complex'>
print(type(restored["propagator"]))  # <class 'numpy.ndarray'>
print(type(restored["step_range"]))  # <class 'range'>
```

### Saving and Loading Parameters

```python
import qlat_utils as q
import numpy as np

params = {
    "lattice_size": [16, 16, 16, 32],
    "beta": np.float32(6.0),
    "mass": np.float32(0.05),
    "n_therm": np.int32(100),
    "n_traj": np.int64(2000),
}

# Write to file
with open("params.json", "w") as f:
    f.write(q.json_dumps(params))

# Read back
with open("params.json") as f:
    params = q.json_loads(f.read())
```

### Compact (No Indentation)

```python
s = q.json_dumps(data, indent=None)
```

### Nested Structures

```python
data = {
    "config_0": {
        "plaquette": np.float64(0.58),
        "correlator": np.array([1.0, 0.8, 0.6]),
    },
    "config_1": {
        "plaquette": np.float64(0.59),
        "correlator": np.array([1.0, 0.79, 0.58]),
    },
}

s = q.json_dumps(data)
restored = q.json_loads(s)
# All numpy types preserved through nesting
```
