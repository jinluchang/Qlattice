# `qlat.fields_io` — Shuffled Binary Field I/O

Source: `qlat/qlat/fields_io.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [`ShuffledFieldsWriter` Class](#shuffledfieldswriter-class)
   - [Constructor](#constructor)
   - [Methods](#writer-methods)
   - [Properties](#writer-properties)
3. [`ShuffledFieldsReader` Class](#shuffledfieldsreader-class)
   - [Constructor](#reader-constructor)
   - [Methods](#reader-methods)
   - [Properties](#reader-properties)
4. [`ShuffledBitSet` Class](#shufflebitset-class)
5. [Module-Level Functions](#module-level-functions)
6. [Examples](#examples)

---

## Overview

`fields_io` provides efficient parallel I/O for lattice fields using a
shuffled binary format. Data is distributed across MPI nodes and stored in a
directory structure where each node writes its portion to a separate file.
This avoids the bottleneck of single-node I/O while maintaining a simple
on-disk layout.

The two main classes are:

- **`ShuffledFieldsWriter`** — writes fields to a directory of shuffled
  binary files. Supports append mode.
- **`ShuffledFieldsReader`** — reads fields from a shuffled binary directory.

Fields are written/read by name (`fn`). Both `Field` and `SelectedField`
(sparse) objects are supported. A `FieldSelection` can be provided to
write/read only selected sites.

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
new_size_node = q.Coordinate([1, 1, 1, 1])

# Write a field
sfw = q.open_fields("/tmp/test_fields", "w", new_size_node)
f = q.Field("Complex", geo, 1)
sfw.write("my_field", f)
sfw.close()

# Read it back
sfr = q.open_fields("/tmp/test_fields", "r")
f2 = q.Field("Complex", geo, 1)
sfr.read("my_field", f2)
sfr.close()

q.end_with_mpi()
```

---

## `ShuffledFieldsWriter` Class

Writes lattice fields to a shuffled binary directory.

### <a id="constructor"></a> Constructor

#### `ShuffledFieldsWriter(path: str, new_size_node: Coordinate, is_append=False)`

Open a writer at `path` with the given I/O node layout.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Directory path for the output files |
| `new_size_node` | `Coordinate` | I/O node grid dimensions for the shuffled layout |
| `is_append` | `bool` | If `True`, append to existing data; otherwise create new |

---

### <a id="writer-methods"></a> Methods

#### `close()`

Close the writer and release resources.

#### `path() -> str`

Return the directory path.

#### `new_size_node() -> Coordinate`

Return the I/O node grid dimensions.

#### `list() -> list`

Return a list of field names currently stored.

#### `has(fn: str) -> bool`

Check if a field with name `fn` exists.

#### `flush()`

Flush buffered data to disk.

#### `write(fn: str, obj)`

Write a field object (`FieldBase` or `SelectedFieldBase`) under the name `fn`.
Internally dispatches to `obj.write_sfw_direct(self, fn)`.

#### `get_cache_sbs(fsel: FieldSelection) -> ShuffledBitSet`

Get or create a cached `ShuffledBitSet` for the given `FieldSelection`. Used
internally for sparse field I/O.

---

### <a id="writer-properties"></a> Properties

#### `__contains__(fn: str) -> bool`

Check membership via the `in` operator.

---

## `ShuffledFieldsReader` Class

Reads lattice fields from a shuffled binary directory.

### <a id="reader-constructor"></a> Constructor

#### `ShuffledFieldsReader(path: str, new_size_node=None)`

Open a reader at `path`.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Directory path to read from |
| `new_size_node` | `Coordinate` or `None` | I/O node grid; if `None`, auto-detect from stored metadata |

---

### <a id="reader-methods"></a> Methods

#### `close()`

Close the reader and release resources.

#### `path() -> str`

Return the directory path.

#### `new_size_node() -> Coordinate`

Return the I/O node grid dimensions.

#### `list() -> list`

Return a list of field names stored in the directory.

#### `has(fn: str) -> bool`

Check if a field with name `fn` exists.

#### `has_duplicates() -> bool`

Check if the stored data contains duplicate entries.

#### `is_sparse_field(fn: str) -> bool`

Return `True` if the named field is stored as a sparse (selected) field.

#### `read_as_char(fn: str)`

Read a field as raw bytes. Returns `SelectedFieldChar` for sparse fields,
`FieldChar` for dense fields, or `None` if the field does not exist.

#### `read(fn: str, obj)`

Read a field into `obj` (`FieldBase` or `SelectedFieldBase`). For
`SelectedField` objects, `obj.fsel` may be `None` before reading and will
be properly populated afterwards.

#### `get_cache_sbs(fsel: FieldSelection) -> ShuffledBitSet`

Get or create a cached `ShuffledBitSet` for the given `FieldSelection`.

---

### <a id="reader-properties"></a> Properties

#### `__contains__(fn: str) -> bool`

Check membership via the `in` operator.

---

## `ShuffledBitSet` Class

Maps a `FieldSelection` to the shuffled I/O node layout. Used internally by
both `ShuffledFieldsWriter` and `ShuffledFieldsReader` to determine which
sites each I/O node is responsible for.

#### `ShuffledBitSet(fsel: FieldSelection, new_size_node: Coordinate)`

| Parameter | Type | Description |
|---|---|---|
| `fsel` | `FieldSelection` | The field selection to map |
| `new_size_node` | `Coordinate` | The I/O node grid dimensions |

---

## Module-Level Functions

#### `open_fields(path: str, mode: str, new_size_node=None)`

Open a shuffled fields directory. Convenience function that creates a
`ShuffledFieldsReader` or `ShuffledFieldsWriter` depending on `mode`.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Directory path (or path to `geon-info.txt` inside it) |
| `mode` | `str` | `"r"` for read, `"w"` for write, `"a"` for append |
| `new_size_node` | `Coordinate` or `None` | Required for `"w"` mode; optional otherwise |

#### `list_fields(path: str, new_size_node=None) -> list`

Return a list of field names stored in the directory. Opens and closes a
reader internally.

#### `fields_build_index(path: str, new_size_node=None)`

Build the internal index for faster field lookup. Opens and closes a reader
internally.

#### `fields_has_duplicates(path: str, new_size_node=None) -> bool`

Check whether the stored data contains duplicate entries.

#### `properly_truncate_fields(path, is_check_all=False, is_only_check=False, new_size_node=None)`

Check and optionally truncate corrupted trailing data from the shuffled
files. Returns a list of successfully stored field names.

| Parameter | Type | Description |
|---|---|---|
| `path` | `str` | Directory path |
| `is_check_all` | `bool` | If `True`, check all fields (not just the last one) |
| `is_only_check` | `bool` | If `True`, only check without modifying |
| `new_size_node` | `Coordinate` | I/O node grid |

#### `truncate_fields(path: str, fns_keep: list, new_size_node=None)`

Truncate the stored data to keep only the fields listed in `fns_keep`.
The names must be in the same order as they appear on disk. All names in
`fns_keep` must already exist in the directory.

#### `check_fields(path: str, is_check_all=True, new_size_node=None) -> list`

Return a list of field names that can be read successfully from the
directory.

#### `check_compressed_eigen_vectors(path: str) -> bool`

Check whether compressed eigenvector data can be read. Returns `True` if
there is a problem, `False` if the data is OK.

#### `eigen_system_repartition(new_size_node, path, path_new="")`

Repartition a stored eigensystem to a new I/O node layout.

#### `show_all_shuffled_fields_writer() -> str`

Return debug information about all currently open `ShuffledFieldsWriter`
instances.

---

## Examples

### Writing and Reading Fields

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
new_size_node = q.Coordinate([1, 1, 1, 1])

# Create fields
f1 = q.Field("Complex", geo, 1)
f2 = q.Field("Complex", geo, 1)

# Write fields
sfw = q.open_fields("/tmp/fields_test", "w", new_size_node)
sfw.write("field_a", f1)
sfw.write("field_b", f2)
sfw.close()

# List stored fields
names = q.list_fields("/tmp/fields_test", new_size_node)
print(f"Stored fields: {names}")  # ["field_a", "field_b"]

# Read fields back
sfr = q.open_fields("/tmp/fields_test", "r")
f1_read = q.Field("Complex", geo, 1)
sfr.read("field_a", f1_read)
print("field_a" in sfr)  # True
sfr.close()

q.end_with_mpi()
```

### Append Mode

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
new_size_node = q.Coordinate([1, 1, 1, 1])

# Initial write
sfw = q.open_fields("/tmp/fields_append", "w", new_size_node)
sfw.write("field_1", q.Field("Complex", geo, 1))
sfw.close()

# Append more fields
sfw = q.open_fields("/tmp/fields_append", "a", new_size_node)
sfw.write("field_2", q.Field("Complex", geo, 1))
sfw.close()

names = q.list_fields("/tmp/fields_append", new_size_node)
print(f"Stored fields: {names}")  # ["field_1", "field_2"]

q.end_with_mpi()
```

### Checking and Truncating Fields

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
new_size_node = q.Coordinate([1, 1, 1, 1])

# Write some fields
sfw = q.open_fields("/tmp/fields_check", "w", new_size_node)
sfw.write("f1", q.Field("Complex", geo, 1))
sfw.write("f2", q.Field("Complex", geo, 1))
sfw.write("f3", q.Field("Complex", geo, 1))
sfw.close()

# Check which fields are readable
ok = q.check_fields("/tmp/fields_check", new_size_node=new_size_node)
print(f"Readable fields: {ok}")

# Truncate to keep only f1 and f2
q.truncate_fields("/tmp/fields_check", ["f1", "f2"], new_size_node)
names = q.list_fields("/tmp/fields_check", new_size_node)
print(f"After truncate: {names}")  # ["f1", "f2"]

q.end_with_mpi()
```

### Using open_fields with geon-info.txt Path

```python
import qlat as q

q.begin_with_mpi([[1, 1, 1, 1]])

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
new_size_node = q.Coordinate([1, 1, 1, 1])

path = "/tmp/fields_geon"

# Write
sfw = q.open_fields(path, "w", new_size_node)
sfw.write("test", q.Field("Complex", geo, 1))
sfw.close()

# Read using geon-info.txt path (the /geon-info.txt suffix is stripped)
sfr = q.open_fields(path + "/geon-info.txt", "r")
print(sfr.list())
sfr.close()

q.end_with_mpi()
```
