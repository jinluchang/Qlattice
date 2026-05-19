# `qlat_utils.qar` — QAR Archive Format, QFile, and File I/O

Source: `qlat-utils/qlat_utils/qar.pyx`

> **Note:** Update this document when updating the source file.

## Outline

1. [Overview](#overview)
2. [QFile (Virtual File Abstraction)](#qfile-virtual-file-abstraction)
3. [QarFile (QAR Archive Reader/Writer)](#qarfile-qar-archive-readerwriter)
4. [QAR Archive Operations](#qar-archive-operations)
5. [Basic File I/O](#basic-file-io)
6. [Single-Node Operations](#single-node-operations)
7. [Sync-Node Operations](#sync-node-operations)
8. [Examples](#examples)

---

## Overview

The `qlat_utils.qar` module provides:

- **QFile** — a virtual file abstraction supporting in-memory (string) and on-disk (CFile) backends, with transparent CRC-32 checksumming.
- **QarFile** — a QAR archive reader/writer that packs multiple named files into a single archive with an index, supporting multi-volume archives.
- **QAR archive utilities** — create, extract, build index, list, and verify archives.
- **File I/O helpers** — `qcat`, `qcat_bytes`, `qtouch`, `qappend`, `qload_datatable`, `compute_crc32`.
- **Single-node variants** — `*_info` functions that only operate on node 0.
- **Sync-node variants** — `*_sync_node` functions that operate on node 0 and broadcast results to all nodes.

---

## QFile (Virtual File Abstraction)

`QFile` is a virtual file object that can wrap either a regular file on disk or an in-memory string buffer.

### Constructor

```python
QFile(ftype="CFile", path="", mode="r")
QFile(ftype, path, mode, content)
QFile(qfile, q_offset_start, q_offset_end)
```

Always use keyword arguments.

| Parameter | Type | Description |
|-----------|------|-------------|
| `ftype` | `str` | File type: `"CFile"` (disk) or `"String"` (in-memory). Default: `"CFile"`. |
| `path` | `str` | File path or name identifier. |
| `mode` | `str` | Mode: `"r"` (read), `"w"` (write), `"a"` (append). Default: `"r"`. |
| `content` | `str` or `bytes` | Initial content for in-memory files (requires `ftype="String"`). |
| `qfile` | `QFile` | Parent `QFile` to create a sub-range view. |
| `q_offset_start` | `int` | Start offset for sub-range view. |
| `q_offset_end` | `int` | End offset for sub-range view (-1 means end of file). |

### Methods

| Method | Description |
|--------|-------------|
| `path()` | Return the file path as a string. |
| `mode()` | Return the file mode as a string. |
| `close()` | Close the file. |
| `null()` | Return `True` if the QFile is null (uninitialised). |
| `eof()` | Return `True` if at end of file. |
| `tell()` | Return the current position in the file. |
| `flush()` | Flush any buffered data. |
| `seek_set(offset)` | Seek to absolute `offset` from the beginning. |
| `seek_end(offset)` | Seek to `offset` relative to the end. |
| `seek_cur(offset)` | Seek to `offset` relative to the current position. |
| `content()` | Return the full file content as a `str`. |
| `content_bytes()` | Return the full file content as `bytes`. |
| `size()` | Return the file size. |
| `remaining_size()` | Return the number of bytes from the current position to the end. |
| `getline()` | Read the next line as a string. |
| `getlines()` | Read all remaining lines as a list of strings. |
| `qcat()` | Read and return the entire file as `str` (alias for `content()`). |
| `qcat_bytes()` | Read and return the entire file as `bytes`. |
| `write(data)` | Write data. Accepts `QFile`, `str`, `bytes`, `list`, or `tuple`. |
| `compute_crc32()` | Compute the CRC-32 checksum of the file content. |

### Creating QFile from Path

```python
open_qfile(path, mode)
open_qfile_str(path, mode, content=None)
```

`open_qfile` opens a disk file. `open_qfile_str` opens an in-memory string file, optionally with initial `content`.

---

## QarFile (QAR Archive Reader/Writer)

`QarFile` provides read/write access to QAR archive files. A QAR archive packs multiple named files into a single file with an index, similar to `tar` but with a simpler format and CRC-32 support.

### Constructor

```python
QarFile(path=path, mode=mode)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `str` | Path to the QAR archive file. |
| `mode` | `str` | Mode: `"r"` (read), `"w"` (write), `"a"` (append). |

### Methods

| Method | Description |
|--------|-------------|
| `path()` | Return the archive path. |
| `mode()` | Return the archive mode. |
| `close()` | Close the archive. |
| `null()` | Return `True` if the QarFile is null. |
| `flush()` | Flush buffered data. |
| `list()` | List all entries in the archive. |
| `has_regular_file(fn)` | Check if entry `fn` is a regular file. |
| `has(fn)` | Check if entry `fn` exists (file or directory). |
| `fn in qarfile` | Use `in` operator to check membership. |
| `read(fn)` | Read entry `fn` and return a `QFile`. |
| `read_data(fn)` | Read entry `fn` and return content as `str`. |
| `read_data_bytes(fn)` | Read entry `fn` and return content as `bytes`. |
| `read_info(fn)` | Read metadata info string for entry `fn`. |
| `verify_index()` | Verify the archive index integrity. |
| `write(fn, info, data, skip_if_exist=False)` | Write entry `fn` with metadata `info` and `data`. |
| `show_index()` | Return a string representation of the archive index. |
| `read_index(qar_index_content)` | Load index from a string. |
| `index_size_saved()` | Return the size of the saved index. |
| `index_size()` | Return the size of the in-memory index. |
| `save_index(max_diff=0)` | Save the index to disk. |

### Opening a QAR Archive

```python
open_qar(path, mode)
open_qar_info(path, mode)
```

`open_qar` opens a QAR archive. `open_qar_info` opens only on node 0 and returns a no-op `Gobble()` on other nodes (for MPI parallelism).

### Configuration

```python
get_qar_multi_vol_max_size()
set_qar_multi_vol_max_size(size)
```

Controls the maximum size of a single QAR volume in bytes. When set, large archives are automatically split into multiple volume files. QAR never splits a single file across volumes.

---

## QAR Archive Operations

```python
properly_truncate_qar_file(path)
does_regular_file_exist_qar(path)
does_file_exist_qar(path)
qar_build_index(path_qar)
qar_create(path_qar, path_folder, is_remove_folder_after=False)
qar_extract(path_qar, path_folder, is_remove_qar_after=False)
qcopy_file(path_src, path_dst)
list_qar(path_qar)
```

| Function | Description |
|----------|-------------|
| `properly_truncate_qar_file` | Truncate a QAR file to remove any trailing garbage. |
| `does_regular_file_exist_qar` | Check if a path exists as a regular file (handles both QAR and regular files). |
| `does_file_exist_qar` | Check if a path exists as a file or directory (handles both QAR and regular files). |
| `qar_build_index` | Create a `.idx` index file for a QAR archive (speeds up random access). |
| `qar_create` | Create a QAR archive from a directory. |
| `qar_extract` | Extract a QAR archive to a directory. |
| `qcopy_file` | Copy a file (handles both QAR and regular files). |
| `list_qar` | List entries in a QAR archive by path. |

---

## Basic File I/O

```python
qcat(path)
qcat_bytes(path)
qtouch(path, content=None)
qappend(path, content)
qload_datatable(path, is_par=False)
compute_crc32(path)
```

| Function | Description |
|----------|-------------|
| `qcat` | Read file contents as `str`. |
| `qcat_bytes` | Read file contents as `bytes`. |
| `qtouch` | Create an empty file, or write `content` (supports `str`, `bytes`, `list`, `tuple`). |
| `qappend` | Append `content` to a file (supports `str`, `bytes`, `list`, `tuple`). |
| `qload_datatable` | Load a text datatable (columns of numbers). |
| `compute_crc32` | Compute CRC-32 checksum of a file. |

---

## Single-Node Operations

These functions operate only on MPI node 0:

```python
qar_build_index_info(path_qar)
qar_create_info(path_qar, path_folder, is_remove_folder_after=False)
qar_extract_info(path_qar, path_folder, is_remove_qar_after=False)
qcopy_file_info(path_src, path_dst)
qtouch_info(path, content=None)
qappend_info(path, content)
check_all_files_crc32_info(path)
```

---

## Sync-Node Operations

These functions execute on node 0 and broadcast the result to all MPI nodes:

```python
does_regular_file_exist_qar_sync_node(path)
does_file_exist_qar_sync_node(path)
qar_create_sync_node(path_qar, path_folder, is_remove_folder_after=False)
qar_extract_sync_node(path_qar, path_folder, is_remove_qar_after=False)
qcopy_file_sync_node(path_src, path_dst)
qcat_sync_node(path)
qcat_bytes_sync_node(path)
qload_datatable_sync_node(path, is_par=False)
```

---

## Examples

### Working with QFile

```python
import qlat_utils as q

# In-memory file with initial content
f = q.open_qfile_str("hello.txt", "r", content="Hello, world!")
print(f.content())  # "Hello, world!"

# Read from a QFile
f.seek_set(0)
line = f.getline()
print(line)  # "Hello, world!"

# Copy a QFile
f2 = f.copy()
print(f2.content())  # "Hello, world!"

# CRC-32 checksum
crc = f.compute_crc32()
print(crc)  # e.g., 222957957
```

### Creating and Reading a QAR Archive

```python
import qlat_utils as q
import os

# Create a temporary directory with files
os.makedirs("tmp_qar/inputs", exist_ok=True)
q.qtouch("tmp_qar/inputs/a.txt", "content a")
q.qtouch("tmp_qar/inputs/b.txt", "content b")

# Create QAR archive
q.qar_create("archive.qar", "tmp_qar")

# List contents
print(q.list_qar("archive.qar"))

# Open as QarFile for random access
qar = q.open_qar("archive.qar", "r")
print(qar.list())
f = qar.read("inputs/b.txt")
print(f.qcat())  # "content b"

# Read data directly as a string
data = qar.read_data("inputs/a.txt")
print(data)  # "content a"
qar.close()
```

### Basic File I/O

```python
import qlat_utils as q

# Write and read
q.qtouch("example.txt", "line1\nline2\n")
print(q.qcat("example.txt"))  # "line1\nline2\n"

# Append
q.qappend("example.txt", "line3\n")
print(q.qcat("example.txt"))  # "line1\nline2\nline3\n"

# CRC-32
crc = q.compute_crc32("example.txt")
print(crc)

# Load datatable
q.qtouch("table.txt", "1.0 2.0\n3.0 4.0\n")
table = q.qload_datatable("table.txt")
print(table)  # [[1.0, 2.0], [3.0, 4.0]]
```

### QAR with Index

```python
import qlat_utils as q

# Build an index for faster random access
q.qar_build_index("archive.qar")

# The index is saved to "archive.qar.idx"
# Subsequent reads use the index automatically
```
