# QAR format

## Introduction

A simple file format which allow combining many files together into a single archive.

## Utilities

```
Usage: qar list path.qar
Usage: qar build-idx path.qar
Usage: qar create path.qar path
Usage: qar extract path.qar path
Usage: qar cp path_src path_dst
Usage: qar cat path1 path2 ...
Usage: qar l path1.qar path2.qar ...
Usage: qar b path1.qar path2.qar ...
Usage: qar lr path1.qar path2.qar ...
Usage: qar c path1 path2 ...
Usage: qar x path1.qar path2.qar ...
Usage: qar cr path1 path2 ...
       Remove folders after qar files created
Usage: qar xr path1.qar path2.qar ...
       Remove qar files after folder extracted
```

```
Usage: qar-glimpse path1.qar path2.qar ... path1 path2 ...
```

## Multiple volumes

Can store the single folder in multiple `qar` archives. For example:

```
folder.qar
folder.qar.v1
folder.qar.v2
folder.qar.v3
```

Each file is a valid `qar` file on its own, but only contains a portion of the content of the original folder.

Environment variable `q_qar_multi_vol_max_size` controls the size of the each `qar` volume file.

## Note

1. Only store file names and content of files (allow record some info about the file in FILE-INFO).
2. '/' is allowed in `FILE-NAME`, and represent directory structure after extraction.
3. Does not store directories in QAR format. During extraction, directories will only be created to allow the files to be extracted to the specified path.
4. Empty directories will not be recorded and recovered after extraction.

## Description

Source code located in `qlat-utils/include/qlat-utils/qar.h`

Example file is located in `docs/contents/qar-sample.qar`.

File format is

```
FILE-FORMAT-LINE

FILE-HEADER
[newline character]
FILE-NAME
[newline character]
FILE-INFO
[newline character]
FILE-DATA
[newline character]
[newline character]

FILE-HEADER
[newline character]
FILE-NAME
[newline character]
FILE-INFO
[newline character]
FILE-DATA
[newline character]
[newline character]

...

```

### FILE-FORMAT-LINE

```
"#!/usr/bin/env qar-glimpse\n\n"
```

NOTE: There are two newline characters in the end of the `FILE-FORMAT-LINE`.
Quote symbols are not part of the `FILE-FORMAT-LINE`.

### FILE-HEADER

```
"QAR-FILE"
[one or more space characters]
[FILE-NAME size in bytes stored in ASCII]
[one or more space characters]
[FILE-INFO size in bytes stored in ASCII]
[one or more space characters]
[FILE-DATA size in bytes stored in ASCII]
```

Quote symbols are not part of the `FILE-HEADER`.

### FILE-INFO

Can store some metadata information about the file. The default is simply empty.

### Segment size

```
[FILE-HEADER size] + 1 + [FILE-NAME size] + 1 + [FILE-INFO size] + 1 + [FILE-DATA size] + 2
```

### Example

```
#!/usr/bin/env qar-glimpse

QAR-FILE 13 0 20
filename1.txt

Contents for file1.


QAR-FILE 13 0 20
filename2.txt

Contents for file2.


QAR-FILE 13 0 20
filename3.txt

Contents for file3.


QAR-FILE 18 0 21
folder1/file-a.txt

Contents for file-a.


QAR-FILE 18 0 21
folder2/file-b.txt

Contents for file-b.


QAR-FILE 18 0 21
folder2/file-c.txt

Contents for file-c.


```
