# QAR format

Source code located in `qutils/qar.h`

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

# FILE-INFO

Can store some metadata information about the file. The default is simply empty.

# Segment size

```
[FILE-HEADER size] + 1 + [FILE-NAME size] + 1 + [FILE-INFO size] + 1 + [FILE-DATA size] + 2
```

# Note

1. Only store file names and content of files (allow record some info about the file in FILE-INFO).
2. '/' is allowed in `FILE-NAME`, and represent directory structure after extraction.
3. Does not store directories in QAR format. During extraction, directories will only be created to allow the files to be extracted to the specified path.
4. Empty directories will not be recorded and recovered after extraction.
