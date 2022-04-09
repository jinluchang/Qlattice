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

NOTE: There are two newline characters in the end of the FILE-FORMAT-LINE.
Quote symbols are not part of the FILE-FORMAT-LINE.

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

Quote symbols are not part of the FILE-HEADER.

# FILE-INFO

Can store some metadata information about the file. Can be simply empty.

# Segment size

```
[FILE-HEADER size] + 1 + [FILE-NAME size] + 1 + [FILE-INFO size] + 1 + [FILE-DATA size] + 2
```
