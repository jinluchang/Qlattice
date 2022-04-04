# QAR format

Source code located in `qutils/qar.h`

File format is

```
FILE-FORMAT-LINE

FILE-HEADER
[newline charactor]
FILE-NAME
[newline charactor]
FILE-INFO
[newline charactor]
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
Quote is not part of the FILE-FORMAT-LINE.

### FILE-HEADER

```
[FILE-NAME size in bytes stored in ASCII]
[space character]
[FILE-INFO size in bytes stored in ASCII]
[space character]
[FILE-DATA size in bytes stored in ASCII]
```

# FILE-INFO

Can store some metadata information about the file. Can be simply empty.

# Segment size

```
[FILE-HEADER size] + [FILE-NAME size] + [File info size] + [File info size] + 5
```
