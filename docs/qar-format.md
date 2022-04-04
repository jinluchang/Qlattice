# QAR format

Source code located in `qutils/qar.h`

File format is

```
FILE-FORMAT-LINE
[File1 name size in bytes stored as int64_t in little endian]
[File1 info size in bytes stored as int64_t in little endian]
[File1 data size in bytes stored as int64_t in little endian]
FILE1-NAME
FILE1-INFO
FILE1-DATA
[File2 name size in bytes stored as int64_t in little endian]
[File2 info size in bytes stored as int64_t in little endian]
[File2 data size in bytes stored as int64_t in little endian]
FILE2-NAME
FILE2-INFO
FILE2-DATA
...
```

### FILE-FORMAT-LINE

```
#!/usr/bin/env qar-glimpse
```

Note that there is no newline character between the data segments.


# FILE-INFO

Can store some metadata information about the file. Can be simply empty.
