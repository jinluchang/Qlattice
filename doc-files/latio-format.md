# LAT-IO format

Source code located in `qutils/lat-io.h`

Example code is `examples/latio/main.C`

File format is

```
FILE-HEADER
BINARY-DATA
```

### FILE-HEADER

Example:

```
#!/usr/bin/env lat-io-glimpse
data_size
128
ndim: 3
"tsep"[4]: "0" "1" "2" "3"
"op"[2]: "0" "1"
"re-im"[2]: "re" "im"
crc32: 77A655DB
END_HEADER

```

### BINARY-DATA

Consist of the data stored in double precision (little endian) in sequential order.
