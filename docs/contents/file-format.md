# File format

## Field (single file)

Text header:

```c++
BEGIN_FIELD_HEADER
field_version = 1.0
total_site[0] = 32
total_site[1] = 32
total_site[2] = 32
total_site[3] = 64
multiplicity = 64
sizeof(M) = 16
field_crc32 = DF96CC3E
END_HEADER
```

Binary data in big endian (double or float depends on size):

```c++
data[t][z][y][x][index_of_elem_on_this_site]
```

Can read or write data in this format via:

```c++
template <class M>
long write_field(const Field<M>& f, const std::string& path,
                 const Coordinate& new_size_node = Coordinate());
```
```c++
template <class M>
long read_field(Field<M>& f, const std::string& path,
                const Coordinate& new_size_node_ = Coordinate());
```
```c++
template <class M>
long write_field_float_from_double(
    const Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate());
```
```c++
template <class M>
long write_field_double(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node = Coordinate());
```
```c++
template <class M>
long read_field_double_from_float(
    Field<M>& f, const std::string& path,
    const Coordinate& new_size_node = Coordinate());
```
```c++
template <class M>
long read_field_double(Field<M>& f, const std::string& path,
                       const Coordinate& new_size_node = Coordinate());
```

## Subset of field (single file)

Selection of the field is should be saved as a field:

```c++
BEGIN_FIELD_HEADER
field_version = 1.0
total_site[0] = 24
total_site[1] = 24
total_site[2] = 24
total_site[3] = 64
multiplicity = 1
sizeof(M) = 8
field_crc32 = FD837EDD
END_HEADER
```

The elements of the field is of type ``long``  in big endian. Value ``-1``  means not selected. Positive number means selected points and the number implies the order they are selected, starting from ``0``. Each time slice is independently treated.

Text header of the selected field file:

```c++
BEGIN_SELECTED_FIELD_HEADER
selected_field_version = 1.0
total_site[0] = 24
total_site[1] = 24
total_site[2] = 24
total_site[3] = 64
n_per_tslice = 864
multiplicity = 288
sizeof(M) = 4
selected_field_crc32 = 9F68E3F2
END_HEADER
```

Data is stored in the same order as the field format but only for the selected sites.

Example on how to read/write the sparse field is ``examples-cpp/selected-field/main.cpp``.
