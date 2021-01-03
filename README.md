# Qlattice

A simple lattice QCD library.

## Installation:

The library itself is only composed of header files located in the
`qlat` directory.

The major dependencies are included in the `distfiles` directory.

There are scripts used to install all the dependencies, which are located
in the `scripts` directory.

The `build.sh` will install everything into the `$prefix` directory
by running all the scripts in the `scripts` directory. For example,
the following command will install everything into `DEST_DIRECTORY`.

`$ prefix=DEST_DIRECTORY ./build.sh`

If `$prefix` is not set, the default value in `setenv.sh` is used,
which is `$HOME/qlat-build/VERSION_NUMBER`.

## Usage:

A sample `Makefile` is provided which can compile and run a simple
program using the library. The `Makefile` assumes that the library
and all its dependencies are installed in their default locations.

There are also example programs provided in the examples directory. Once
the library is completed installed, one can run the following command
to compile all the examples:

`$ make -C examples run`

## Structure

```c++
const int DIMN = 4;
struct Coordinate : public std::array<int, DIMN> {};
```

### GeometryNode

```c++
struct GeometryNode {
  bool initialized;
  int num_node;
  int id_node;
  Coordinate size_node;
  Coordinate coor_node;
};
```

### Geometry

```c++
struct Geometry {
  bool initialized;
  GeometryNode geon;
  int eo;  // 0:full; 1:odd; 2:even
  int multiplicity;
  Coordinate node_site;
  Coordinate expansion_left;
  Coordinate expansion_right;
  Coordinate node_site_expanded;
};
```

### Field

```c++
template <class M>
struct Field {
  bool initialized;
  box<Geometry> geo;
  vector<M> field;
};
```

### FieldM

```c++
template <class M, int multiplicity>
struct FieldM : Field<M> {};
```

## Random number generator description:

The state of the generator is effectively composed of the history of the
generator encoded as a string.

To generate random numbers, one computes the SHA-256 hash of the string.
The hash result is viewed as a `8` 32-bit unsigned integers.

The `8` 32-bit unsigned integers are merged into `4` 64-bit unsigned
integers. These `4` numbers are treated as the random numbers generated
by this random number generator.

Relevant source files: `qutils/rng-state.h`.

Relevant examples: `examples/rng-state-tests`, `examples/field-rng-tests`.

### Environment variables

- ``q_end_time``

  Program finish time in seconds since epoch.

  Used when ``check_time_limit()``.

  Default is empty.

  ``q_time_limit``

  Total running time of program in seconds.

  Used when ``check_time_limit()``.

  Default is empty.

- ``q_budget``

  Default budget time in seconds.

  Used when ``check_time_limit()``.

  Default is ``15*60``.

- ``q_field_init``

  Control how field objects' data are initialized.

  Choices are ``fast`` (default), ``zero``, ``random``.

- ``q_mem_cache_max_size``

  Memory cache size in MB (per processes) for ``qlat::vector`` allocation.

  Default is ``512 MB``.

- ``q_acc_num_threads``

  Number of ``qacc`` threads.

  Default is ``32``.
