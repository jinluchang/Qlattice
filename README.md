# Qlattice

A simple lattice QCD library.

## Installation:

The library itself is only composed of header files located in the `qlat` directory.

The major dependencies can be downloaded with

`$ ./scripts/download.sh`

These dependencies will be downloaded to the `distfiles` directory.

There are scripts in the `scripts` directory used to install the dependencies.

To install, we first set the `prefix` environment variable to specify the target directory.

`$ export prefix=DEST_DIRECTORY`

If `$prefix` is not set, the default value in `conf.sh` is used, which is `$HOME/qlat-build/default`.

The `build.sh` will install everything into the `$prefix` directory by running all the scripts in the `scripts` directory. For example, the following command will install everything into `DEST_DIRECTORY`.

`$ ./build.sh`

On some systems, one may need to specify the system type:

`$ ./build.sh TARGET`

`TARGET` can be `default`, `uconn`, `qcdserver`, `bnlknl`, `summit`, ... This script will call `./scripts/build.TARGET.sh` to perform relevant builds.

The environment variables needed to use the library can be set with the following command:

`$ source DEST_DIRECTORY/setenv.sh`

There are few different scripts to build the `Grid` library. Choose one best suit the machine (or create a custom one).

After the first complete install, one can re-install individual components by running the specific script. For example, to just re-install the `Qlattice` header files and python library:

`$ ./scripts/qlat.sh`

The most time consuming part of the above command is the compilation the `cqlat.so` dynamic library. To avoid the compilation of `cqlat.so`, but re-install the header files and the `qlat` python files:

`$ ./scripts/qlat-header.sh`

It is also possible to build `Grid` and `gpt`. For examples:

`$ ./scripts/grid.avx2.sh`

`$ ./scripts/gpt.sh`

It may be convenient to create a symbolic link `$HOME/qlat-build/default`, which points to the actual directory `DEST_DIRECTORY`. The `prefix` environment variable can be empty if the symbolic link is created.

### Install on UCONN HPC

First obtain an interactive session:

`$ fisbatch -n 2 -p generalsky`

The default build and installation directory is `$HOME/qlat-build/default`. If a different directory is needed:

`$ export prefix=DEST_DIRECTORY`

Run the build script:

`$ ./build.sh uconn`

End the interactive session.

`$ exit`

## Install on new machines

If there is no existing `./scripts/build.TARGET.sh` for the machine. You may need to create a new one. Also, you may need to create a `./scripts/setenv.TARGET.sh` script.

The purpose of the `./scripts/setenv.TARGET.sh` is to create the `DEST_DIRECTORY/setenv.sh` scripts, which can be sourced to set the environmental variables. The `DEST_DIRECTORY/setenv.sh` scripts will be sourced automatically in all the other installation scripts (e.g. `./scripts/qlat.sh`).

The `./scripts/compiler-wrappers.sh` script will create the following wrappers:

`DEST_DIRECTORY/bin/CC.sh`

`DEST_DIRECTORY/bin/CXX.sh`

`DEST_DIRECTORY/bin/MPICC.sh`

`DEST_DIRECTORY/bin/MPICXX.sh`

These wrappers will simply try to call the appropriate corresponding compiler.

The environment variables, `CC`, `CXX`, `MPICC`, `MPICXX`, affect most installation scripts in `./scripts`. If these environment variables are empty, then the default `setenv.sh` will set them to use the wrappers created by `./scripts/compiler-wrappers.sh`. This behavior can be overridden by setting these environment variables to the desired non-empty value in `./scripts/setenv.TARGET.sh` or `./scripts/build.TARGET.sh`.

Once `python3` is available, the compiler wrappers will support the `--wrapper-remove-arg=` option. One can set `CFLAGS` and/or `CXXFLAGS` with this option to remove some unwanted flags generated by the configure script.

## Usage:

A sample `Makefile` is provided which can compile and run a simple program using the library. The `Makefile` assumes that the library and all its dependencies are installed in their default locations.

There are also example programs provided in the examples directory. Once the library is completed installed, one can run the following command to compile all the examples:

`$ make -C examples run`

## Structure

```c++
const int DIMN = 4;
struct Coordinate : public array<int, DIMN> {};
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

The state of the generator is effectively composed of the history of the generator encoded as a string.

To generate random numbers, one computes the SHA-256 hash of the string.  The hash result is viewed as a `8` 32-bit unsigned integers.

The `8` 32-bit unsigned integers are merged into `4` 64-bit unsigned integers. These `4` numbers are treated as the random numbers generated by this random number generator.

Relevant source files: `qutils/rng-state.h`

### Environment variables

- `q_end_time`

  Program finish time in seconds since epoch.

  Used when `check_time_limit()`.

  Default is empty.

  `q_time_limit`

  Total running time of program in seconds.

  Used when `check_time_limit()`, but will be override if `q_end_time` is set.

  Default is `43200`.

- `q_budget`

  Default budget time in seconds.

  Used when `check_time_limit()`.

  Default is `900`.

- `q_field_init`

  Control how field objects' data are initialized.

  Choices are `fast` (default), `zero`, `random`.

- `q_mem_cache_max_size`

  Memory cache size in MB (per processes) for `qlat::vector` allocation.

  Default is `512 MB`.

- `q_num_threads`

  Number of OpenMP threads (will be override by `OMP_NUM_THREADS`).

  Default is `2`.

- `q_acc_num_threads`

  Number of `qacc` threads.

  Default is `32`.

- `q_verbose`

  Level of verbosity. Need to be more than `0` for the timing info to be shown automatically.

  Default is `0`.

- `q_timer_mini_auto_display`

  Minimum time between auto-display of timer information summary.

  Default is `5.0 * 60.0`.

- `q_timer_mini_auto_show`

  Minimum run time for a function for its information to be shown when it start or stop.

  Default is `1.0`.

- `q_timer_max_always_show`

  Maximum number of times to always show function start or stop.

  Default is `10`.

- `q_timer_max_func_name_len`

  Maximum length for a function name before truncation.

  Default is `50`.
