# Qlattice

A simple lattice QCD library.

1. Fast Fourier transformation is implemented based on FFTW3.

2. A random number generator is implemented based on the SHA-256 hash.

## Installation:

The library itself is only composed of header files located in the
`qlat` directory.

The major dependencies are included in the `distfiles` directory.

Among the dependencies, two special libraries are provided.

1. `timer.h` https://github.com/waterret/Timer

2. `hash-cpp` http://create.stephan-brumme.com/hash-library

There are scripts used to install all the dependencies, which are located
in the `scripts` directory.

The `build.sh` will install everything into the `$prefix` directory
by running all the scripts in the `scripts` directory. For example,
the following command will install everything into `DEST_DIRECTORY`.

`$ prefix=DEST_DIRECTORY ./build.sh`

If `$prefix` is not set, the default value in `setenv.sh` is used,
which is `$HOME/qlat-builds/VERSION_NUMBER`.

## Usage:

A sample `Makefile` is provided which can compile and run a simple
program using the library. The `Makefile` assumes that the library
and all its dependencies are installed in their default locations.

There are also example programs provided in the examples directory. Once
the library is completed installed, one can run the following command
to compile all the examples:

`$ make -C examples run`

## Random number generator description (outdated):

The state of the generator is a 512-bit unsigned integer stored in little
endian format.

To generate random numbers, one computes the SHA-256 hash of the
state. The hash result is viewed as a 256-bit unsigned integer stored
in little endian format.

The 256-bit unsigned integer is split into `4` 64-bit unsigned
integers. These `4` numbers are treated as the random numbers generated
by this random number generator.

The state of the generator is then increased by `1`, and ready to produce
more random numbers.

The seeding strategy:

The random number generator state is initialized by four 64-bit unsigned
integers: `seed`, `type`, `traj`, `index`.

`(seed * (2^64)^3 + type * (2^64)^2 + traj * (2^64) + index) * 2^256`

For a lattice application, one can easily initialize a generator for each
lattice site by giving them different `index`. For each configuration,
one can use its `traj` to initialize a set of generators. One can even
initialize independent sets of generators for different purposes with the
`type` integer.

Relavent source files: `qlat/rng-state.h`, `qlat/field-rng.h`.

Relavent examples: `examples/rng-state-tests`, `examples/field-rng-tests`.
