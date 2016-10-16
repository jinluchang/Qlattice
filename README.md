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

2. `rng-state.h` https://github.com/waterret/RngState

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

## Random number generator description:

The state of the generator is effectively composed of the history of the
generator encoded as a string.

To generate random numbers, one computes the SHA-256 hash of the string.
The hash result is viewed as a `8` 32-bit unsigned integers.

The `8` 32-bit unsigned integers are merged into `4` 64-bit unsigned
integers. These `4` numbers are treated as the random numbers generated
by this random number generator.

Relavent source files: `distfiles/rng-state.h`, `qlat/field-rng.h`.

Relavent examples: `examples/rng-state-tests`, `examples/field-rng-tests`.

Also see the RNG's own repository

https://github.com/waterret/RngState
