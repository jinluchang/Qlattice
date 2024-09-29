# Installation

```bash
export prefix=DEST_DIRECTORY
```

If `prefix` environmental variable is not set (or empty), the default `DEST_DIRECTORY` is `"$HOME/qlat-build/$(cat VERSION)"`.

Optionally, you can also specify a temp directory for the installation procedure:

```bash
export temp_dir=TEMP_DIRECTORY
```

If `temp_dir` environmental variable is not set (or empty), the default `TEMP_DIRECTORY` is `$prefix/temp`. An efficient choice for `TEMP_DIRECTORY` can be `/dev/shm/$(whoami)/temp`.

## Install on Ubuntu

```bash
sudo apt-get install -y python3-full
sudo apt-get install -y python3-dev python3-venv python3-wheel
sudo apt-get install -y git bzip2 autoconf unzip
sudo apt-get install -y libopenmpi-dev ninja-build patchelf
sudo apt-get install -y libeigen3-dev libgsl-dev libopenblas-dev
sudo apt-get install -y zlib1g-dev libssl-dev libmpfr-dev
sudo apt-get install -y gnuplot texlive-metapost poppler-utils
sudo apt-get install -y libfftw3-dev
```

### Qlattice only

```bash
./build.sh
source DEST_DIRECTORY/setenv.sh
```

### Qlattice and Grid/GPT

```bash
./scripts/download-core.sh
./build.sh default-gpt
source DEST_DIRECTORY/setenv.sh
```

## Install on Mac

```bash
brew install llvm autoconf automake
brew install coreutils flock findutils pkg-config
brew install open-mpi ninja patchelf
brew install eigen gsl
brew install zlib openssl@3 mpfr
brew install gnuplot texlive poppler
brew install fftw
```

### Qlattice only

```bash
./build.sh
source DEST_DIRECTORY/setenv.sh
```

### Qlattice and Grid/GPT

```bash
./scripts/download-core.sh
./build.sh default-gpt
source DEST_DIRECTORY/setenv.sh
```

## Install on UCONN HPC

First, download dependencies downloaded to the `distfiles` directory.

```bash
./scripts/download.sh
```

Note you may need to download on your local computer and copy the data to the server.

Run the build script:

```bash
./build.sh uconn
source DEST_DIRECTORY/setenv.sh
```

## General instructions

The library itself is only composed of header files located in the `qlat` directory.

The major dependencies can be downloaded with

```bash
./scripts/download.sh
```

These dependencies will be downloaded to the `distfiles` directory.

There are scripts in the `scripts` directory used to install the dependencies.

To install, we first set the `prefix` environment variable to specify the target directory.

```bash
export prefix=DEST_DIRECTORY
```

If `$prefix` is not set, the default value in `conf.sh` is used, which is `$HOME/qlat-build/default`.

The `build.sh` will install everything into the `$prefix` directory by running all the scripts in the `scripts` directory. For example, the following command will install everything into `DEST_DIRECTORY`.

```bash
./build.sh
```

On some systems, one may need to specify the system type:

```bash
./build.sh TARGET
```

`TARGET` can be `default`, `default-gpt`, `uconn`, `qcdserver`, `bnlknl`, `summit`, ... This script will call `./scripts/build.TARGET.sh` to perform relevant builds.

The environment variables needed to use the library can be set with the following command:

```bash
source DEST_DIRECTORY/setenv.sh
```

There are few different scripts to build the `Grid` library. Choose one best suit the machine (or create a custom one).

After the first complete install, one can re-install individual components by running the specific script. For example, to just re-install the `Qlattice` header files and python library:

```bash
./scripts/qlat-utils.sh
./scripts/qlat.sh
```

It is also possible to build `CPS`, `Grid`, and `GPT`. For examples:

```bash
./scripts/qmp.sh
./scripts/qio.sh
./scripts/cps.sh
./scripts/qlat-cps.sh

./scripts/c-lime.sh
./scripts/grid-clehner.avx2.sh
./scripts/qlat-grid.sh

./scripts/gpt.sh
```

It may be convenient to create a symbolic link `$HOME/qlat-build/default`, which points to the actual directory `DEST_DIRECTORY`. The `prefix` environment variable can be empty if the symbolic link is created.

## Custom installation

If there is no existing `./scripts/build.TARGET.sh` for the machine. You may need to create a new one. Also, you may need to create a `./scripts/setenv.TARGET.sh` script.

The purpose of the `./scripts/setenv.TARGET.sh` is to create the `DEST_DIRECTORY/setenv.sh` scripts, which can be sourced to set the environmental variables. The `DEST_DIRECTORY/setenv.sh` scripts will be sourced automatically in all the other installation scripts (e.g. `./scripts/qlat.sh`).

The `./scripts/qcore.sh` script will create the following wrappers:

```
DEST_DIRECTORY/qcore/bin/CC.sh
DEST_DIRECTORY/qcore/bin/CXX.sh
DEST_DIRECTORY/qcore/bin/MPICC.sh
DEST_DIRECTORY/qcore/bin/MPICXX.sh
```

These wrappers will simply try to call the appropriate corresponding compiler.

The environment variables, `CC`, `CXX`, `MPICC`, `MPICXX`, affect most installation scripts in `./scripts`. If these environment variables are empty, then the default `setenv.sh` will use `qcore/setenv.sh` to set them to use the wrappers created by `./scripts/qcore.sh`. This behavior can be overridden by setting these environment variables to the desired non-empty value in `./scripts/setenv.TARGET.sh` or `./scripts/build.TARGET.sh`.

Once `python3` is available, the compiler wrappers will support the `--wrapper-remove-arg=` option. One can set `CFLAGS` and/or `CXXFLAGS` with this option to remove some unwanted flags generated by the configure script.

## Usage

The environment variables needed to use the library can be set with the following command:

```bash
source DEST_DIRECTORY/setenv.sh
```

There are example programs provided in the examples directory. Once the library is completed installed, one can run the following command to compile and run all the examples:

```bash
./scripts/qlat-examples-cpp.sh
./scripts/qlat-examples-py.sh
./scripts/qlat-examples-py-gpt.sh # require GPT
./scripts/qlat-examples-py-cps.sh # require CPS and qlat-cps
./scripts/qlat-examples-cpp-grid.sh # require Grid and qlat-grid
```
