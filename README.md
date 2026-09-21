# Qlattice

A simple lattice QCD library.

Qlattice is a C++17 lattice QCD library with Cython bindings and a Python
interface. It is built around MPI: fields are decomposed over the lattice
geometry, and the high level operations (lattice I/O, propagator solving,
jackknife analysis, contractions, ...) are collective operations, so the same
script runs on a laptop and on a large cluster.

Contributions and pull requests are very welcome.

## Packages

| Package | Import | Contents |
|---|---|---|
| `qlat-utils` | `qlat_utils` | Core utilities: parallel-safe RNG (`RngState`), timers, coordinates, the QAR chunked file format, caches, `LatData`, data analysis (jackknife, fits, plotting via gnuplot) |
| `qlat` | `qlat` | Lattice QCD: geometry, gauge and fermion fields, HMC, propagators, contractions (`qlat/auto_contractor`), topology, Wilson flow, lattice I/O, job scripts |
| `qlat-grid` | `qlat_grid` | Interface to the [Grid](https://github.com/paboyle/Grid) framework |
| `qlat-cps` | `qlat_cps` | Interface to [CPS](https://github.com/RBC-UKQCD/cps) |
| (part of `qlat`) | `qlat_gpt` | Interface to [GPT](https://github.com/lehner/gpt) / Grid |

`qlat` re-exports everything from `qlat_utils`, so `import qlat as q` is
usually enough. `qlat-utils` can also be used on its own; it does not require
MPI.

## Documentation

<https://jinluchang.github.io/Qlattice>

The sources of the documentation are in [`docs/`](docs/): `docs/contents/`
holds the guides (installation, environment variables, command line
arguments, ...), `docs/qlat-utils/` and `docs/qlat/` document the Python
modules, and `jackknife.rst`, `auto-contractor.rst`, `qlat_gpt.md`, ... cover
the larger subsystems.

## Quick start with Nix (recommended)

Nix builds the whole dependency chain reproducibly (including OpenMPI, and
Grid, GPT and CPS for the corresponding interfaces):

```bash
./nixpkgs/install-py-local-kernel-with-nix.py    # creates ./result-py-local and a Jupyter kernel
source result-py-local/bin/setenv-qlat.sh        # PATH, PYTHONPATH, LD_LIBRARY_PATH, ...
```

Run a script with MPI (from a scratch directory, e.g. `tmp/`):

```bash
mkdir -p tmp && cd tmp
mpiexec -n 2 --oversubscribe --bind-to none python3 -m mpi4py ./my-script.py
```

```python
import qlat as q

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
gf.set_rand(q.RngState("gf").split("set_rand"), sigma=0.3, n_step=2)
gf.show_info()

ga = q.GaugeAction(2.13, -0.331)
rs = q.RngState("hmc")
for traj in range(10):
    q.run_hmc_pure_gauge(gf, ga, traj, rs.split("run_hmc_pure_gauge"))
    gf.show_info()

q.timer_display()
q.end_with_mpi()
```

Build variants (`--variant cuda`, `--cudasupport`, `--cu`, `--clang`,
`--pypi`, or `--all-variants`) and all options are listed by
`./nixpkgs/install-py-local-kernel-with-nix.py --help`.

## Installation without Nix

`qlat-utils` and `qlat` are standard Python packages built with
[meson-python](https://github.com/mesonbuild/meson-python) and are published
to PyPI (see [`release.py`](release.py)):

```bash
pip install qlat-utils -Ccompile-args="-j2"
pip install qlat -Ccompile-args="-j2"
```

`qlat-utils` only needs a C++ compiler, Python, NumPy and Cython; `qlat`
additionally needs MPI (`mpi4py`). The `qlat-grid`, `qlat-cps` and `qlat_gpt`
interfaces need Grid, CPS and GPT respectively.

The repository also keeps the original, non-Nix install scripts
(`./build.sh`, `./scripts/*.sh`, `./qcore/` compiler wrappers). See
[`docs/contents/install.md`](docs/contents/install.md) for the system
dependencies and for the per-machine installation instructions.

## Tests and examples

The `examples-*` directories are the canonical, tested reference for current
API usage — keep them up to date when the API changes:

| Directory | Contents |
|---|---|
| `examples-py/` | Python examples and tests (run with `mpiexec`) |
| `examples-py-gpt/` | Python examples using the GPT/Grid interface |
| `examples-py-cps/` | Python examples using the CPS interface |
| `examples-cpp/` | C++ examples and tests |
| `examples-cpp-grid/` | C++ examples using the Grid interface |

Each example prints `CHECK:` lines that are compared against a committed
reference log, so the tests double as regression tests. Run a single test with
the Nix environment, or the whole suite:

```bash
./nixpkgs/run-one-example-py.py utils        # one Python test
./nixpkgs/build-many-qlat-pkgs.py --group tests   # every test
```

## Contributing

- Keep the `examples-*` directories and the `docs/` files in sync with the
  code; each module docstring names the documentation file that belongs to it.
- Follow the style guides in [`AGENTS.md`](AGENTS.md) (C++, Python, Cython,
  formatting, lint, testing and commit conventions).

## License

GPL-3.0-or-later — see [`LICENSE`](LICENSE).
