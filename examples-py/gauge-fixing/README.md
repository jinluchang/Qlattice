# Gauge fixing example project

A `py-qlat-ext`-style Python extension project (meson + Cython + C++) that
implements plain local steepest-descent **Landau** and **Coulomb** gauge
fixing on a distributed `qlat` gauge field.  The local update is

    g(x) = exp(alpha * D(x))

with `D(x)` the anti-hermitian traceless part of the staple-sum matrix; one
sweep is a Jacobi sweep over all sites.  There is no Fourier acceleration and
no FFT.

## Test

The example is registered in `examples-py/Makefile` as the
`gauge-fixing-ext.log` target and is run by the examples-py test suite (the
extension project cannot be run by `nixpkgs/run-one-example-py.py`, which only
copies a single top-level `.py` script):

    make -B gauge-fixing-ext.log

The test is `gauge-fixing-ext.py`; it runs with 2 MPI processes on a
`4^4` lattice, fixes Landau and Coulomb gauges, and checks convergence, the
gauge condition, unitarity and plaquette (gauge) invariance.

## Cross-check against `qlat-ref`

`cross_check.py` fixes the same deterministic random field with this C++
implementation and with the pure-NumPy reference `qlat-ref`, and compares the
converged functional, the RMS divergence, the unitarity error and the fixed
link matrices.  It is not part of the CI suite because the nix example build
only packages `examples-py/`, so run it manually from a checkout after
building the extension project:

    source result-py-local/bin/setenv-qlat.sh
    PYTHONPATH=<build-install>/lib/python3/dist-packages:$PYTHONPATH \
      mpiexec -n 2 --oversubscribe --bind-to none \
      python3 -m mpi4py examples-py/gauge-fixing/cross_check.py
