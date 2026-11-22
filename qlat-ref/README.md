# `qlat-ref`

Pure-NumPy reference implementations used to cross-check the C++/Cython code
in `qlat`.

The gauge-fixing reference is `qlat_ref.gauge`; its documentation is
[`docs/qlat-ref/qlat_gauge.md`](../docs/qlat-ref/qlat_gauge.md).

```python
import numpy as np
import qlat_ref as qr

U = qr.get_random_gauge_field((2, 2, 2, 2), seed=0)
U_fixed, info = qr.fix_landau_gauge(U)
print(info["functional"], info["divergence"])
```

The package is pure Python (only `numpy` is required), single process, single
thread, and vectorized over all lattice sites.  It has no build step: add the
`qlat-ref` directory to `PYTHONPATH` (or run from inside it) and import
`qlat_ref`.

Self-checks:

```bash
cd qlat-ref
python3 -m qlat_ref.tests
```

The checks use tiny lattices and finish in a couple of seconds.
