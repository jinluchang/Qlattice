"""
Module ``qlat.tempita``
=======================

Tempita template substitution utility for code generation.

Documentation: ``docs/qlat/qlat_tempita.md``

.. note:: Update the documentation when updating this source file.

References:
  - https://github.com/cython/cython/blob/master/Cython/Tempita/_tempita.py
  - https://github.com/cython/cython/blob/023d4af351042a3b6241dbe06cbb003b3ce1fb58/Cython/Tempita/_tempita.py
"""

import sys
from pathlib import Path
from Cython.Tempita import sub

template = Path(sys.argv[1]).read_text("utf8")
output = sub(template)
Path(sys.argv[2]).write_text(output, "utf8")
