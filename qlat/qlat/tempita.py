import sys
from pathlib import Path
from Cython.Tempita import sub

"""
Link to some docs
https://github.com/cython/cython/blob/master/Cython/Tempita/_tempita.py
As of writing the version is
https://github.com/cython/cython/blob/023d4af351042a3b6241dbe06cbb003b3ce1fb58/Cython/Tempita/_tempita.py
"""

template = Path(sys.argv[1]).read_text('utf8')
output = sub(template)
Path(sys.argv[2]).write_text(output, 'utf8')
