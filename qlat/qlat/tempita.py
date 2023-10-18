import sys
from pathlib import Path
from Cython.Tempita import sub

template = Path(sys.argv[1]).read_text('utf8')
output = sub(template)
Path(sys.argv[2]).write_text(output, 'utf8')
