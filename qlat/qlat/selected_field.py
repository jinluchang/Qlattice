from qlat_utils import *
import qlat.c as c
from qlat.c import SelectedField, SelectedFieldBase

from qlat.field_selection_utils import *
from qlat.utils_io import *

@timer
def set_selected_field(sf, f):
    # deprecated use @=
    displayln_info("set_selected_field: deprecated")
    assert isinstance(sf, SelectedFieldBase)
    if isinstance(f, FieldBase):
        c.set_sfield_field(sf, f)
    elif isinstance(f, SelectedFieldBase):
        c.set_sfield_sfield(sf, f)
    else:
        raise Exception("set_selected_field")
