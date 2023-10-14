from qlat_utils import *
import qlat.c as c

from qlat.geometry import *
from qlat.field_selection import *
from qlat.field import *

from qlat.c import SelectedPoints, SelectedPointsBase

@timer
def set_selected_points(sp, f):
    # deprecated use @=
    from qlat.selected_field import SelectedField
    assert isinstance(sp, SelectedPointsBase)
    if isinstance(f, FieldBase):
        c.set_spfield_field(sp, f)
    elif isinstance(f, SelectedFieldBase):
        c.set_spfield_sfield(sp, f)
    else:
        raise Exception("set_selected_points")
