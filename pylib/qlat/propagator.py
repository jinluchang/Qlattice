import cqlat as c

from qlat.field import *
from qlat.rng_state import *

from qlat.selected_field import *
from qlat.selected_points import *

class Propagator4d(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "WilsonMatrix", geo, 1)

class Prop(Propagator4d):

    pass

class SelProp(SelectedField):

    def __init__(self, fsel):
        SelectedField.__init__(self, "WilsonMatrix", fsel, 1)

class PselProp(SelectedPoints):

    def __init__(self, psel):
        SelectedPoints.__init__(self, "WilsonMatrix", psel, 1)
