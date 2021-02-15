import cqlat as c

from qlat.timer import *

from qlat.lat_io import *

@timer
def contract_pion_field(prop, tslice):
    ld = LatData()
    c.contract_pion_field(ld, prop, tslice)
    return ld
