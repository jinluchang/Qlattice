import qlat.c as c

from qlat.c import Geometry

def geo_reform(geo, multiplicity = 1, expansion_left = None, expansion_right = None):
    assert isinstance(geo, Geometry)
    if expansion_left is None:
        expansion_left = [ 0, 0, 0, 0, ]
    elif isinstance(expansion_left, int):
        e = expansion_left
        expansion_left = [ e, e, e, e, ]
    if expansion_right is None:
        expansion_right = [ 0, 0, 0, 0, ]
    elif isinstance(expansion_right, int):
        e = expansion_right
        expansion_right = [ e, e, e, e, ]
    geo_new = geo.copy()
    c.set_geo_reform(geo_new, multiplicity, expansion_left, expansion_right)
    return geo_new

def geo_eo(geo, eo = 0):
    assert isinstance(geo, Geometry)
    geo_new = Geometry([ 0, 0, 0, 0, ])
    geo_new = geo.copy()
    c.set_geo_eo(geo_new, eo)
    return geo_new
