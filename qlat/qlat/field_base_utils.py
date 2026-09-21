"""
Module ``qlat.field_base_utils``
=================================\n
Pure-Python field helpers that do not call C++ functions directly: the
``Field``/``SelectedField``/``SelectedPoints`` factories, index validation,
and the split/merge/accumulate helpers that delegate to the ``_cc_*``
methods of the Cython field classes.\n
"""

class q:
    from qlat_utils import (
        Coordinate,
        CoordinateD,
    )
    from .field_base import (
        FieldBase,
        SelectedFieldBase,
        SelectedPointsBase,
    )
    from .field_type_dict import (
        field_type_dict,
        selected_field_type_dict,
        selected_points_type_dict,
    )

### -------------------------------------------------------------------

def Field(ctype, geo=None, multiplicity=0):
    assert ctype in q.field_type_dict
    FieldType = q.field_type_dict[ctype]
    field = FieldType(geo, multiplicity)
    return field

def SelectedField(ctype, fsel, multiplicity=0):
    """
    SelectedField(ctype, fsel) with the default multiplicity == 0 creates an
    *empty*, uninitialized field that keeps fsel: it is meant to be filled
    later, e.g. with load_double / float_from_double (see examples-py
    selected-convert-io.py).  Pass a positive multiplicity to allocate now.
    """
    assert ctype in q.field_type_dict
    FieldType = q.selected_field_type_dict[ctype]
    field = FieldType(fsel, multiplicity)
    return field

def SelectedPoints(ctype, psel, multiplicity=0):
    """
    SelectedPoints(ctype, psel) with the default multiplicity == 0 creates an
    *empty*, uninitialized field that keeps psel; pass a positive multiplicity
    to allocate now.
    """
    assert ctype in q.field_type_dict
    FieldType = q.selected_points_type_dict[ctype]
    field = FieldType(psel, multiplicity)
    return field

def field_check_key(idx):
    """
    Validate a field index and return the NumPy index to use.\n
    Field buffers are ``(local_volume, multiplicity, *elem_shape)``,
    C-contiguous, indexed by the **flat local site index** with the first
    coordinate varying fastest, matching ``geo.coordinate_from_index``.\n
    A ``Coordinate`` (or a tuple/list containing one) is a common mistake
    because it looks like the C++ ``get_elem`` API; reject it with an
    actionable message instead of silently doing the wrong thing.
    """
    if isinstance(idx, (q.Coordinate, q.CoordinateD)):
        raise TypeError(
            f"field indices are flat local site indices, not {type(idx).__name__}"
            "; use get_elem_xg(xg, m) for global coordinates, or "
            "geo.index_from_coordinate(xl) to convert a local coordinate to a "
            "flat local index"
        )
    if isinstance(idx, (tuple, list)):
        for key in idx:
            if isinstance(key, (q.Coordinate, q.CoordinateD)):
                raise TypeError(
                    "field indices are flat local site indices; a tuple "
                    "containing a Coordinate would be interpreted by NumPy as a "
                    "fancy index over the site axis. Use "
                    "get_elem_xg(xg, m) for global coordinates, or "
                    "geo.index_from_coordinate(xl)"
                )
    return idx

### -------------------------------------------------------------------

def split_fields(fs, f):
    nf = len(fs)
    assert nf >= 1
    ctype = f.ctype
    for i in range(nf):
        if not isinstance(fs[i], q.FieldBase):
            fs[i] = Field(ctype)
        else:
            assert fs[i].ctype is ctype
    f._cc_split_fields(fs)

def merge_fields(f, fs):
    nf = len(fs)
    assert nf >= 1
    assert isinstance(f, q.FieldBase)
    assert f.ctype is fs[0].ctype
    f._cc_merge_fields(fs)

def merge_fields_ms(f, fms):
    """
    fms = [ (f0, m0,), (f1, m1,), ... ]
    f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    """
    multiplicity = len(fms)
    assert multiplicity >= 1
    assert isinstance(f, q.FieldBase)
    assert f.ctype is fms[0][0].ctype
    fs, ms = zip(*fms)
    f._cc_merge_fields_ms(fs, ms)

def mk_merged_fields_ms(fms):
    """
    fms = [ (f0, m0,), (f1, m1,), ... ]
    f.get_elem(x, m) = fms[m][0].get_elem(x, fms[m][1])
    return f
    """
    multiplicity = len(fms)
    assert multiplicity >= 1
    for m in range(multiplicity):
        assert isinstance(fms[m][0], q.FieldBase)
        assert isinstance(fms[m][1], int)
    ctype = fms[0][0].ctype
    for m in range(multiplicity):
        assert ctype is fms[m][0].ctype
    f = Field(ctype)
    merge_fields_ms(f, fms)
    return f

### -------------------------------------------------------------------
### low-level cqlat-compatible entry points

def get_mview_field(field):
    """
    Return a flat, writable memoryview of the field data.
    """
    assert isinstance(field, q.FieldBase)
    return field.mview()

def set_add_sfield(f_new, f):
    """
    ``f_new += f`` for two SelectedField objects with the same FieldSelection.
    """
    assert isinstance(f_new, q.SelectedFieldBase)
    assert isinstance(f, q.SelectedFieldBase)
    f_new._cc_iadd(f)

def set_mul_double_sfield(f, factor):
    """
    ``f *= factor`` for a SelectedField.
    """
    assert isinstance(f, q.SelectedFieldBase)
    f._cc_imul_double(float(factor))

def acc_field_sfield(f, f1):
    """
    Accumulate a SelectedField into a Field: ``f += f1``.
    """
    assert isinstance(f, q.FieldBase)
    assert isinstance(f1, q.SelectedFieldBase)
    assert f1.ctype is f.ctype
    f._cc_acc_field_sfield(f1, f1.fsel)

def acc_field_spfield(f, f1, geo=None, psel=None):
    """
    Accumulate a SelectedPoints into a Field: ``f += f1``.
    """
    assert isinstance(f, q.FieldBase)
    assert isinstance(f1, q.SelectedPointsBase)
    assert f1.ctype is f.ctype
    if psel is None:
        psel = f1.psel
    if geo is None:
        geo = psel.geo
    f._cc_acc_field_spfield(f1, geo, psel)

def glb_sum_tslice_long_sfield(sp, f, t_dir=3):
    """
    Global-sum a SelectedField over the spatial sites of each time slice
    into the SelectedPoints ``sp``.
    """
    assert isinstance(f, q.SelectedFieldBase)
    f._cc_glb_sum_tslice(sp, f.fsel, t_dir)
