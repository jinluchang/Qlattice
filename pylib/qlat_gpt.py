import gpt as g
import qlat as q

def get_gpt_default_size_node():
    l_size = 8 * 3 * 5
    t_size = l_size * 2
    grid = g.grid([l_size, l_size, l_size, t_size], g.double)
    return grid.mpi

def begin_with_gpt():
    size_node = get_gpt_default_size_node()
    id_node = g.rank()
    q.begin(id_node, size_node)

def mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag):
    return ctype + "," + str(list(total_site)) + "," + str(multiplicity) + "," + tag

def mk_grid(geo):
    return g.grid(list(geo.total_site()), g.double)

def mk_gpt_field(ctype, geo):
    if ctype == "ColorMatrix":
        return g.mcolor(mk_grid(geo))
    elif ctype == "WilsonMatrix":
        return g.mspincolor(mk_grid(geo))
    else:
        raise Exception("make_gpt_field")

@q.timer_verbose
def mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    geo = q.Geometry(total_site, multiplicity)
    f_gpt = mk_gpt_field(ctype, geo)
    f_qlat = q.Field(ctype, geo)
    lexicographic_coordinates = g.coordinates(f_gpt)
    buffer = f_qlat.mview()
    if tag == "qlat_from_gpt":
        qlat_from_gpt = g.copy_plan(buffer, f_gpt)
        qlat_from_gpt.destination += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buffer, 0, buffer.nbytes]])
        qlat_from_gpt.source += f_gpt.view[lexicographic_coordinates]
        qlat_from_gpt = qlat_from_gpt(local_only = True)
        return qlat_from_gpt
    elif tag == "gpt_from_qlat":
        gpt_from_qlat = g.copy_plan(f_gpt, buffer)
        gpt_from_qlat.source += g.global_memory_view(
            f_gpt.grid,
            [[f_gpt.grid.processor, buffer, 0, buffer.nbytes]])
        gpt_from_qlat.destination += f_gpt.view[lexicographic_coordinates]
        gpt_from_qlat = gpt_from_qlat(local_only = True)
        return gpt_from_qlat
    else:
        print(tag)
        raise Exception("mk_qlat_gpt_copy_plan")

qlat_gpt_copy_plan_cache = {}

def get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag):
    key = mk_qlat_gpt_copy_plan_key(ctype, total_site, multiplicity, tag)
    if key in qlat_gpt_copy_plan_cache:
        return qlat_gpt_copy_plan_cache[key]
    else:
        plan = mk_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
        qlat_gpt_copy_plan_cache[key] = plan
        return plan

@q.timer
def qlat_from_gpt_gauge_field(gpt_gf):
    assert len(gpt_gf) == 4
    ctype = "ColorMatrix"
    total_site = gpt_gf[0].grid.fdimensions
    multiplicity = 1
    tag = "qlat_from_gpt"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    geo = q.Geometry(total_site, 1)
    fs = [ q.Field(ctype, geo) for i in range(4)]
    assert len(fs) == 4
    for i in range(4):
        plan(fs[i].mview(), gpt_gf[i])
    gf = q.GaugeField()
    q.merge_fields(gf, fs)
    return gf

@q.timer
def gpt_from_qlat_gauge_field(gf):
    assert isinstance(gf, q.GaugeField)
    geo = gf.geo()
    ctype = "ColorMatrix"
    total_site = geo.total_site()
    multiplicity = 1
    tag = "gpt_from_qlat"
    plan = get_qlat_gpt_copy_plan(ctype, total_site, multiplicity, tag)
    fs = [None] * 4
    q.split_fields(fs, gf)
    assert len(fs) == 4
    grid = mk_grid(geo)
    gpt_gf = [g.mcolor(grid) for i in range(4)]
    for i in range(4):
        plan(gpt_gf[i], fs[i].mview())
    return gpt_gf

