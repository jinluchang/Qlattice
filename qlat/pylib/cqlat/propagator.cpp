#include "lib.h"

EXPORT(set_point_src_prop, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_geo = NULL;
  PyObject* p_xg = NULL;
  Complex value = 1.0;
  if (!PyArg_ParseTuple(args, "OOO|D", &p_prop, &p_geo, &p_xg, &value)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  Coordinate xg;
  py_convert(xg, p_xg);
  set_point_src(prop, geo, xg, value);
  Py_RETURN_NONE;
})

EXPORT(set_wall_src_prop, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_geo = NULL;
  int tslice = -1;
  PyObject* p_lmom = NULL;
  if (!PyArg_ParseTuple(args, "OOi|O", &p_prop, &p_geo, &tslice, &p_lmom)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  CoordinateD lmom;
  py_convert(lmom, p_lmom);
  set_wall_src(prop, geo, tslice, lmom);
  Py_RETURN_NONE;
})

EXPORT(set_rand_u1_src_psel, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_fu1 = NULL;
  PyObject* p_psel = NULL;
  PyObject* p_geo = NULL;
  PyObject* p_rs = NULL;
  if (!PyArg_ParseTuple(args, "OOOOO", &p_prop, &p_fu1, &p_psel, &p_geo, &p_rs)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  prop.init();
  FieldM<Complex, 1>& fu1 = py_convert_type_field<Complex, 1>(p_fu1);
  fu1.init();
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  const Geometry& geo = py_convert_type<Geometry>(p_geo);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  set_rand_u1_src_psel(prop, fu1, psel, geo, rs);
  Py_RETURN_NONE;
})

EXPORT(set_rand_u1_sol_psel, {
  using namespace qlat;
  PyObject* p_sp_prop = NULL;
  PyObject* p_prop = NULL;
  PyObject* p_fu1 = NULL;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sp_prop, &p_prop, &p_fu1, &p_psel)) {
    return NULL;
  }
  SelectedPoints<WilsonMatrix>& sp_prop =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_sp_prop);
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const FieldM<Complex, 1>& fu1 = py_convert_type_field<Complex, 1>(p_fu1);
  pqassert(fu1.geo().multiplicity == 1);
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  set_rand_u1_sol_psel(sp_prop, prop, fu1, psel);
  Py_RETURN_NONE;
})

EXPORT(set_rand_u1_src_fsel, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_fu1 = NULL;
  PyObject* p_fsel = NULL;
  PyObject* p_rs = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_prop, &p_fu1, &p_fsel, &p_rs)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  prop.init();
  FieldM<Complex, 1>& fu1 = py_convert_type_field<Complex, 1>(p_fu1);
  fu1.init();
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  set_rand_u1_src_fsel(prop, fu1, fsel, rs);
  Py_RETURN_NONE;
})

EXPORT(set_rand_u1_sol_fsel, {
  using namespace qlat;
  PyObject* p_sf_prop = NULL;
  PyObject* p_prop = NULL;
  PyObject* p_fu1 = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sf_prop, &p_prop, &p_fu1, &p_fsel)) {
    return NULL;
  }
  SelectedField<WilsonMatrix>& sf_prop =
      py_convert_type<SelectedField<WilsonMatrix> >(p_sf_prop);
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const FieldM<Complex, 1>& fu1 = py_convert_type_field<Complex, 1>(p_fu1);
  pqassert(fu1.geo().multiplicity == 1);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  set_rand_u1_sol_fsel(sf_prop, prop, fu1, fsel);
  Py_RETURN_NONE;
})

EXPORT(free_invert_prop, {
  using namespace qlat;
  PyObject* p_prop_sol = NULL;
  PyObject* p_prop_src = NULL;
  double mass = 0.0;
  double m5 = 1.0;
  PyObject* p_momtwist = NULL;
  if (!PyArg_ParseTuple(args, "OOd|dO", &p_prop_sol, &p_prop_src, &mass, &m5,
                        &p_momtwist)) {
    return NULL;
  }
  Propagator4d& prop_sol = py_convert_type<Propagator4d>(p_prop_sol);
  const Propagator4d& prop_src = py_convert_type<Propagator4d>(p_prop_src);
  CoordinateD momtwist;
  py_convert(momtwist, p_momtwist);
  free_invert(prop_sol, prop_src, mass, m5, momtwist);
  Py_RETURN_NONE;
})

EXPORT(convert_wm_from_mspincolor_prop, {
  using namespace qlat;
  PyObject* p_prop_wm = NULL;
  PyObject* p_prop_msc = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_wm, &p_prop_msc)) {
    return NULL;
  }
  Propagator4d& prop_wm = py_convert_type<Propagator4d>(p_prop_wm);
  const Propagator4d& prop_msc = py_convert_type<Propagator4d>(p_prop_msc);
  convert_wm_from_mspincolor(prop_wm, prop_msc);
  Py_RETURN_NONE;
})

EXPORT(convert_mspincolor_from_wm_prop, {
  using namespace qlat;
  PyObject* p_prop_msc = NULL;
  PyObject* p_prop_wm = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_msc, &p_prop_wm)) {
    return NULL;
  }
  Propagator4d& prop_msc = py_convert_type<Propagator4d>(p_prop_msc);
  const Propagator4d& prop_wm = py_convert_type<Propagator4d>(p_prop_wm);
  convert_mspincolor_from_wm(prop_msc, prop_wm);
  Py_RETURN_NONE;
})

EXPORT(convert_wm_from_mspincolor_sp_prop, {
  using namespace qlat;
  PyObject* p_prop_wm = NULL;
  PyObject* p_prop_msc = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_wm, &p_prop_msc)) {
    return NULL;
  }
  SelectedPoints<WilsonMatrix>& prop_wm =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop_wm);
  const SelectedPoints<WilsonMatrix>& prop_msc =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop_msc);
  convert_wm_from_mspincolor(prop_wm, prop_msc);
  Py_RETURN_NONE;
})

EXPORT(convert_mspincolor_from_wm_sp_prop, {
  using namespace qlat;
  PyObject* p_prop_msc = NULL;
  PyObject* p_prop_wm = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_msc, &p_prop_wm)) {
    return NULL;
  }
  SelectedPoints<WilsonMatrix>& prop_msc =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop_msc);
  const SelectedPoints<WilsonMatrix>& prop_wm =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop_wm);
  convert_mspincolor_from_wm(prop_msc, prop_wm);
  Py_RETURN_NONE;
})

EXPORT(convert_wm_from_mspincolor_s_prop, {
  using namespace qlat;
  PyObject* p_prop_wm = NULL;
  PyObject* p_prop_msc = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_wm, &p_prop_msc)) {
    return NULL;
  }
  SelectedField<WilsonMatrix>& prop_wm =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop_wm);
  const SelectedField<WilsonMatrix>& prop_msc =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop_msc);
  convert_wm_from_mspincolor(prop_wm, prop_msc);
  Py_RETURN_NONE;
})

EXPORT(convert_mspincolor_from_wm_s_prop, {
  using namespace qlat;
  PyObject* p_prop_msc = NULL;
  PyObject* p_prop_wm = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop_msc, &p_prop_wm)) {
    return NULL;
  }
  SelectedField<WilsonMatrix>& prop_msc =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop_msc);
  const SelectedField<WilsonMatrix>& prop_wm =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop_wm);
  convert_mspincolor_from_wm(prop_msc, prop_wm);
  Py_RETURN_NONE;
})

EXPORT(free_scalar_invert_mom_cfield, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double mass = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_field, &mass)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  pqassert(pf.ctype == "Complex");
  Field<Complex>& f = *(Field<Complex>*)pf.cdata;
  const Geometry& geo = f.geo();
  const Coordinate total_site = geo.total_site();
  const CoordinateD momtwist;
  const double m_pi_sq = 4.0 * sqr(std::sinh(mass / 2.0));
  qacc_for(index, geo.local_volume(), {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk, ks;
    double s2 = 0.0;
    for (int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      ks[i] = 2.0 * std::sin(kk[i] / 2.0);
      s2 += sqr(ks[i]);
    }
    const double fac = 1.0 / (m_pi_sq + s2);
    Vector<Complex> v = f.get_elems(kl);
    for (int i = 0; i < v.size(); ++i) {
      v[i] *= fac;
    }
  })
  Py_RETURN_NONE;
})

EXPORT(flip_tpbc_with_tslice_sp_prop, {
  using namespace qlat;
  PyObject* p_sp_prop = NULL;
  int tslice_flip_tpbc = -1;
  if (!PyArg_ParseTuple(args, "Oi", &p_sp_prop, &tslice_flip_tpbc)) {
    return NULL;
  }
  SelectedPoints<WilsonMatrix>& sp_prop =
      py_convert_type_spoints<WilsonMatrix>(p_sp_prop);
  const PointSelection& psel =
      py_convert_type<PointSelection>(p_sp_prop, "psel");
  const Geometry& geo = py_convert_type<Geometry>(p_sp_prop, "psel", "geo");
  const int t_size = geo.total_site()[3];
  flip_tpbc_with_tslice(sp_prop, psel, tslice_flip_tpbc, t_size);
  Py_RETURN_NONE;
})

EXPORT(flip_tpbc_with_tslice_s_prop, {
  using namespace qlat;
  PyObject* p_s_prop = NULL;
  int tslice_flip_tpbc = -1;
  if (!PyArg_ParseTuple(args, "Oi", &p_s_prop, &tslice_flip_tpbc)) {
    return NULL;
  }
  SelectedField<WilsonMatrix>& s_prop =
      py_convert_type_sfield<WilsonMatrix>(p_s_prop);
  const FieldSelection& fsel =
      py_convert_type<FieldSelection>(p_s_prop, "fsel");
  flip_tpbc_with_tslice(s_prop, fsel, tslice_flip_tpbc);
  Py_RETURN_NONE;
})
