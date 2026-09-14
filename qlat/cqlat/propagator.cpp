#include "lib.h"

EXPORT(set_rand_u1_src_psel, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_fu1 = NULL;
  PyObject* p_psel = NULL;
  PyObject* p_geo = NULL;
  PyObject* p_rs = NULL;
  if (!PyArg_ParseTuple(args, "OOOOO", &p_prop, &p_fu1, &p_psel, &p_geo,
                        &p_rs)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  prop.init();
  FieldM<ComplexD, 1>& fu1 = py_convert_type_field<ComplexD, 1>(p_fu1);
  fu1.init();
  const PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
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
      py_convert_type<SelectedPoints<WilsonMatrix>>(p_sp_prop);
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const FieldM<ComplexD, 1>& fu1 = py_convert_type_field<ComplexD, 1>(p_fu1);
  qassert(fu1.multiplicity == 1);
  const PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
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
  FieldM<ComplexD, 1>& fu1 = py_convert_type_field<ComplexD, 1>(p_fu1);
  fu1.init();
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  QLAT_DIAGNOSTIC_POP;
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
      py_convert_type<SelectedField<WilsonMatrix>>(p_sf_prop);
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const FieldM<ComplexD, 1>& fu1 = py_convert_type_field<ComplexD, 1>(p_fu1);
  qassert(fu1.multiplicity == 1);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  QLAT_DIAGNOSTIC_POP;
  set_rand_u1_sol_fsel(sf_prop, prop, fu1, fsel);
  Py_RETURN_NONE;
})

EXPORT(flip_tpbc_with_tslice_sp_prop, {
  using namespace qlat;
  PyObject* p_sp_prop = NULL;
  Int tslice_flip_tpbc = -1;
  if (!PyArg_ParseTuple(args, "Oi", &p_sp_prop, &tslice_flip_tpbc)) {
    return NULL;
  }
  SelectedPoints<WilsonMatrix>& sp_prop =
      py_convert_type_spoints<WilsonMatrix>(p_sp_prop);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const PointsSelection& psel =
      py_convert_type<PointsSelection>(p_sp_prop, "psel");
  const Geometry& geo = py_convert_type<Geometry>(p_sp_prop, "psel", "geo");
  QLAT_DIAGNOSTIC_POP;
  const Int t_size = geo.total_site()[3];
  flip_tpbc_with_tslice(sp_prop, psel, tslice_flip_tpbc, t_size);
  Py_RETURN_NONE;
})

EXPORT(flip_tpbc_with_tslice_s_prop, {
  using namespace qlat;
  PyObject* p_s_prop = NULL;
  Int tslice_flip_tpbc = -1;
  if (!PyArg_ParseTuple(args, "Oi", &p_s_prop, &tslice_flip_tpbc)) {
    return NULL;
  }
  SelectedField<WilsonMatrix>& s_prop =
      py_convert_type_sfield<WilsonMatrix>(p_s_prop);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel =
      py_convert_type<FieldSelection>(p_s_prop, "fsel");
  QLAT_DIAGNOSTIC_POP;
  flip_tpbc_with_tslice(s_prop, fsel, tslice_flip_tpbc);
  Py_RETURN_NONE;
})
