#pragma once

#include <qlat/core.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <ctime>
#include <fstream>
#include <vector>

namespace qlat
{  //

inline CoordinateD lattice_mom_mult(const Coordinate& total_site)
{
  return 2 * PI / CoordinateD(total_site);
}

inline CoordinateD lattice_mom_mult(const Geometry& geo)
{
  return lattice_mom_mult(geo.total_site());
}

void set_mom_phase_field(FieldM<ComplexD, 1>& f, const CoordinateD& mom);

void set_phase_field(FieldM<ComplexD, 1>& f, const CoordinateD& lmom);

// --------------------

template <class M>
bool is_initialized(const Field<M>& f)
{
  return f.initialized;
}

template <class M>
const Field<M>& operator+=(Field<M>& f, const Field<M>& f1)
{
  TIMER("field_operator+=");
  if (not is_initialized(f)) {
    f = f1;
    return f;
  }
  Qassert(is_matching_geo(f.geo(), f1.geo()));
  Qassert(f.multiplicity == f1.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      v[m] += v1[m];
    }
  });
  return f;
}

template <class M>
const Field<M>& operator-=(Field<M>& f, const Field<M>& f1)
{
  TIMER("field_operator-=");
  if (not is_initialized(f)) {
    f.init(f1.geo(), f1.multiplicity);
    set_zero(f);
    f -= f1;
    return f;
  }
  Qassert(is_matching_geo(f.geo(), f1.geo()));
  if (f.multiplicity != f1.multiplicity) {
    qerr(fname +
         ssprintf(": f.mult=%d ; f1.mult=%d", f.multiplicity, f1.multiplicity));
  }
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      v[m] -= v1[m];
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Field<RealD>& f_factor)
{
  TIMER("field_operator*=(F,FD)");
  Qassert(is_matching_geo(f.geo(), f_factor.geo()));
  const Int multiplicity_f = f_factor.multiplicity;
  Qassert(multiplicity_f == 1 or multiplicity_f == f.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const RealD fac = f_factor.get_elem(xl);
      for (Int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == f.multiplicity);
      Vector<RealD> fac = f_factor.get_elems_const(xl);
      for (Int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac[m];
      }
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Field<ComplexD>& f_factor)
{
  TIMER("field_operator*=(F,FC)");
  Qassert(is_matching_geo(f.geo(), f_factor.geo()));
  const Int multiplicity_f = f_factor.multiplicity;
  Qassert(multiplicity_f == 1 or multiplicity_f == f.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const ComplexD fac = f_factor.get_elem(xl);
      for (Int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == f.multiplicity);
      const Vector<ComplexD> fac = f_factor.get_elems_const(xl);
      for (Int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac[m];
      }
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const RealD factor)
{
  TIMER("field_operator*=(F,D)");
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const ComplexD& factor)
{
  TIMER("field_operator*=(F,C)");
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
RealD qnorm(const Field<M>& f)
{
  TIMER("qnorm(f)");
  const Geometry& geo = f.geo();
  RealD sum = 0.0;
#pragma omp parallel
  {
    RealD psum = 0.0;
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate x = geo.coordinate_from_index(index);
      const Vector<M> fx = f.get_elems_const(x);
      for (Int m = 0; m < f.multiplicity; ++m) {
        psum += qnorm(fx[m]);
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        sum += psum;
      }
    }
  }
  glb_sum(sum);
  return sum;
}

template <class M>
void qnorm_field(Field<RealD>& f, const Field<M>& f1)
{
  TIMER("qnorm_field(f,f1)");
  f.init(f1.geo(), 1);
  qacc_for(index, f1.geo().local_volume(), {
    // const Geometry& geo = f1.geo();
    const Vector<M> f1v = f1.get_elems_const(index);
    f.get_elem(index) = qnorm(f1v);
  });
}

template <class M>
RealD qnorm_double(const Field<M>& f1, const Field<M>& f2)
{
  const Geometry& geo = f1.geo();
  Qassert(geo.is_only_local);
  Qassert(geo == f2.geo());
  RealD sum = qnorm_double(get_data(f1), get_data(f2));
  glb_sum(sum);
  return sum;
}

template <class M>
qacc Long get_data_size(const Field<M>& f)
// NOT including the expended parts, only local volume data size
// only size on one node
{
  return f.geo().local_volume() * f.multiplicity * sizeof(M);
}

void set_sqrt_field(Field<RealD>& f, const Field<RealD>& f1);

// --------------------

template <class M, class N>
void assign(Field<N>& f, const Field<M>& f1)
{
  const Geometry& geo = f1.geo();
  Qassert(geo.is_only_local);
  f.init();
  f.init(geo, f1.multiplicity * sizeof(M) / sizeof(N));
  Qassert(f.multiplicity * sizeof(N) == f1.multiplicity * sizeof(M));
  qacc_for(index, geo.local_volume(), {
    const Vector<M> v1 = f1.get_elems_const(index);
    Vector<N> v = f.get_elems(index);
    assign(v, v1);
  });
}

template <class M>
std::vector<M> field_sum(const Field<M>& f)
// length = multiplicity
{
  TIMER("field_sum");
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec(multiplicity);
  set_zero(vec);
#pragma omp parallel
  {
    std::vector<M> pvec(multiplicity);
    set_zero(pvec);
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<M> fvec = f.get_elems_const(xl);
      for (Int m = 0; m < multiplicity; ++m) {
        pvec[m] += fvec[m];
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        for (Int m = 0; m < multiplicity; ++m) {
          vec[m] += pvec[m];
        }
      }
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_sum_tslice(const Field<M>& f, const Int t_dir = 3)
// length = t_size * multiplicity
{
  TIMER("field_sum_tslice");
  const Geometry& geo = f.geo();
  const Int t_size = geo.total_site()[t_dir];
  const Int t_size_local = geo.node_site[t_dir];
  const Int t_shift = t_size_local * geo.geon.coor_node[t_dir];
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec(t_size * multiplicity);
  set_zero(vec);
#pragma omp parallel
  {
    std::vector<M> pvec(t_size_local * multiplicity);
    set_zero(pvec);
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Int tl = xl[t_dir];
      qassert(0 <= tl and tl < t_size_local);
      const Vector<M> fvec = f.get_elems_const(xl);
      for (Int m = 0; m < multiplicity; ++m) {
        pvec[tl * multiplicity + m] += fvec[m];
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        for (Int tl = 0; tl < t_size_local; ++tl) {
          const Int tg = tl + t_shift;
          qassert(0 <= tg and tg < t_size);
          for (Int m = 0; m < multiplicity; ++m) {
            vec[tg * multiplicity + m] += pvec[tl * multiplicity + m];
          }
        }
      }
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_glb_sum(const Field<M>& f)
{
  TIMER("field_glb_sum");
  std::vector<M> vec = field_sum(f);
  glb_sum(vec);
  return vec;
}

template <class M>
std::vector<M> field_glb_sum_tslice(const Field<M>& f, const Int t_dir = 3)
// length = t_size * multiplicity
{
  TIMER("field_glb_sum_tslice");
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum(get_data(vec));
  return vec;
}

template <class M>
std::vector<M> field_max(const Field<M>& f)
// length = multiplicity
{
  TIMER("field_max");
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec(multiplicity);
  assign(vec, f.get_elems_const(0));
#pragma omp parallel
  {
    std::vector<M> pvec(multiplicity);
    assign(pvec, vec);
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<M> fvec = f.get_elems_const(xl);
      for (Int m = 0; m < multiplicity; ++m) {
        pvec[m] = std::max(pvec[m], fvec[m]);
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        for (Int m = 0; m < multiplicity; ++m) {
          vec[m] = std::max(vec[m], pvec[m]);
        }
      }
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_glb_max(const Field<M>& f)
{
  TIMER("field_glb_max");
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_max(f);
  Qassert((Long)vec.size() == (Long)multiplicity);
  glb_max(vec);
  return vec;
}

template <class M>
std::vector<M> field_min(const Field<M>& f)
// length = multiplicity
{
  TIMER("field_min");
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec(multiplicity);
  assign(vec, f.get_elems_const(0));
#pragma omp parallel
  {
    std::vector<M> pvec(multiplicity);
    assign(pvec, vec);
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Vector<M> fvec = f.get_elems_const(xl);
      for (Int m = 0; m < multiplicity; ++m) {
        pvec[m] = std::min(pvec[m], fvec[m]);
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        for (Int m = 0; m < multiplicity; ++m) {
          vec[m] = std::min(vec[m], pvec[m]);
        }
      }
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_glb_min(const Field<M>& f)
{
  TIMER("field_glb_min");
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_min(f);
  Qassert((Long)vec.size() == (Long)multiplicity);
  glb_min(vec);
  return vec;
}

template <class M>
std::vector<M> field_project_mom(const Field<M>& f, const CoordinateD& mom)
// mom is in lattice unit (1/a)
// project to component with momentum 'mom'
// use glb_sum_double_vec to perform glb_sum
{
  TIMER("field_project_mom");
  const Geometry& geo = f.geo();
  std::vector<M> ret(f.multiplicity);
  set_zero(ret);
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    RealD phase = 0;
    for (Int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    const ComplexD factor = qpolar(1.0, -phase);
    const Vector<M> v = f.get_elems_const(xl);
    for (Int m = 0; m < f.multiplicity; ++m) {
      M x = v[m];
      x *= factor;
      ret[m] += x;
    }
  }
  glb_sum(ret);
  return ret;
}

template <class M>
std::vector<M> field_get_elems(const Field<M>& f, const Coordinate& xg)
// xg is same on all the nodes
{
  const Geometry& geo = f.geo();
  const Coordinate xg_r = mod(xg, geo.total_site());
  const Coordinate xl = geo.coordinate_l_from_g(xg_r);
  std::vector<M> ret(f.multiplicity);
  if (geo.is_local(xl)) {
    assign(ret, f.get_elems_const(xl));
  } else {
    set_zero(ret);
  }
  glb_sum(get_data_char(ret));
  return ret;
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg, const Int m)
// xg is same on all the nodes
{
  const Geometry& geo = f.geo();
  const Coordinate xg_r = mod(xg, geo.total_site());
  const Coordinate xl = geo.coordinate_l_from_g(xg_r);
  M ret;
  if (geo.is_local(xl)) {
    ret = f.get_elem(xl, m);
  } else {
    set_zero(ret);
  }
  glb_sum(get_data_char(ret));
  return ret;
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg)
// xg is same on all the nodes
{
  Qassert(f.multiplicity == 1);
  return field_get_elem(f, xg, 0);
}

template <class M>
void field_set_elems(Field<M>& f, const Coordinate& xg, const Vector<M> val)
// xg do not need to be the same on all the nodes
{
  const Geometry& geo = f.geo();
  const Coordinate xg_r = mod(xg, geo.total_site());
  const Coordinate xl = geo.coordinate_l_from_g(xg_r);
  if (geo.is_local(xl)) {
    assign(f.get_elems(xl), val);
  }
}

template <class M>
void field_set_elem(Field<M>& f, const Coordinate& xg, const Int m, const M& val)
// xg do not need to be the same on all the nodes
{
  const Geometry& geo = f.geo();
  const Coordinate xg_r = mod(xg, geo.total_site());
  const Coordinate xl = geo.coordinate_l_from_g(xg_r);
  if (geo.is_local(xl)) {
    f.get_elem(xl, m) = val;
  }
}

template <class M>
void field_set_elem(Field<M>& f, const Coordinate& xg, const M& val)
// xg do not need to be the same on all the nodes
{
  const Geometry& geo = f.geo();
  const Coordinate xg_r = mod(xg, geo.total_site());
  const Coordinate xl = geo.coordinate_l_from_g(xg_r);
  if (geo.is_local(xl)) {
    f.get_elem(xl) = val;
  }
}

template <class M>
void split_fields(std::vector<Handle<Field<M> > >& vec, const Field<M>& f)
// fields in vector will be reinitialized to have the same geo and multiplicity
{
  TIMER("split_fields");
  Qassert(vec.size() >= 1);
  Qassert(is_initialized(f));
  const Long nf = vec.size();
  const Geometry& geo = f.geo();
  Qassert(geo.is_only_local);
  const Int multiplicity = f.multiplicity;
  const Int multiplicity_v = multiplicity / nf;
  Qassert(multiplicity_v * nf == multiplicity);
  for (Long i = 0; i < nf; ++i) {
    Field<M>& f1 = vec[i]();
    f1.init();
    f1.init(geo, multiplicity_v);
    const Int m_offset = i * multiplicity_v;
    qacc_for(index, geo.local_volume(), {
      const Vector<M> fv = f.get_elems_const(index);
      Vector<M> f1v = f1.get_elems(index);
      for (Int m = 0; m < multiplicity_v; ++m) {
        f1v[m] = fv[m + m_offset];
      }
    });
  }
}

template <class M>
void merge_fields(Field<M>& f, const std::vector<ConstHandle<Field<M> > >& vec)
// fields in vector should have the same geo and multiplicity
{
  TIMER("merge_fields");
  Qassert(vec.size() >= 1);
  Qassert(not vec[0].null());
  Qassert(is_initialized(vec[0]()));
  const Long nf = vec.size();
  const Geometry& geo = vec[0]().geo();
  Qassert(geo.is_only_local);
  const Int multiplicity_v = vec[0]().multiplicity;
  const Int multiplicity = nf * multiplicity_v;
  for (Long i = 1; i < nf; ++i) {
    Qassert(geo == vec[i]().geo());
    Qassert(multiplicity_v == vec[i]().multiplicity);
  }
  f.init(geo, multiplicity);
  for (Long i = 0; i < nf; ++i) {
    const Field<M>& f1 = vec[i]();
    const Int m_offset = i * multiplicity_v;
    qacc_for(index, geo.local_volume(), {
      Vector<M> fv = f.get_elems(index);
      const Vector<M> f1v = f1.get_elems_const(index);
      for (Int m = 0; m < multiplicity_v; ++m) {
        fv[m + m_offset] = f1v[m];
      }
    });
  }
}

template <class M>
void merge_fields_ms(Field<M>& f, const std::vector<ConstHandle<Field<M> > >& vec, const std::vector<Int>& m_vec)
// f.get_elem(x, m) = vec[m].get_elem(x, m_vec[m])
{
  TIMER("merge_fields_ms");
  Qassert(vec.size() >= 1);
  Qassert(not vec[0].null());
  Qassert(is_initialized(vec[0]()));
  const Long multiplicity = vec.size();
  Qassert(multiplicity == (Long)m_vec.size());
  const Geometry geo = vec[0]().geo();
  f.init(geo, multiplicity);
  for (Long m = 0; m < multiplicity; ++m) {
    const Field<M>& f1 = vec[m]();
    const Geometry& geo_v = f1.geo();
    Qassert(geo_v.is_only_local);
    check_matching_geo(geo_v, geo);
  }
  qthread_for(index, geo.local_volume(), {
    Vector<M> fv = f.get_elems(index);
    for (Int m = 0; m < multiplicity; ++m) {
      const Field<M>& f1 = vec[m]();
      const Vector<M> f1v = f1.get_elems_const(index);
      fv[m] = f1v[m_vec[m]];
    }
  });
}

template <class M>
void field_shift_dir(Field<M>& f, const Field<M>& f1, const Int dir,
                     const Int shift)
// shift f1 in 'dir' direction for 'shift' steps
// roughly f[xg + shift * nvec] = f1[xg]
// where
// Coordinate nvec;
// nvec[dir] = 1;
{
  TIMER("field_shift_dir");
  Qassert(0 <= dir and dir < 4);
  const Geometry geo = geo_resize(f1.geo());
  const Int multiplicity = f1.multiplicity;
  const Coordinate total_site = geo.total_site();
  f.init(geo, multiplicity);
  Qassert(is_matching_geo(f.geo(), f1.geo()));
  Coordinate nvec;
  nvec[dir] = 1;
  Field<M> tmp, tmp1;
  tmp.init(geo, multiplicity);
  tmp1.init(geo, multiplicity);
  tmp1 = f1;
  for (Int i = 0; i < geo.geon.size_node[dir]; ++i) {
#pragma omp parallel for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Coordinate xg1 =
          mod(xg - (shift + i * geo.node_site[dir]) * nvec, total_site);
      const Coordinate xl1 = geo.coordinate_l_from_g(xg1);
      if (geo.is_local(xl1)) {
        assign(f.get_elems(xl), tmp1.get_elems_const(xl1));
      }
    }
    if (i < geo.geon.size_node[dir] - 1) {
      get_data_plus_mu(get_data(tmp), get_data(tmp1), dir);
      qswap(tmp1, tmp);
    }
  }
}

template <class M>
void field_shift_steps(Field<M>& f, const Field<M>& f1, const Coordinate& shift)
// shift f1 with 'shift'
{
  TIMER("field_shift_steps");
  Field<M> tmp, tmp1;
  field_shift_dir(tmp, f1, 0, shift[0]);
  field_shift_dir(tmp1, tmp, 1, shift[1]);
  field_shift_dir(tmp, tmp1, 2, shift[2]);
  field_shift_dir(f, tmp, 3, shift[3]);
}

template <typename M1, typename M2>
void vector_field_cast(std::vector<M1 >& res, const std::vector<M2 >& src)
{
  res.resize(0);
  const Long num_field = src.size();
  if(num_field == 0){return ;}
  res.resize(num_field);
  qfor(id_field, num_field, {
    res[id_field].set_view_cast(src[id_field]);
  });
}

template <class M>
void vector_to_acc(vector<M >& res, const std::vector<M >& src)
{
  res.resize(0);
  if(src.size() == 0){
    return ;
  }
  const Long num_field = src.size();
  res.resize_zero(num_field, MemType::Cpu);
  qfor(id_field, num_field, {
    res[id_field].set_view(src[id_field]);
  });
  res.set_mem_type(src[0].get_mem_type());
}

template <class M>
void field_shift_local(Field<M>& fr, const Field<M>& fs, const Coordinate& shift)
// fr and fs must be different
// node process local shift
// roughly fr[(xl + shift) % node_site] == fs[xl]
{
  TIMER("field_shift_no_comm");
  const Geometry& geo = fs.geo();
  if (shift == Coordinate()) {
    fr = fs;
    return;
  }
  Field<M> f0;f0.set_view(fs);
  const MemType mem = fs.field.mem_type;
  qmem_for(index, geo.local_volume(), mem, {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xl_s = mod(xl + shift, geo.node_site);
    Vector<M> fv = fr.get_elems(xl_s);
    const Vector<M> f0v = f0.get_elems_const(xl);
    assign(fv, f0v);
  })
}

template <class M>
void field_shift_directT(std::vector<Field<M> >& fr, const std::vector<Field<M> >& fs,
                        const Coordinate& shift, std::vector<vector<M> >& to_bufL)
// fr and fs can be the same
// shift fs with 'shift'
// roughly f[(xg + shift) % total_site] == fs[xg]
// use the fact that the ordering does not change
// UNLESS in some direction there is only one node,
// THEN periodic boundary condition shall mess up the order
// JUST do NOT shift in such direction
// shift it afterwards (in the final step of this function)
// to_bufL : buffers which should have the mem_type for communication
{
  TIMER("field_shift_direct");
  if (fs.size() == 0) {
    fr.resize(0);
    return;
  }
  const unsigned int Nvec = fs.size();
  //
  const Geometry& geo = fs[0].geo();
  Qassert(geo.is_only_local);
  const MemType mem_type = fs[0].get_mem_type();
  const MemType mem_comm = get_comm_mem_type(mem_type);  // for buffer mem_type
  /*
    mem_type transfer for to_bufL currently not supported
    May try to have another buffer to hold the memeory for communications
  */
  Qassert(is_same_comm_mem_type(mem_type, mem_comm));
  //
  Qassert(fr.size() == Nvec);
  const Long MULTI = fs[0].multiplicity;
  for (unsigned int iv = 0; iv < Nvec; iv++) {
    Qassert(fs[iv].geo() == fs[0].geo() and
            fs[iv].field.mem_type == fs[0].field.mem_type);
    Qassert(fs[iv].multiplicity == fs[0].multiplicity);
    Qassert(fr[iv].initialized);  // in case a bad reference ?
    Qassert(fr[iv].geo() == fs[0].geo() and
            fr[iv].multiplicity == fs[iv].multiplicity);
  }
  //
  const Int num_node = geo.geon.num_node;
  const Coordinate& node_site = geo.node_site;
  const Coordinate& size_node = geo.geon.size_node;
  const Coordinate total_site = geo.total_site();
  Coordinate shift_corrected = shift;
  Coordinate shift_remain;
  for (Int mu = 0; mu < 4; ++mu) {
    Qassert(size_node[mu] >= 1);
    shift_corrected[mu] = shift[mu] % total_site[mu];
    if (size_node[mu] == 1) {
      shift_remain[mu] = shift_corrected[mu];
      shift_corrected[mu] = 0;
    }
  }
  if (shift_corrected == Coordinate()) {
    for (unsigned int iv = 0; iv < Nvec; iv++) {
      if (get_data(fr[iv]).data() != get_data(fs[iv]).data()) {
        field_shift_local(fr[iv], fs[iv], shift);
      } else {
        const Long Nd = MULTI * geo.local_volume();
        if (to_bufL.size() == 0) {
          to_bufL.resize(1);
          to_bufL[0].resize(Nd, mem_comm);
        }
        Qassert(to_bufL[0].size() >= Nd);
        Qassert(is_same_comm_mem_type(to_bufL[0].mem_type, mem_comm));
        Vector<M> v(&to_bufL[0][0], Nd);
        //
        Field<M> buf;
        buf.initialized = true;
        buf.geo.set_view(geo);
        buf.multiplicity = MULTI;
        buf.mem_order = fs[iv].mem_order;
        buf.field.set_view(v);
        buf = fs[iv];
        field_shift_local(fr[iv], buf, shift);
      }
    }
    return;
  }
  std::vector<Long> to_send_size(num_node, 0);
  std::vector<Long> to_recv_size(num_node, 0);
  //  Uvm
  FieldM<Long, 2> f_send_idx, f_recv_idx;  // id_node, idx_for_that_node
  {
    TIMER("field_shift_direct setup");
    f_send_idx.init(geo);
    f_recv_idx.init(geo);
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Coordinate xg_send = mod(xg + shift_corrected, total_site);
      const Coordinate xg_recv = mod(xg - shift_corrected, total_site);
      const Long id_node_send =
          index_from_coordinate(xg_send / node_site, size_node);
      const Long id_node_recv =
          index_from_coordinate(xg_recv / node_site, size_node);
      Vector<Long> fsv = f_send_idx.get_elems(index);
      Vector<Long> frv = f_recv_idx.get_elems(index);
      fsv[0] = id_node_send;
      frv[0] = id_node_recv;
      fsv[1] = to_send_size[id_node_send];
      frv[1] = to_recv_size[id_node_recv];
      to_send_size[id_node_send] += 1;
      to_recv_size[id_node_recv] += 1;
    }
  }
  Long N_send = 0;
  Long N_recv = 0;
  // Uvm
  vector<Long> to_send_offset(num_node);
  vector<Long> to_recv_offset(num_node);
  {
    TIMER("field_shift_direct setup");
    for (Int i = 0; i < num_node; ++i) {
      to_send_offset[i] = N_send;
      to_recv_offset[i] = N_recv;
      const Long size_s = to_send_size[i];
      const Long size_r = to_recv_size[i];
      N_send += size_s;
      N_recv += size_r;
    }
  }
  //
  if (to_bufL.size() != 2 * Nvec) {
    to_bufL.resize(2 * Nvec);
  }
  for (unsigned int iv = 0; iv < to_bufL.size(); iv++) {
    const Long Nd = iv < Nvec ? N_send : N_recv;
    if (to_bufL[iv].mem_type != mem_comm or to_bufL[iv].size() < Nd * MULTI) {
      to_bufL[iv].resize_zero(Nd * MULTI,
                              mem_comm);  // may detach from the given buffer
    }
    Qassert(is_same_comm_mem_type(to_bufL[0].mem_type, mem_comm));
  }
  vector<vector<M>> bufh;
  vector_to_acc(bufh, to_bufL);
  // copy to send
  {
    TIMER("field_shift_direct copy");
    vector<Field<M>> fh1;
    vector_to_acc(fh1, fs);
    qmem_for(index, geo.local_volume(), mem_comm, {
      const Vector<Long> fsv = f_send_idx.get_elems_const(index);
      const Int id_node = fsv[0];
      const Long offset = fsv[1] * MULTI;
      const Long off_send = to_send_offset[id_node] * MULTI;
      for (unsigned int iv = 0; iv < Nvec; iv++) {
        const Vector<M> fv = fh1[iv].get_elems_const(index);
        Vector<M> to_send(&bufh[iv][off_send + offset], MULTI);
        assign(to_send, fv);
      }
    });
  }
  // comm send to recv
  {
    TIMER_FLOPS("field_shift_direct-commT");
    std::vector<MPI_Request> reqs_send;
    std::vector<MPI_Request> reqs_recv;
    const Int rank = qlat::get_id_node();
    for (unsigned int iv = 0; iv < Nvec; iv++) {
      timer.flops += geo.local_volume() * (Long)MULTI * (Long)sizeof(M);
      const Int mpi_tag = 11 + iv;
      for (Int id_node = 0; id_node < num_node; id_node++) {
        const Long off_send = to_send_offset[id_node] * MULTI;
        const Long off_recv = to_recv_offset[id_node] * MULTI;
        if (to_recv_size[id_node] > 0 and id_node != rank) {
          mpi_irecv(&to_bufL[Nvec + iv][off_recv],
                    to_recv_size[id_node] * MULTI * sizeof(M), MPI_BYTE,
                    id_node, mpi_tag, get_comm(), reqs_recv);
        }
        //
        if (to_send_size[id_node] > 0 and id_node != rank) {
          mpi_isend(&to_bufL[iv][off_send],
                    to_send_size[id_node] * MULTI * sizeof(M), MPI_BYTE,
                    id_node, mpi_tag, get_comm(), reqs_send);
        }
      }
    }
    // within node copy
    {
      const Int id_node = rank;
      const size_t Nd = to_recv_size[id_node] * MULTI * sizeof(M);
      if (Nd > 0) {
        for (unsigned int iv = 0; iv < Nvec; iv++) {
          TIMER("field_shift_direct node copy");
          const Long off_send = to_send_offset[id_node] * MULTI;
          const Long off_recv = to_recv_offset[id_node] * MULTI;
          Qassert(Nd == to_send_size[id_node] * MULTI * sizeof(M));
          copy_mem(&to_bufL[Nvec + iv][off_recv], mem_comm,
                   &to_bufL[iv][off_send], mem_comm, Nd);
        }
      }
    }
    //
    mpi_waitall(reqs_recv);  ////receive done and write
    mpi_waitall(reqs_send);
  }
  // copy to f
  {
    TIMER("field_shift_direct copy");
    vector<Field<M>> fh0;
    vector_to_acc(fh0, fr);
    qmem_for(index, geo.local_volume(), mem_comm, {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xl_s = mod(xl + shift_remain, total_site);
      const Long index_s = geo.index_from_coordinate(xl_s);
      const Vector<Long> frv = f_recv_idx.get_elems_const(index);
      const Int id_node = frv[0];
      const Long offset = frv[1] * MULTI;
      const Long off_recv = to_recv_offset[id_node] * MULTI;
      for (unsigned int iv = 0; iv < Nvec; iv++) {
        Vector<M> fv = fh0[iv].get_elems(index_s);
        const Vector<M> to_recv(&bufh[Nvec + iv][off_recv + offset], MULTI);
        assign(fv, to_recv);
      }
    });
  }
}

template <class M>
void field_shift_direct(std::vector<Field<M> >& fr, const std::vector<Field<M> >& fs,
                        const Coordinate& shift, std::vector<vector<M> >& to_bufL)
{
  if (sizeof(M) % sizeof(Long) == 0) {
    std::vector<Field<Long>> f0;
    std::vector<Field<Long>> f1;
    std::vector<vector<Long>> tb;
    vector_field_cast(f0, fr);
    vector_field_cast(f1, fs);
    vector_field_cast(tb, to_bufL);
    field_shift_directT(f0, f1, shift, tb);
  } else {
    std::vector<Field<Char>> f0;
    std::vector<Field<Char>> f1;
    std::vector<vector<Char>> tb;
    vector_field_cast(f0, fr);
    vector_field_cast(f1, fs);
    vector_field_cast(tb, to_bufL);
    field_shift_directT(f0, f1, shift, tb);
  }
}

template <class M>
void field_shift_direct(Field<M>& fr, const Field<M>& fs,
                        const Coordinate& shift)
{
  if(!fr.initialized){fr.init(fs.geo(), fs.multiplicity);}
  std::vector<Field<M> > frL;
  std::vector<Field<M> > fsL;
  frL.resize(1);
  fsL.resize(1);
  frL[0].set_view(fr);
  fsL[0].set_view(fs);
  std::vector<vector<M> > to_bufL;
  field_shift_direct(frL, fsL, shift, to_bufL);
}

template <class M>
void field_shift(Field<M>& f, const Field<M>& f1, const Coordinate& shift)
{
  // field_shift_steps(f, f1, shift);
  // field_shift_shuffle(f, f1, shift);
  field_shift_direct(f, f1, shift);
}

// --------------------

void set_xg_field(Field<Int>& f, const Geometry& geo_);

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                        \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator+=                      \
      <TYPENAME>(Field<TYPENAME> & f, const Field<TYPENAME>& f1);             \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator-=                      \
      <TYPENAME>(Field<TYPENAME> & f, const Field<TYPENAME>& f1);             \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                     \
      Field<TYPENAME>& f, const Field<RealD>& f_factor);                     \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                     \
      Field<TYPENAME>& f, const Field<ComplexD>& f_factor);                   \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                      \
      <TYPENAME>(Field<TYPENAME> & f, const double factor);                   \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                      \
      <TYPENAME>(Field<TYPENAME> & f, const ComplexD& factor);                \
                                                                              \
  QLAT_EXTERN template RealD qnorm<TYPENAME>(const Field<TYPENAME>& f);       \
                                                                              \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(Field<RealD> & f,           \
                                                  const Field<TYPENAME>& f1); \
                                                                              \
  QLAT_EXTERN template RealD qnorm_double<TYPENAME>(                          \
      const Field<TYPENAME>& f1, const Field<TYPENAME>& f2);                  \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum<TYPENAME>(             \
      const Field<TYPENAME>& f);                                              \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum_tslice<TYPENAME>(      \
      const Field<TYPENAME>& f, const Int t_dir);                             \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum<TYPENAME>(         \
      const Field<TYPENAME>& f);                                              \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum_tslice<TYPENAME>(  \
      const Field<TYPENAME>& f, const Int t_dir);                             \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_project_mom<TYPENAME>(     \
      const Field<TYPENAME>& f, const CoordinateD& mom);                      \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_get_elems<TYPENAME>(       \
      const Field<TYPENAME>& f, const Coordinate& xg);                        \
                                                                              \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                     \
      const Field<TYPENAME>& f, const Coordinate& xg, const Int m);           \
                                                                              \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                     \
      const Field<TYPENAME>& f, const Coordinate& xg);                        \
                                                                              \
  QLAT_EXTERN template void field_set_elems<TYPENAME>(                        \
      Field<TYPENAME> & f, const Coordinate& xg, const Vector<TYPENAME> val); \
                                                                              \
  QLAT_EXTERN template void field_set_elem<TYPENAME>(                         \
      Field<TYPENAME> & f, const Coordinate& xg, const Int m,                 \
      const TYPENAME& val);                                                   \
                                                                              \
  QLAT_EXTERN template void field_set_elem<TYPENAME>(                         \
      Field<TYPENAME> & f, const Coordinate& xg, const TYPENAME& val);        \
                                                                              \
  QLAT_EXTERN template void split_fields<TYPENAME>(                           \
      std::vector<Handle<Field<TYPENAME>>> & vec, const Field<TYPENAME>& f);  \
                                                                              \
  QLAT_EXTERN template void merge_fields<TYPENAME>(                           \
      Field<TYPENAME> & f,                                                    \
      const std::vector<ConstHandle<Field<TYPENAME>>>& vec);                  \
                                                                              \
  QLAT_EXTERN template void merge_fields_ms<TYPENAME>(                        \
      Field<TYPENAME> & f,                                                    \
      const std::vector<ConstHandle<Field<TYPENAME>>>& vec,                   \
      const std::vector<Int>& m_vec);                                         \
                                                                              \
  QLAT_EXTERN template void field_shift_dir<TYPENAME>(                        \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1, const Int dir,          \
      const Int shift);                                                       \
                                                                              \
  QLAT_EXTERN template void field_shift_steps<TYPENAME>(                      \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1,                         \
      const Coordinate& shift);                                               \
                                                                              \
  QLAT_EXTERN template void field_shift_direct<TYPENAME>(                     \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1,                         \
      const Coordinate& shift);                                               \
                                                                              \
  QLAT_EXTERN template void field_shift<TYPENAME>(                            \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1, const Coordinate& shift)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#define QLAT_EXTERN_TEMPLATE_2(TYPENAME1, TYPENAME2)    \
  QLAT_EXTERN template void assign(Field<TYPENAME1>& f, \
                                   const Field<TYPENAME2>& f1)

QLAT_CALL_WITH_TYPES_2(QLAT_EXTERN_TEMPLATE_2);
#undef QLAT_EXTERN_TEMPLATE_2

QLAT_EXTERN template std::vector<RealD> field_max(const Field<RealD>& f);
QLAT_EXTERN template std::vector<RealF> field_max(const Field<RealF>& f);
QLAT_EXTERN template std::vector<Long> field_max(const Field<Long>& f);
QLAT_EXTERN template std::vector<Int> field_max(const Field<Int>& f);

QLAT_EXTERN template std::vector<RealD> field_min(const Field<RealD>& f);
QLAT_EXTERN template std::vector<RealF> field_min(const Field<RealF>& f);
QLAT_EXTERN template std::vector<Long> field_min(const Field<Long>& f);
QLAT_EXTERN template std::vector<Int> field_min(const Field<Int>& f);

QLAT_EXTERN template std::vector<RealD> field_glb_max(const Field<RealD>& f);
QLAT_EXTERN template std::vector<RealF> field_glb_max(const Field<RealF>& f);
QLAT_EXTERN template std::vector<Long> field_glb_max(const Field<Long>& f);
QLAT_EXTERN template std::vector<Int> field_glb_max(const Field<Int>& f);

QLAT_EXTERN template std::vector<RealD> field_glb_min(const Field<RealD>& f);
QLAT_EXTERN template std::vector<RealF> field_glb_min(const Field<RealF>& f);
QLAT_EXTERN template std::vector<Long> field_glb_min(const Field<Long>& f);
QLAT_EXTERN template std::vector<Int> field_glb_min(const Field<Int>& f);

#undef QLAT_EXTERN

}  // namespace qlat
