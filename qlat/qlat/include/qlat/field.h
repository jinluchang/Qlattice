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
  qassert(is_matching_geo(f.geo(), f1.geo()));
  qassert(f.multiplicity == f1.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < f.multiplicity; ++m) {
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
  qassert(is_matching_geo(f.geo(), f1.geo()));
  if (f.multiplicity != f1.multiplicity) {
    qerr(fname +
         ssprintf(": f.mult=%d ; f1.mult=%d", f.multiplicity, f1.multiplicity));
  }
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < f.multiplicity; ++m) {
      v[m] -= v1[m];
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Field<double>& f_factor)
{
  TIMER("field_operator*=(F,FD)");
  qassert(is_matching_geo(f.geo(), f_factor.geo()));
  const int multiplicity_f = f_factor.multiplicity;
  qassert(multiplicity_f == 1 or multiplicity_f == f.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const double fac = f_factor.get_elem(xl);
      for (int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == f.multiplicity);
      Vector<double> fac = f_factor.get_elems_const(xl);
      for (int m = 0; m < f.multiplicity; ++m) {
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
  qassert(is_matching_geo(f.geo(), f_factor.geo()));
  const int multiplicity_f = f_factor.multiplicity;
  qassert(multiplicity_f == 1 or multiplicity_f == f.multiplicity);
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const ComplexD fac = f_factor.get_elem(xl);
      for (int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == f.multiplicity);
      const Vector<ComplexD> fac = f_factor.get_elems_const(xl);
      for (int m = 0; m < f.multiplicity; ++m) {
        v[m] *= fac[m];
      }
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const double factor)
{
  TIMER("field_operator*=(F,D)");
  qacc_for(index, f.geo().local_volume(), {
    const Geometry& geo = f.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < f.multiplicity; ++m) {
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
    for (int m = 0; m < f.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
double qnorm(const Field<M>& f)
{
  TIMER("qnorm(f)");
  const Geometry& geo = f.geo();
  double sum = 0.0;
#pragma omp parallel
  {
    double psum = 0.0;
#pragma omp for nowait
    for (Long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate x = geo.coordinate_from_index(index);
      const Vector<M> fx = f.get_elems_const(x);
      for (int m = 0; m < f.multiplicity; ++m) {
        psum += qnorm(fx[m]);
      }
    }
    for (int i = 0; i < omp_get_num_threads(); ++i) {
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
double qnorm_double(const Field<M>& f1, const Field<M>& f2)
{
  const Geometry& geo = f1.geo();
  qassert(geo.is_only_local);
  qassert(geo == f2.geo());
  double sum = qnorm_double(get_data(f1), get_data(f2));
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
  qassert(geo.is_only_local);
  f.init();
  f.init(geo, f1.multiplicity * sizeof(M) / sizeof(N));
  qassert(f.multiplicity * sizeof(N) == f1.multiplicity * sizeof(M));
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
  const int multiplicity = f.multiplicity;
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
      for (int m = 0; m < multiplicity; ++m) {
        pvec[m] += fvec[m];
      }
    }
    for (Int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
        for (int m = 0; m < multiplicity; ++m) {
          vec[m] += pvec[m];
        }
      }
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_sum_tslice(const Field<M>& f, const int t_dir = 3)
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
  glb_sum_vec(get_data(vec));
  return vec;
}

template <class M>
std::vector<M> field_glb_sum_tslice(const Field<M>& f, const int t_dir = 3)
// length = t_size * multiplicity
{
  TIMER("field_glb_sum_tslice");
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum_vec(get_data(vec));
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
    double phase = 0;
    for (int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    const ComplexD factor = qpolar(1.0, -phase);
    const Vector<M> v = f.get_elems_const(xl);
    for (int m = 0; m < f.multiplicity; ++m) {
      M x = v[m];
      x *= factor;
      ret[m] += x;
    }
  }
  glb_sum_double_vec(get_data(ret));
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
  glb_sum_byte_vec(get_data(ret));
  return ret;
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg, const int m)
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
  glb_sum_byte_vec(get_data_one_elem(ret));
  return ret;
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg)
// xg is same on all the nodes
{
  qassert(f.multiplicity == 1);
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
void field_set_elem(Field<M>& f, const Coordinate& xg, const int m, const M& val)
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
  qassert(vec.size() >= 1);
  qassert(is_initialized(f));
  const Long nf = vec.size();
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  const int multiplicity = f.multiplicity;
  const int multiplicity_v = multiplicity / nf;
  qassert(multiplicity_v * nf == multiplicity);
  for (Long i = 0; i < nf; ++i) {
    Field<M>& f1 = vec[i]();
    f1.init();
    f1.init(geo, multiplicity_v);
    const int m_offset = i * multiplicity_v;
    qacc_for(index, geo.local_volume(), {
      const Vector<M> fv = f.get_elems_const(index);
      Vector<M> f1v = f1.get_elems(index);
      for (int m = 0; m < multiplicity_v; ++m) {
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
  qassert(vec.size() >= 1);
  qassert(not vec[0].null());
  qassert(is_initialized(vec[0]()));
  const Long nf = vec.size();
  const Geometry& geo = vec[0]().geo();
  qassert(geo.is_only_local);
  const int multiplicity_v = vec[0]().multiplicity;
  const int multiplicity = nf * multiplicity_v;
  for (Long i = 1; i < nf; ++i) {
    qassert(geo == vec[i]().geo());
    qassert(multiplicity_v == vec[i]().multiplicity);
  }
  f.init(geo, multiplicity);
  for (Long i = 0; i < nf; ++i) {
    const Field<M>& f1 = vec[i]();
    const int m_offset = i * multiplicity_v;
    qacc_for(index, geo.local_volume(), {
      Vector<M> fv = f.get_elems(index);
      const Vector<M> f1v = f1.get_elems_const(index);
      for (int m = 0; m < multiplicity_v; ++m) {
        fv[m + m_offset] = f1v[m];
      }
    });
  }
}

template <class M>
void merge_fields_ms(Field<M>& f, const std::vector<ConstHandle<Field<M> > >& vec, const std::vector<int> m_vec)
// f.get_elem(x, m) = vec[m].get_elem(x, m_vec[m])
{
  TIMER("merge_fields_ms");
  qassert(vec.size() >= 1);
  qassert(not vec[0].null());
  qassert(is_initialized(vec[0]()));
  const Long multiplicity = vec.size();
  qassert(multiplicity == (Long)m_vec.size());
  const Geometry geo = vec[0]().geo();
  f.init(geo, multiplicity);
  for (Long m = 0; m < multiplicity; ++m) {
    const Field<M>& f1 = vec[m]();
    const Geometry& geo_v = f1.geo();
    qassert(geo_v.is_only_local);
    check_matching_geo(geo_v, geo);
  }
  qthread_for(index, geo.local_volume(), {
    Vector<M> fv = f.get_elems(index);
    for (int m = 0; m < multiplicity; ++m) {
      const Field<M>& f1 = vec[m]();
      const Vector<M> f1v = f1.get_elems_const(index);
      fv[m] = f1v[m_vec[m]];
    }
  });
}

template <class M>
void field_shift_dir(Field<M>& f, const Field<M>& f1, const int dir,
                     const int shift)
// shift f1 in 'dir' direction for 'shift' steps
{
  TIMER("field_shift_dir");
  qassert(0 <= dir and dir < 4);
  const Geometry geo = geo_resize(f1.geo());
  const Int multiplicity = f1.multiplicity;
  const Coordinate total_site = geo.total_site();
  f.init(geo, multiplicity);
  qassert(is_matching_geo(f.geo(), f1.geo()));
  Coordinate nvec;
  nvec[dir] = 1;
  Field<M> tmp, tmp1;
  tmp.init(geo, multiplicity);
  tmp1.init(geo, multiplicity);
  tmp1 = f1;
  for (int i = 0; i < geo.geon.size_node[dir]; ++i) {
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

template <class M>
void field_shift_local(Field<M>& f, const Field<M>& f1, const Coordinate& shift)
// node process local shift
// roughly f[xl + shift % node_site] == f1[xl]
{
  TIMER("field_shift_no_comm");
  if (shift == Coordinate()) {
    f = f1;
    return;
  }
  const Geometry geo = geo_resize(f1.geo());
  Field<M> f0;
  f0.init(geo, f1.multiplicity);
  f0 = f1;
  f.init(geo, f1.multiplicity);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xl_s = mod(xl + shift, geo.node_site);
    Vector<M> fv = f.get_elems(xl_s);
    const Vector<M> f0v = f0.get_elems_const(xl);
    for (int m = 0; m < f.multiplicity; ++m) {
      fv[m] = f0v[m];
    }
  });
}

template <class M>
void field_shift_direct(Field<M>& f, const Field<M>& f1,
                        const Coordinate& shift)
// shift f1 with 'shift'
// use the fact that the ordering does not change
// UNLESS in some direction there is only one node,
// THEN periodic boundary condition shall mess up the order
// JUST do NOT shift in such direction
// shift it afterwards (in the final step of this function)
{
  TIMER("field_shift_direct");
  const Geometry& geo = f1.geo();
  qassert(geo.is_only_local);
  const int num_node = geo.geon.num_node;
  const Coordinate& node_site = geo.node_site;
  const Coordinate& size_node = geo.geon.size_node;
  const Coordinate total_site = geo.total_site();
  Coordinate shift_corrected = shift;
  Coordinate shift_remain;
  for (int mu = 0; mu < 4; ++mu) {
    qassert(size_node[mu] >= 1);
    shift_corrected[mu] = shift[mu] % total_site[mu];
    if (size_node[mu] == 1) {
      shift_remain[mu] = shift_corrected[mu];
      shift_corrected[mu] = 0;
    }
  }
  if (shift_corrected == Coordinate()) {
    field_shift_local(f, f1, shift);
    return;
  }
  std::vector<Long> to_send_size(num_node, 0);
  std::vector<Long> to_recv_size(num_node, 0);
  FieldM<Long, 2> f_send_idx, f_recv_idx;  // id_node, idx_for_that_node
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
  int n_send = 0;
  int n_recv = 0;
  std::vector<vector<M>> to_send(num_node);
  std::vector<vector<M>> to_recv(num_node);
  for (int i = 0; i < num_node; ++i) {
    const Long size_s = to_send_size[i];
    if (size_s > 0) {
      to_send[i].resize(size_s * f1.multiplicity);
      n_send += 1;
    }
    const Long size_r = to_recv_size[i];
    if (size_r > 0) {
      to_recv[i].resize(size_r * f1.multiplicity);
      n_recv += 1;
    }
  }
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Vector<M> fv = f1.get_elems_const(index);
    const Vector<Long> fsv = f_send_idx.get_elems_const(index);
    const int id_node = fsv[0];
    const Long offset = fsv[1] * f1.multiplicity;
    vector<M>& to_send_v = to_send[id_node];
    for (int m = 0; m < f1.multiplicity; ++m) {
      to_send_v[offset + m] = fv[m];
    }
  }
  {
    TIMER_FLOPS("field_shift_direct-comm");
    timer.flops +=
        geo.local_volume() * (Long)f1.multiplicity * (Long)sizeof(M);
    const Long max_elem = 1 + get_max_field_shift_direct_msg_size() / sizeof(M);
    std::vector<MPI_Request> reqs(n_recv + n_send);
    Vector<MPI_Request> recv_reqs(reqs.data(), n_recv);
    Vector<MPI_Request> send_reqs(reqs.data() + n_recv, n_send);
    const int mpi_tag = 11;
    int i_send = 0;
    int i_recv = 0;
    for (int i = 0; i < num_node; ++i) {
      Long size_r = to_recv[i].size();
      if (size_r > 0) {
        const int id_node = i;
        qassert(i_recv < n_recv);
        vector<M>& to_recv_v = to_recv[id_node];
        Long offset = 0;
        while (size_r > max_elem) {
          MPI_Request req;
          MPI_Irecv(&to_recv_v[offset], max_elem * sizeof(M), MPI_BYTE, id_node,
                    mpi_tag, get_comm(), &req);
          reqs.push_back(req);
          offset += max_elem;
          size_r -= max_elem;
        }
        MPI_Irecv(&to_recv_v[offset], size_r * sizeof(M), MPI_BYTE, id_node,
                  mpi_tag, get_comm(), &recv_reqs[i_recv]);
        i_recv += 1;
      }
      Long size_s = to_send[i].size();
      if (size_s > 0) {
        const int id_node = i;
        qassert(i_send < n_send);
        vector<M>& to_send_v = to_send[id_node];
        Long offset = 0;
        while (size_s > max_elem) {
          MPI_Request req;
          MPI_Isend(&to_send_v[offset], max_elem * sizeof(M), MPI_BYTE, id_node,
                    mpi_tag, get_comm(), &req);
          reqs.push_back(req);
          offset += max_elem;
          size_s -= max_elem;
        }
        MPI_Isend(&to_send_v[offset], size_s * sizeof(M), MPI_BYTE, id_node,
                  mpi_tag, get_comm(), &send_reqs[i_send]);
        i_send += 1;
      }
    }
    if (reqs.size() > 0) {
      MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUS_IGNORE);
    }
    sync_node();
  }
  f.init(geo, f1.multiplicity);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xl_s = mod(xl + shift_remain, total_site);
    const Long index_s = geo.index_from_coordinate(xl_s);
    Vector<M> fv = f.get_elems(index_s);
    const Vector<Long> frv = f_recv_idx.get_elems_const(index);
    const int id_node = frv[0];
    const Long offset = frv[1] * f1.multiplicity;
    vector<M>& to_recv_v = to_recv[id_node];
    for (int m = 0; m < f1.multiplicity; ++m) {
      fv[m] = to_recv_v[offset + m];
    }
  }
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
      <TYPENAME>(Field<TYPENAME>& f, const Field<TYPENAME>& f1);              \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator-=                      \
      <TYPENAME>(Field<TYPENAME>& f, const Field<TYPENAME>& f1);              \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                     \
      Field<TYPENAME>& f, const Field<double>& f_factor);                     \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                     \
      Field<TYPENAME>& f, const Field<ComplexD>& f_factor);                   \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                      \
      <TYPENAME>(Field<TYPENAME>& f, const double factor);                    \
                                                                              \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                      \
      <TYPENAME>(Field<TYPENAME>& f, const ComplexD& factor);                 \
                                                                              \
  QLAT_EXTERN template RealD qnorm<TYPENAME>(const Field<TYPENAME>& f);      \
                                                                              \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(Field<RealD> & f,           \
                                                  const Field<TYPENAME>& f1); \
                                                                              \
  QLAT_EXTERN template RealD qnorm_double<TYPENAME>(                         \
      const Field<TYPENAME>& f1, const Field<TYPENAME>& f2);                  \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum<TYPENAME>(             \
      const Field<TYPENAME>& f);                                              \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum_tslice<TYPENAME>(      \
      const Field<TYPENAME>& f, const int t_dir);                             \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum<TYPENAME>(         \
      const Field<TYPENAME>& f);                                              \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum_tslice<TYPENAME>(  \
      const Field<TYPENAME>& f, const int t_dir);                             \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_project_mom<TYPENAME>(     \
      const Field<TYPENAME>& f, const CoordinateD& mom);                      \
                                                                              \
  QLAT_EXTERN template std::vector<TYPENAME> field_get_elems<TYPENAME>(       \
      const Field<TYPENAME>& f, const Coordinate& xg);                        \
                                                                              \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                     \
      const Field<TYPENAME>& f, const Coordinate& xg, const int m);           \
                                                                              \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                     \
      const Field<TYPENAME>& f, const Coordinate& xg);                        \
                                                                              \
  QLAT_EXTERN template void field_set_elems<TYPENAME>(                        \
      Field<TYPENAME> & f, const Coordinate& xg, const Vector<TYPENAME> val); \
                                                                              \
  QLAT_EXTERN template void field_set_elem<TYPENAME>(                         \
      Field<TYPENAME> & f, const Coordinate& xg, const int m,                 \
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
      const std::vector<int> m_vec);                                          \
                                                                              \
  QLAT_EXTERN template void field_shift_dir<TYPENAME>(                        \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1, const int dir,          \
      const int shift);                                                       \
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

#undef QLAT_EXTERN

}  // namespace qlat
