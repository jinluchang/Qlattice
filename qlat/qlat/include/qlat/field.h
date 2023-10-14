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

void set_mom_phase_field(FieldM<Complex, 1>& f, const CoordinateD& mom);

void set_phase_field(FieldM<Complex, 1>& f, const CoordinateD& lmom);

// --------------------

template <class M, class N>
Field<M>& qcast(Field<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  if (x.initialized) {
    const int size = x.geo().multiplicity * sizeof(N);
    x.geo().multiplicity = size / sizeof(M);
    qassert(x.geo().multiplicity * sizeof(M) == size);
  }
  return (Field<M>&)x;
}

template <class M, class N>
const Field<M>& qcast_const(const Field<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  return qcast<M,N>((Field<N>&)x);
}

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
  qassert(is_matching_geo_mult(f.geo(), f1.geo()));
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
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
    f.init(f1.geo());
    set_zero(f);
    f -= f1;
    return f;
  }
  qassert(is_matching_geo_mult(f.geo(), f1.geo()));
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v1 = f1.get_elems_const(xl);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
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
  const Geometry& geo = f.geo();
  const int multiplicity_f = f_factor.geo().multiplicity;
  qassert(multiplicity_f == 1 or multiplicity_f == geo.multiplicity);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const double fac = f_factor.get_elem(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == geo.multiplicity);
      Vector<double> fac = f_factor.get_elems_const(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        v[m] *= fac[m];
      }
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Field<Complex>& f_factor)
{
  TIMER("field_operator*=(F,FC)");
  qassert(is_matching_geo(f.geo(), f_factor.geo()));
  const Geometry& geo = f.geo();
  const int multiplicity_f = f_factor.geo().multiplicity;
  qassert(multiplicity_f == 1 or multiplicity_f == geo.multiplicity);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    if (multiplicity_f == 1) {
      const Complex fac = f_factor.get_elem(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        v[m] *= fac;
      }
    } else {
      qassert(multiplicity_f == geo.multiplicity);
      const Vector<Complex> fac = f_factor.get_elems_const(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
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
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
const Field<M>& operator*=(Field<M>& f, const Complex& factor)
{
  TIMER("field_operator*=(F,C)");
  const Geometry& geo = f.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = f.get_elems(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      v[m] *= factor;
    }
  });
  return f;
}

template <class M>
double qnorm(const Field<M>& f)
{
  const Geometry& geo = f.geo();
  double sum = 0.0;
#pragma omp parallel
  {
    double psum = 0.0;
#pragma omp for nowait
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate x = geo.coordinate_from_index(index);
      const Vector<M> fx = f.get_elems_const(x);
      for (int m = 0; m < geo.multiplicity; ++m) {
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
qacc long get_data_size(const Field<M>& f)
// NOT including the expended parts, only local volume data size
// only size on one node
{
  return f.geo().local_volume() * f.geo().multiplicity * sizeof(M);
}

template <class M>
void qswap(Field<M>& f1, Field<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  qswap(f1.geo, f2.geo);
  qswap(f1.field, f2.field);
}

// --------------------

template <class M, class N>
void assign(Field<N>& f, const Field<M>& f1)
{
  const Geometry& geo1 = f1.geo();
  qassert(geo1.is_only_local);
  const Geometry geo =
      geo_remult(geo1, geo1.multiplicity * sizeof(M) / sizeof(N));
  f.init();
  f.init(geo);
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
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec(multiplicity);
  set_zero(vec);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> fvec = f.get_elems_const(xl);
    for (int m = 0; m < multiplicity; ++m) {
      vec[m] += fvec[m];
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
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec(t_size * multiplicity);
  set_zero(vec);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<M> fvec = f.get_elems_const(xl);
    for (int m = 0; m < multiplicity; ++m) {
      vec[xg[t_dir] * multiplicity + m] += fvec[m];
    }
  }
  return vec;
}

template <class M>
std::vector<M> field_glb_sum_double(const Field<M>& f)
{
  TIMER("field_glb_sum_double");
  std::vector<M> vec = field_sum(f);
  glb_sum_double_vec(get_data(vec));
  return vec;
}

template <class M>
std::vector<M> field_glb_sum_long(const Field<M>& f)
{
  TIMER("field_glb_sum_long");
  std::vector<M> vec = field_sum(f);
  glb_sum_long_vec(get_data(vec));
  return vec;
}

template <class M>
std::vector<std::vector<M> > field_glb_sum_tslice_double(const Field<M>& f,
                                                         const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_double");
  const Geometry& geo = f.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum_double_vec(get_data(vec));
  std::vector<std::vector<M> > vecs(t_size);
  for (int t = 0; t < t_size; ++t) {
    vecs[t].resize(multiplicity);
    for (int m = 0; m < multiplicity; ++m) {
      vecs[t][m] = vec[t * multiplicity + m];
    }
  }
  return vecs;
}

template <class M>
std::vector<std::vector<M> > field_glb_sum_tslice_long(const Field<M>& f,
                                                       const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_long");
  const Geometry& geo = f.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum_long_vec(get_data(vec));
  std::vector<std::vector<M> > vecs(t_size);
  for (int t = 0; t < t_size; ++t) {
    vecs[t].resize(multiplicity);
    for (int m = 0; m < multiplicity; ++m) {
      vecs[t][m] = vec[t * multiplicity + m];
    }
  }
  return vecs;
}

template <class M>
std::vector<M> field_project_mom(const Field<M>& f, const CoordinateD& mom)
// mom is in lattice unit (1/a)
// project to component with momentum 'mom'
// use glb_sum_double_vec to perform glb_sum
{
  TIMER("field_project_mom");
  const Geometry& geo = f.geo();
  std::vector<M> ret(geo.multiplicity);
  set_zero(ret);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0;
    for (int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    const Complex factor = qpolar(1.0, -phase);
    const Vector<M> v = f.get_elems_const(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
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
  std::vector<M> ret(geo.multiplicity);
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
  const Geometry& geo = f.geo();
  qassert(geo.multiplicity == 1);
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
  const long nf = vec.size();
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  const int multiplicity = geo.multiplicity;
  const int multiplicity_v = multiplicity / nf;
  qassert(multiplicity_v * nf == multiplicity);
  const Geometry geo_v = geo_reform(geo, multiplicity_v);
  for (long i = 0; i < nf; ++i) {
    Field<M>& f1 = vec[i]();
    f1.init();
    f1.init(geo_v);
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
  const long nf = vec.size();
  const Geometry& geo_v = vec[0]().geo();
  qassert(geo_v.is_only_local);
  const int multiplicity_v = geo_v.multiplicity;
  const int multiplicity = nf * multiplicity_v;
  const Geometry geo = geo_reform(geo_v, multiplicity);
  for (long i = 1; i < nf; ++i) {
    qassert(geo_v == vec[i]().geo());
  }
  f.init(geo);
  for (long i = 0; i < nf; ++i) {
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
  const long multiplicity = vec.size();
  qassert(multiplicity == (long)m_vec.size());
  const Geometry geo = geo_reform(vec[0]().geo(), multiplicity);
  f.init(geo);
  for (long m = 0; m < multiplicity; ++m) {
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
  const Coordinate total_site = geo.total_site();
  f.init(geo);
  qassert(is_matching_geo_mult(f.geo(), f1.geo()));
  Coordinate nvec;
  nvec[dir] = 1;
  Field<M> tmp, tmp1;
  tmp.init(geo);
  tmp1.init(geo);
  tmp1 = f1;
  for (int i = 0; i < geo.geon.size_node[dir]; ++i) {
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
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
  for (int mu = 0; mu < 4; ++mu) {
    qassert(size_node[mu] >= 1);
    if (size_node[mu] == 1) {
      shift_corrected[mu] = 0;
    }
  }
  const Coordinate shift_remain = shift - shift_corrected;
  std::vector<long> to_send_size(num_node, 0);
  std::vector<long> to_recv_size(num_node, 0);
  FieldM<long, 2> f_send_idx, f_recv_idx;  // id_node, idx_for_that_node
  f_send_idx.init(geo);
  f_recv_idx.init(geo);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xg_send = mod(xg + shift_corrected, total_site);
    const Coordinate xg_recv = mod(xg - shift_corrected, total_site);
    const long id_node_send =
        index_from_coordinate(xg_send / node_site, size_node);
    const long id_node_recv =
        index_from_coordinate(xg_recv / node_site, size_node);
    Vector<long> fsv = f_send_idx.get_elems(index);
    Vector<long> frv = f_recv_idx.get_elems(index);
    fsv[0] = id_node_send;
    frv[0] = id_node_recv;
    fsv[1] = to_send_size[id_node_send];
    frv[1] = to_recv_size[id_node_recv];
    to_send_size[id_node_send] += 1;
    to_recv_size[id_node_recv] += 1;
  }
  int n_send = 0;
  int n_recv = 0;
  std::vector<std::vector<M> > to_send(num_node);
  std::vector<std::vector<M> > to_recv(num_node);
  for (int i = 0; i < num_node; ++i) {
    const long size_s = to_send_size[i];
    if (size_s > 0) {
      to_send[i].resize(size_s * geo.multiplicity);
      n_send += 1;
    }
    const long size_r = to_recv_size[i];
    if (size_r > 0) {
      to_recv[i].resize(size_r * geo.multiplicity);
      n_recv += 1;
    }
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Vector<M> fv = f1.get_elems_const(index);
    const Vector<long> fsv = f_send_idx.get_elems_const(index);
    const int id_node = fsv[0];
    const long offset = fsv[1] * geo.multiplicity;
    std::vector<M>& to_send_v = to_send[id_node];
    for (int m = 0; m < geo.multiplicity; ++m) {
      to_send_v[offset + m] = fv[m];
    }
  }
  {
    TIMER_FLOPS("field_shift_direct-comm");
    timer.flops +=
        geo.local_volume() * (long)geo.multiplicity * (long)sizeof(M);
    const long max_elem = 1 + get_max_field_shift_direct_msg_size() / sizeof(M);
    std::vector<MPI_Request> reqs(n_recv + n_send);
    Vector<MPI_Request> recv_reqs(reqs.data(), n_recv);
    Vector<MPI_Request> send_reqs(reqs.data() + n_recv, n_send);
    const int mpi_tag = 11;
    int i_send = 0;
    int i_recv = 0;
    for (int i = 0; i < num_node; ++i) {
      long size_r = to_recv[i].size();
      if (size_r > 0) {
        const int id_node = i;
        qassert(i_recv < n_recv);
        std::vector<M>& to_recv_v = to_recv[id_node];
        long offset = 0;
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
      long size_s = to_send[i].size();
      if (size_s > 0) {
        const int id_node = i;
        qassert(i_send < n_send);
        std::vector<M>& to_send_v = to_send[id_node];
        long offset = 0;
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
  f.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xl_s = mod(xl + shift_remain, total_site);
    const long index_s = geo.index_from_coordinate(xl_s);
    Vector<M> fv = f.get_elems(index_s);
    const Vector<long> frv = f_recv_idx.get_elems_const(index);
    const int id_node = frv[0];
    const long offset = frv[1] * geo.multiplicity;
    std::vector<M>& to_recv_v = to_recv[id_node];
    for (int m = 0; m < geo.multiplicity; ++m) {
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

template <class M>
void qnorm_field(FieldM<double, 1>& f, const Field<M>& f1)
{
  TIMER("qnorm_field");
  const Geometry geo = geo_reform(f1.geo());
  f.init();
  f.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> f1v = f1.get_elems_const(xl);
    f.get_elem(index) = qnorm(f1v);
  });
}

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                         \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator+=                       \
      <TYPENAME>(Field<TYPENAME>& f, const Field<TYPENAME>& f1);               \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator-=                       \
      <TYPENAME>(Field<TYPENAME>& f, const Field<TYPENAME>& f1);               \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                      \
      Field<TYPENAME>& f, const Field<double>& f_factor);                      \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=(                      \
      Field<TYPENAME>& f, const Field<Complex>& f_factor);                     \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                       \
      <TYPENAME>(Field<TYPENAME>& f, const double factor);                     \
                                                                               \
  QLAT_EXTERN template const Field<TYPENAME>& operator*=                       \
      <TYPENAME>(Field<TYPENAME>& f, const Complex& factor);                   \
                                                                               \
  QLAT_EXTERN template double qnorm<TYPENAME>(const Field<TYPENAME>& f);       \
                                                                               \
  QLAT_EXTERN template double qnorm_double<TYPENAME>(                          \
      const Field<TYPENAME>& f1, const Field<TYPENAME>& f2);                   \
                                                                               \
  QLAT_EXTERN template void qswap<TYPENAME>(Field<TYPENAME> & f1,              \
                                            Field<TYPENAME> & f2);             \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum<TYPENAME>(              \
      const Field<TYPENAME>& f);                                               \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_sum_tslice<TYPENAME>(       \
      const Field<TYPENAME>& f, const int t_dir);                              \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum_double<TYPENAME>(   \
      const Field<TYPENAME>& f);                                               \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_glb_sum_long<TYPENAME>(     \
      const Field<TYPENAME>& f);                                               \
                                                                               \
  QLAT_EXTERN template std::vector<std::vector<TYPENAME> >                     \
  field_glb_sum_tslice_double<TYPENAME>(const Field<TYPENAME>& f,              \
                                        const int t_dir);                      \
                                                                               \
  QLAT_EXTERN template std::vector<std::vector<TYPENAME> >                     \
  field_glb_sum_tslice_long<TYPENAME>(const Field<TYPENAME>& f,                \
                                      const int t_dir);                        \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_project_mom<TYPENAME>(      \
      const Field<TYPENAME>& f, const CoordinateD& mom);                       \
                                                                               \
  QLAT_EXTERN template std::vector<TYPENAME> field_get_elems<TYPENAME>(        \
      const Field<TYPENAME>& f, const Coordinate& xg);                         \
                                                                               \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                      \
      const Field<TYPENAME>& f, const Coordinate& xg, const int m);            \
                                                                               \
  QLAT_EXTERN template TYPENAME field_get_elem<TYPENAME>(                      \
      const Field<TYPENAME>& f, const Coordinate& xg);                         \
                                                                               \
  QLAT_EXTERN template void field_set_elems<TYPENAME>(                         \
      Field<TYPENAME> & f, const Coordinate& xg, const Vector<TYPENAME> val);  \
                                                                               \
  QLAT_EXTERN template void field_set_elem<TYPENAME>(                          \
      Field<TYPENAME> & f, const Coordinate& xg, const int m,                  \
      const TYPENAME& val);                                                    \
                                                                               \
  QLAT_EXTERN template void field_set_elem<TYPENAME>(                          \
      Field<TYPENAME> & f, const Coordinate& xg, const TYPENAME& val);         \
                                                                               \
  QLAT_EXTERN template void split_fields<TYPENAME>(                            \
      std::vector<Handle<Field<TYPENAME> > > & vec, const Field<TYPENAME>& f); \
                                                                               \
  QLAT_EXTERN template void merge_fields<TYPENAME>(                            \
      Field<TYPENAME> & f,                                                     \
      const std::vector<ConstHandle<Field<TYPENAME> > >& vec);                 \
                                                                               \
  QLAT_EXTERN template void merge_fields_ms<TYPENAME>(                         \
      Field<TYPENAME> & f,                                                     \
      const std::vector<ConstHandle<Field<TYPENAME> > >& vec,                  \
      const std::vector<int> m_vec);                                           \
                                                                               \
  QLAT_EXTERN template void field_shift_dir<TYPENAME>(                         \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1, const int dir,           \
      const int shift);                                                        \
                                                                               \
  QLAT_EXTERN template void field_shift_steps<TYPENAME>(                       \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1,                          \
      const Coordinate& shift);                                                \
                                                                               \
  QLAT_EXTERN template void field_shift_direct<TYPENAME>(                      \
      Field<TYPENAME> & f, const Field<TYPENAME>& f1,                          \
      const Coordinate& shift);                                                \
                                                                               \
  QLAT_EXTERN template void field_shift<TYPENAME>(Field<TYPENAME> & f,         \
                                                  const Field<TYPENAME>& f1,   \
                                                  const Coordinate& shift);    \
                                                                               \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(FieldM<double, 1> & f,       \
                                                  const Field<TYPENAME>& f1)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);

#define QLAT_EXTERN_TEMPLATE_2(TYPENAME1, TYPENAME2)    \
  QLAT_EXTERN template void assign(Field<TYPENAME1>& f, \
                                   const Field<TYPENAME2>& f1)

QLAT_CALL_WITH_TYPES_2(QLAT_EXTERN_TEMPLATE_2);

#undef QLAT_EXTERN
#undef QLAT_EXTERN_TEMPLATE
#undef QLAT_EXTERN_CLASS

}  // namespace qlat
