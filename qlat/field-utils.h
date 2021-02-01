#pragma once

#include <qlat/config.h>
#include <qlat/coordinate-d.h>
#include <qlat/field-shuffle.h>
#include <qlat/field.h>
#include <qlat/matrix.h>

namespace qlat
{  //

// field.h and field-utils.h are including each other. Need this forward
// declaration. template <class M> struct Field; End of forward declaration.

template <class M>
void set_zero(Field<M>& f)
{
  TIMER("set_zero(Field)");
  set_zero(f.field);
}

template <class M>
void set_unit(Field<M>& f, const Complex& coef = 1.0)
{
  for (long offset = 0; offset < f.field.size(); ++offset) {
    set_unit(f.get_elem(offset), coef);
  }
}

template <class M>
std::vector<M> field_sum(const Field<M>& f)
{
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
std::vector<M> field_glb_sum_double(const Field<M>& f)
{
  std::vector<M> vec = field_sum(f);
  glb_sum_double_vec(Vector<M>(vec));
  return vec;
}

template <class M>
std::vector<M> field_glb_sum_long(const Field<M>& f)
{
  std::vector<M> vec = field_sum(f);
  glb_sum_long_vec(Vector<M>(vec));
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
  std::vector<M> ret(geo.multiplicity);
  set_zero(ret);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0;
    for (int k = 0; k < DIMN; ++k) {
      phase += mom[k] * xg[k];
    }
    const Vector<M> v = f.get_elems_const(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      ret[m] += std::polar(1.0, -phase) * v[m];
    }
  }
  glb_sum_double_vec(get_data(ret));
  return ret;
}

template <class M>
std::vector<M> field_get_elems(const Field<M>& f, const Coordinate& xg)
{
  const Geometry& geo = f.geo();
  const Coordinate xl = geo.coordinate_l_from_g(xg);
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
{
  return field_get_elems(f, xg)[m];
}

template <class M>
M field_get_elem(const Field<M>& f, const Coordinate& xg)
{
  const Geometry& geo = f.geo();
  qassert(geo.multiplicity == 1);
  return field_get_elem(f, xg, 0);
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

inline long& get_max_field_shift_direct_msg_size()
{
  static long size = 1024L * 1024L * 1024L;
  // static long size = 16;
  return size;
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
  qassert(geo.is_only_local());
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
    std::vector<MPI_Request> send_reqs(n_send), recv_reqs(n_recv);
    const int mpi_tag = 11;
    int i_send = 0;
    int i_recv = 0;
    for (int i = 0; i < num_node; ++i) {
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
          recv_reqs.push_back(req);
          offset += max_elem;
          size_s -= max_elem;
        }
        MPI_Isend(&to_send_v[offset], size_s * sizeof(M), MPI_BYTE, id_node,
                  mpi_tag, get_comm(), &send_reqs[i_send]);
        i_send += 1;
      }
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
          recv_reqs.push_back(req);
          offset += max_elem;
          size_r -= max_elem;
        }
        MPI_Irecv(&to_recv_v[offset], size_r * sizeof(M), MPI_BYTE, id_node,
                  mpi_tag, get_comm(), &recv_reqs[i_recv]);
        i_recv += 1;
      }
    }
    MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
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
void set_u_rand_double(Field<M>& f, const RngState& rs,
                       const double upper = 1.0, const double lower = -1.0)
{
  TIMER("set_u_rand_double");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = f.get_elems(xl);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, 1.0, -1.0);
    }
  }
}

template <class M>
void set_u_rand_float(Field<M>& f, const RngState& rs, const double upper = 1.0,
                      const double lower = -1.0)
{
  TIMER("set_u_rand_float");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = f.get_elems(xl);
    Vector<float> dv((float*)v.data(), v.data_size() / sizeof(float));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, 1.0, -1.0);
    }
  }
}

}  // namespace qlat
