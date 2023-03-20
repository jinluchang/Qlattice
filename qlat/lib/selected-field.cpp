#define QLAT_INSTANTIATE_SELECTED_FIELD

#include <qlat/selected-field.h>

namespace qlat
{  //

void add_field_selection(FieldM<int64_t, 1>& f_rank, const PointSelection& psel,
                         const long rank_psel)
// interface function
// add psel points to f_rank. (only lower rank if already selected)
{
  TIMER_VERBOSE("add_field_selection(psel)");
  const Geometry& geo = f_rank.geo();
#pragma omp parallel for
  for (long i = 0; i < (long)psel.size(); ++i) {
    const Coordinate xl = geo.coordinate_l_from_g(psel[i]);
    if (geo.is_local(xl)) {
      int64_t& rank = f_rank.get_elem(xl);
      if (rank < 0 or rank > rank_psel) {
        rank = rank_psel;
      }
    }
  }
}

void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                        const Coordinate& total_site, const int64_t val)
// interface function
// select everything with val
// default val = 0 ; means selection everything
// val = -1 deselection everything
{
  TIMER_VERBOSE("mk_field_selection(f_rank,total_site)");
  Geometry geo;
  geo.init(total_site, 1);
  f_rank.init();
  f_rank.init(geo);
  qassert(f_rank.geo().is_only_local);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    int64_t& rank = f_rank.get_elem(index);
    rank = val;
  }
}

void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                        const Coordinate& total_site,
                        const std::vector<Coordinate>& xgs, const long rank_xgs)
{
  TIMER_VERBOSE("mk_field_selection(xgs)");
  mk_field_selection(f_rank, total_site, -1);
  add_field_selection(f_rank, xgs, rank_xgs);
}

void select_rank_range(FieldM<int64_t, 1>& f_rank, const long rank_start,
                       const long rank_stop)
// keep rank info if rank_start <= rank and (rank < rank_stop or rank_stop ==
// -1) otherwise rank = -1 default parameter does not change selection but will
// erase the rank information for points not selected (rank = -1)
{
  TIMER_VERBOSE("select_rank_range");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  // const Coordinate total_site = geo.total_site();
  qacc_for(index, geo.local_volume(), {
    int64_t& rank = f_rank.get_elem(index);
    if (not(rank_start <= rank and (rank < rank_stop or rank_stop == -1))) {
      rank = -1;
    }
  });
}

void select_t_range(FieldM<int64_t, 1>& f_rank, const long t_start,
                    const long t_stop)
// keep rank info if t_start <= t and (t < t_stop or t_stop == -1)
// otherwise rank = -1
// default parameter does not change selection
// but will erase the rank information for points not selected (rank = -1)
{
  TIMER_VERBOSE("select_t_range");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  // const Coordinate total_site = geo.total_site();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int t = xg[3];
    if (not(t_start <= t and (t < t_stop or t_stop == -1))) {
      int64_t& rank = f_rank.get_elem(index);
      rank = -1;
    }
  });
}

void set_n_per_tslice(FieldM<int64_t, 1>& f_rank, const long n_per_tslice)
// will erase the rank information for points not selected (rank = -1)
//
// if n_per_tslice == -1 then all points are selected regardless of rank
// (if a point was not selected before (rank < 0), then rank will be set to be
// rank = spatial_vol).
//
// otherwise: n_per_tslice is not enforced but only serve as an limit for f_rank
//
// 0 <= rank < n_per_tslice
//
// if point is not selected, rank = -1
{
  TIMER_VERBOSE("set_n_per_tslice");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  const Coordinate total_site = geo.total_site();
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  qassert(n_per_tslice == -1 or
          (0 <= n_per_tslice and n_per_tslice <= spatial_vol));
  qacc_for(index, geo.local_volume(), {
    int64_t& rank = f_rank.get_elem(index);
    if (n_per_tslice == -1 and rank < 0) {
      rank = spatial_vol;
    } else if (not(0 <= rank and rank < n_per_tslice)) {
      rank = -1;
    }
  });
}

void update_field_selection(FieldSelection& fsel)
// interface function
// update fsel based only on f_rank
// do not touch n_per_tslice and prob at all
{
  TIMER_VERBOSE("update_field_selection");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  fsel.f_local_idx.init();
  fsel.f_local_idx.init(geo);
  long n_elems = 0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const int64_t& rank = fsel.f_rank.get_elem(index);
    long& idx = fsel.f_local_idx.get_elem(index);
    if (0 <= rank) {
      idx = n_elems;
      n_elems += 1;
    } else {
      idx = -1;
    }
  }
  fsel.n_elems = n_elems;
  fsel.ranks.resize(fsel.n_elems);
  fsel.indices.resize(fsel.n_elems);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const long idx = fsel.f_local_idx.get_elems_const(index)[0];
    if (idx >= 0) {
      const long rank = fsel.f_rank.get_elem(index);
      fsel.ranks[idx] = rank;
      fsel.indices[idx] = index;
    }
  }
}

void update_field_selection(FieldSelection& fsel, const long n_per_tslice_)
// interface function
// only adjust parameter, do not change contents
{
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  const Coordinate total_site = geo.total_site();
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  qassert(n_per_tslice_ == -1 or
          (0 <= n_per_tslice_ and n_per_tslice_ <= spatial_vol));
  fsel.n_per_tslice = n_per_tslice_ == -1 ? spatial_vol : n_per_tslice_;
  fsel.prob = (double)fsel.n_per_tslice / (double)spatial_vol;
}

void set_field_selection(FieldSelection& fsel, const FieldM<int64_t, 1>& f_rank,
                         const long n_per_tslice_, const bool is_limit_on_rank)
// call set_n_per_tslice if is_limit_on_rank = true
// otherwise will strictly follow f_rank without constraint of n_per_tslice
{
  TIMER_VERBOSE("set_field_selection(fsel,f_rank,n_per_tslice)");
  fsel.init();
  fsel.f_rank = f_rank;
  if (is_limit_on_rank) {
    set_n_per_tslice(fsel.f_rank, n_per_tslice_);
  }
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice_);
}

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site)
// select everything with rank = 0
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site);
  update_field_selection(fsel);
  update_field_selection(fsel, -1);  // select all points
}

bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2)
// only check selection, does not check rank or parameter
{
  const long n_elems = fsel1.n_elems;
  if (n_elems != fsel2.n_elems) {
    return false;
  }
  bool is_same = true;
  const vector_acc<long>& indices1 = fsel1.indices;
  const vector_acc<long>& indices2 = fsel2.indices;
  qassert(indices1.size() == n_elems);
  qassert(indices2.size() == n_elems);
  qthread_for(idx, n_elems, {
    if (indices1[idx] != indices2[idx]) {
      is_same = false;
    }
  });
  return is_same;
}

PointSelection psel_from_fsel(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel")
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  long n_elems = fsel.n_elems;
  long total_n_elems = n_elems;
  glb_sum(total_n_elems);
  // const int num_node = geo.geon.num_node;
  const int id_node = geo.geon.id_node;
  vector<long> vec(geo.geon.num_node, 0);
  all_gather(get_data(vec), get_data_one_elem(n_elems));
  long idx_offset = 0;
  for (int i = 0; i < id_node; ++i) {
    idx_offset += vec[i];
  }
  qassert(idx_offset <= total_n_elems);
  vector<long> vec_gindex(total_n_elems, 0);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = index_from_coordinate(xg, total_site);
    vec_gindex[idx_offset + idx] = gindex;
  });
  glb_sum(get_data(vec_gindex));
  PointSelection psel(total_n_elems);
  qthread_for(idx, (long)psel.size(), {
    long gindex = vec_gindex[idx];
    psel[idx] = coordinate_from_index(gindex, total_site);
  });
  return psel;
}

PointSelection psel_from_fsel_local(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel_local")
  const Geometry& geo = fsel.f_rank.geo();
  // const Coordinate total_site = geo.total_site();
  long n_elems = fsel.n_elems;
  PointSelection psel(n_elems);
  qthread_for(idx, (long)psel.size(), {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    psel[idx] = xg;
  });
  return psel;
}

void set_selected_gindex(SelectedField<long>& sfgi, const FieldSelection& fsel)
{
  TIMER_VERBOSE("set_selected_gindex");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  sfgi.init(fsel, 1);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = index_from_coordinate(xg, total_site);
    sfgi.get_elem(idx) = gindex;
  });
}

}  // namespace qlat
