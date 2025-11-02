#define QLAT_INSTANTIATE_SELECTED_FIELD

#include <qlat/selected-field-io.h>

namespace qlat
{  //

void mk_field_selection(FieldRank& f_rank, const Geometry& geo_,
                        const int64_t val)
// interface function
// select everything with val
// default val = 0 ; means selection everything
// val = -1 deselection everything
{
  TIMER_VERBOSE("mk_field_selection(f_rank,geo,val)");
  const Geometry geo = geo_resize(geo_);
  f_rank.init();
  f_rank.init(geo, 1);
  qassert(f_rank.geo().is_only_local);
  qacc_for(index, geo.local_volume(), {
    int64_t& rank = f_rank.get_elem(index);
    rank = val;
  });
}

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const int64_t val)
// interface function
// select everything with val
// default val = 0 ; means selection everything
// val = -1 deselection everything
{
  TIMER_VERBOSE("mk_field_selection(f_rank,total_site,val)");
  Geometry geo;
  geo.init(total_site);
  mk_field_selection(f_rank, geo, val);
}

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const PointsSelection& psel, const Long rank_xgs)
{
  TIMER_VERBOSE("mk_field_selection(psel)");
  mk_field_selection(f_rank, total_site, -1);
  add_field_selection(f_rank, psel, rank_xgs);
}

void set_selected_gindex(SelectedField<Long>& sfgi, const FieldSelection& fsel)
{
  TIMER_VERBOSE("set_selected_gindex");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  sfgi.init(fsel, 1);
  qthread_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = index_from_coordinate(xg, total_site);
    sfgi.get_elem(idx) = gindex;
  });
}

void mk_grid_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                             const Long n_per_tslice_, const RngState& rs)
// each time slice has "n_per_tslice = spatial_vol / ratio" points been
// selected. not selected points has value = -1; selected points has value from
// 0 to "n_per_tslice - 1" in random order (different per time slice)
{
  TIMER_VERBOSE("mk_grid_field_selection");
  const Long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  const Long n_per_tslice = n_per_tslice_ == -1 ? spatial_vol : n_per_tslice_;
  qassert(n_per_tslice > 0);
  qassert(spatial_vol % n_per_tslice == 0);
  const Long ratio = spatial_vol / n_per_tslice;
  qassert(ratio >= 0);
  RngState rs_shift(rs, "random_shift");
  const Coordinate random_shift = mod(
      Coordinate(rand_gen(rs_shift), rand_gen(rs_shift), rand_gen(rs_shift), 0),
      total_site);
  Geometry geo;
  geo.init(total_site);
  f_rank.init();
  f_rank.init(geo, 1);
  qassert(f_rank.geo().is_only_local);
  qthread_for(index, geo.local_volume(), { f_rank.get_elem(index) = -1; });
  std::vector<Field<int64_t>> fs;
  const Coordinate new_size_node = get_default_serial_new_size_node(geo);
  shuffle_field(fs, f_rank, new_size_node);
  qassert(fs.size() <= 1);
  if (fs.size() == 1) {
    // require each tslice is on one node.
    Field<int64_t>& nf = fs[0];
    const Geometry& ngeo = nf.geo();
    qassert(new_size_node == ngeo.geon.size_node);
    const Int t_start = ngeo.node_site[3] * ngeo.geon.coor_node[3];
    const Int t_end = t_start + ngeo.node_site[3];
#pragma omp parallel for
    for (Int t = t_start; t < t_end; ++t) {
      RngState rst = rs.split(t);
      std::vector<Long> ranks(n_per_tslice);
      for (Long i = 0; i < n_per_tslice; ++i) {
        ranks[i] = i;
      }
      random_permute(ranks, rst);
      Long idx = 0;
      for (Long index = 0; index < ngeo.local_volume(); ++index) {
        const Coordinate xl = ngeo.coordinate_from_index(index);
        const Coordinate xg = ngeo.coordinate_g_from_l(xl);
        if (xg[3] != t) {
          continue;
        }
        const Coordinate x = mod(xg + random_shift, total_site);
        bool check = true;
        switch (ratio) {
          case 1:
            check = true;
            break;
          case 2:
            check = check and (x[0] + x[1] + x[2]) % 2 == 0;
            break;
          case 4:
            check = check and (x[0] + x[1]) % 2 == 0;
            check = check and (x[1] + x[2]) % 2 == 0;
            check = check and (x[0] + x[2]) % 2 == 0;
            break;
          case 8:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            break;
          case 16:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            check = check and (x[0] + x[1] + x[2]) % 4 == 0;
            break;
          case 32:
            check = check and x[0] % 2 == 0;
            check = check and x[1] % 2 == 0;
            check = check and x[2] % 2 == 0;
            check = check and (x[0] + x[1]) % 4 == 0;
            check = check and (x[1] + x[2]) % 4 == 0;
            check = check and (x[0] + x[2]) % 4 == 0;
            break;
          case 64:
            check = check and x[0] % 4 == 0;
            check = check and x[1] % 4 == 0;
            check = check and x[2] % 4 == 0;
            break;
          default:
            displayln_info(fname + ssprintf(": ERROR: ratio=%d.", ratio));
            qassert(false);
            break;
        }
        if (check) {
          qassert(idx < n_per_tslice);
          nf.get_elem(index) = ranks[idx];
          idx += 1;
        } else {
          nf.get_elem(index) = -1;
        }
      }
      qassert(idx == n_per_tslice);
    }
  }
  shuffle_field_back(f_rank, fs, new_size_node);
}

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const Long n_per_tslice, const RngState& rs)
// interface function
// not selected points has value = -1;
// random select n_per_tslice points based on ranks from mk_grid_field_selection
{
  TIMER_VERBOSE("mk_field_selection(n_per_tslice,rs)");
  const Long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  mk_grid_field_selection(f_rank, total_site, spatial_vol, rs);
  if ((n_per_tslice == -1) or (n_per_tslice == spatial_vol)) {
    return;
  }
  set_n_per_tslice(f_rank, n_per_tslice);
}

void add_field_selection(FieldRank& f_rank, const PointsSelection& psel,
                         const Long rank_psel)
// interface function
// f_rank needs to be already initialized
// add psel points to f_rank. (only lower rank if already selected)
{
  TIMER_VERBOSE("add_field_selection(f_rank,psel,rank_psel)");
  const Geometry& geo = f_rank.geo();
  qthread_for(i, (Long)psel.size(), {
    const Coordinate xl = geo.coordinate_l_from_g(psel[i]);
    if (geo.is_local(xl)) {
      int64_t& rank = f_rank.get_elem(xl);
      if ((rank < 0) or (rank > rank_psel)) {
        rank = rank_psel;
      }
    }
  });
}

void add_field_selection(FieldRank& f_rank, const FieldSelection& fsel)
// interface function
// f_rank needs to be already initialized
// add fsel points to f_rank. (only lower rank if already selected)
{
  TIMER_VERBOSE("add_field_selection(f_rank,fsel)");
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const int64_t rank_fsel = fsel.ranks[idx];
    int64_t& rank = f_rank.get_elem(index);
    if ((rank < 0) or (rank > rank_fsel)) {
      rank = rank_fsel;
    }
  });
}

void select_rank_range(FieldRank& f_rank, const Long rank_start,
                       const Long rank_stop)
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
    if (not((rank_start <= rank) and
            ((rank < rank_stop) or (rank_stop == -1)))) {
      rank = -1;
    }
  });
}

void select_t_range(FieldRank& f_rank, const Long t_start, const Long t_stop)
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
    const Int t = xg[3];
    if (not((t_start <= t) and ((t < t_stop) or (t_stop == -1)))) {
      int64_t& rank = f_rank.get_elem(index);
      rank = -1;
    }
  });
}

void set_n_per_tslice(FieldRank& f_rank, const Long n_per_tslice)
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
  const Long spatial_vol = total_site[0] * total_site[1] * total_site[2];
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
{
  TIMER("update_field_selection");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  fsel.f_local_idx.init();
  fsel.f_local_idx.init(geo);
  Long n_elems = 0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const int64_t& rank = fsel.f_rank.get_elem(index);
    Long& idx = fsel.f_local_idx.get_elem(index);
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
  qthread_for(index, geo.local_volume(), {
    const Long idx = fsel.f_local_idx.get_elem(index);
    if (idx >= 0) {
      const Long rank = fsel.f_rank.get_elem(index);
      fsel.ranks[idx] = rank;
      fsel.indices[idx] = index;
    }
  });
}

void set_grid_field_selection(FieldSelection& fsel,
                              const Coordinate& total_site,
                              const Long n_per_tslice, const RngState& rs)
{
  TIMER_VERBOSE("set_grid_field_selection(fsel,total_site,n_per_tslice,rs)");
  fsel.init();
  mk_grid_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  update_field_selection(fsel);
}

void set_field_selection(FieldSelection& fsel, const FieldRank& f_rank)
// strictly follow f_rank on selection
{
  TIMER("set_field_selection(fsel,f_rank)");
  fsel.init();
  fsel.f_rank = f_rank;
  update_field_selection(fsel);
}

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site)
// select everything with rank = 0
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site);
  update_field_selection(fsel);  // select all points
}

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const Long n_per_tslice, const RngState& rs)
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site,n_per_tslice,rs)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  update_field_selection(fsel);
}

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const Long n_per_tslice, const RngState& rs,
                         const PointsSelection& psel)
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site,n_per_tslice,rs,psel)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
  add_field_selection(fsel.f_rank, psel);
  update_field_selection(fsel);
}

bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2)
// only check selection, does not check rank or parameter
{
  const Long n_elems = fsel1.n_elems;
  if (n_elems != fsel2.n_elems) {
    return false;
  }
  bool is_same = true;
  const vector<Long>& indices1 = fsel1.indices;
  const vector<Long>& indices2 = fsel2.indices;
  qassert(indices1.size() == n_elems);
  qassert(indices2.size() == n_elems);
  qthread_for(idx, n_elems, {
    if (indices1[idx] != indices2[idx]) {
      is_same = false;
    }
  });
  return is_same;
}

bool is_containing(const PointsSelection& psel,
                   const FieldSelection& fsel_small)
// local checking
{
  TIMER("is_containing(psel,fsel_small)");
  const Geometry& geo = fsel_small.f_rank.geo();
  const Long n_points0 = psel.size();
  Long n_missing_points = 0;
  Long idx_last = -1;
  qfor(idx, fsel_small.indices.size(), {
    const Long index = fsel_small.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    bool is_found = false;
    for (Long i = 0; i < n_points0; ++i) {
      idx_last += 1;
      if (idx_last >= n_points0) {
        idx_last = idx_last % n_points0;
      }
      const Coordinate& xg0 = psel[idx_last];
      if (xg0 == xg) {
        is_found = true;
        break;
      }
    }
    if (not is_found) {
      n_missing_points += 1;
      break;
    }
  });
  glb_sum(n_missing_points);
  if (n_missing_points == 0) {
    return true;
  } else {
    qassert(n_missing_points > 0);
    return false;
  }
}

bool is_containing(const PointsSelection& psel,
                   const PointsSelection& psel_small)
// local checking
{
  TIMER("is_containing(psel,psel_small)");
  if (psel.points_dist_type != psel_small.points_dist_type) {
    return false;
  }
  const Long n_points = psel_small.size();
  const Long n_points0 = psel.size();
  Long n_missing_points = 0;
  Long idx_last = -1;
  qfor(idx, n_points, {
    const Coordinate& xg = psel_small[idx];
    bool is_found = false;
    for (Long i = 0; i < n_points0; ++i) {
      idx_last += 1;
      if (idx_last >= n_points0) {
        idx_last = idx_last % n_points0;
      }
      const Coordinate& xg0 = psel[idx_last];
      if (xg0 == xg) {
        is_found = true;
        break;
      }
    }
    if (not is_found) {
      n_missing_points += 1;
      break;
    }
  });
  glb_sum(n_missing_points);
  if (n_missing_points == 0) {
    return true;
  } else {
    qassert(n_missing_points > 0);
    return false;
  }
}

bool is_containing(const FieldSelection& fsel, const FieldSelection& fsel_small)
// local checking
{
  TIMER("is_containing(fsel,fsel_small)");
  Long n_missing_points = 0;
  qthread_for(idx, fsel_small.indices.size(), {
    const Long index = fsel_small.indices[idx];
    const int64_t rank = fsel.f_rank.get_elem(index);
    if (rank < 0) {
      // May miscount due to race condition.
      // But we only need to know if it is zero or not.
      // Should be fine.
      n_missing_points += 1;
    }
  });
  return n_missing_points == 0;
}

bool is_containing(const FieldSelection& fsel,
                   const PointsSelection& psel_small)
// local checking
{
  TIMER("is_containing(fsel,psel_small)");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  Long n_missing_points = 0;
  qthread_for(i, (Long)psel_small.size(), {
    const Coordinate xl = geo.coordinate_l_from_g(psel_small[i]);
    if (geo.is_local(xl)) {
      const int64_t rank = fsel.f_rank.get_elem(xl);
      if (rank < 0) {
        // May miscount due to race condition.
        // But we only need to know if it is zero or not.
        // Should be fine.
        n_missing_points += 1;
      }
    }
  });
  return n_missing_points == 0;
}

void intersect_with(FieldSelection& fsel, const FieldSelection& fsel1)
{
  TIMER("intersect_with(fsel,fsel1)");
  const Geometry& geo = fsel.f_rank.geo();
  Qassert(geo == fsel1.f_rank.geo());
  qthread_for(idx, fsel.indices.size(), {
    const Long index = fsel.indices[idx];
    const int64_t rank = fsel1.f_rank.get_elem(index);
    if (rank < 0) {
      fsel.f_rank.get_elem(index) = -1;
    }
  });
  update_field_selection(fsel);
}

PointsSelection intersect(const FieldSelection& fsel,
                          const PointsSelection& psel)
{
  TIMER("intersect(fsel,psel)");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  qassert(geo.total_site() == psel.total_site);
  qassert(psel.points_dist_type == PointsDistType::Global);
  vector<Int> is_psel_in_fsel(psel.size());
  set_zero(is_psel_in_fsel);
  qthread_for(i, (Long)psel.size(), {
    const Coordinate xl = geo.coordinate_l_from_g(psel[i]);
    if (geo.is_local(xl)) {
      const int64_t rank = fsel.f_rank.get_elem(xl);
      if (rank < 0) {
        is_psel_in_fsel[i] = 1;
      }
    }
  });
  glb_sum(is_psel_in_fsel);
  Long n_points = 0;
  qfor(i, is_psel_in_fsel.size(), {
    if (is_psel_in_fsel[i] == 0) {
      n_points += 1;
    }
  });
  PointsSelection psel_new(psel.total_site, n_points);
  n_points = 0;
  qfor(i, is_psel_in_fsel.size(), {
    if (is_psel_in_fsel[i] == 0) {
      qassert(n_points < (Long)psel_new.size());
      psel_new[n_points] = psel[i];
      n_points += 1;
    }
  });
  return psel_new;
}

PointsSelection psel_from_fsel(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel")
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const Long n_elems = fsel.n_elems;
  Long total_n_elems = n_elems;
  glb_sum(total_n_elems);
  // const Int num_node = geo.geon.num_node;
  const Int id_node = geo.geon.id_node;
  vector<Long> vec(geo.geon.num_node);
  set_zero(vec);
  all_gather(get_data(vec), get_data_one_elem(n_elems));
  Long idx_offset = 0;
  for (Int i = 0; i < id_node; ++i) {
    idx_offset += vec[i];
  }
  qassert(idx_offset <= total_n_elems);
  vector<Long> vec_gindex(total_n_elems);
  set_zero(vec_gindex);
  qthread_for(idx, n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = index_from_coordinate(xg, total_site);
    vec_gindex[idx_offset + idx] = gindex;
  });
  glb_sum(get_data(vec_gindex));
  PointsSelection psel(total_site, total_n_elems);
  qthread_for(idx, (Long)psel.size(), {
    Long gindex = vec_gindex[idx];
    psel[idx] = coordinate_from_index(gindex, total_site);
  });
  return psel;
}

PointsSelection psel_from_fsel_local(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel_local")
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const Long n_elems = fsel.n_elems;
  PointsSelection psel(total_site, n_elems);
  psel.points_dist_type = PointsDistType::Local;
  qthread_for(idx, (Long)psel.size(), {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    psel[idx] = xg;
  });
  return psel;
}

Long write_field_selection(const FieldSelection& fsel, const std::string& path)
{
  TIMER_VERBOSE("write_field_selection");
  return write_field_64(fsel.f_rank, path);
}

Long read_field_selection(FieldSelection& fsel, const std::string& path)
{
  TIMER_VERBOSE("read_field_selection");
  fsel.init();
  FieldRank f_rank;
  const Long total_bytes = read_field_64(f_rank, path);
  if (total_bytes > 0) {
    set_field_selection(fsel, f_rank);
  }
  return total_bytes;
}

Long idx_from_xg(const Coordinate& xg, const FieldSelection& fsel)
// idx for xg for this node.
// idx >= 0 if in the local_volume and selected in fsel.
// idx=-2 if not in the local volume.
// idx=-1 if not selected in fsel.
{
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  if (not geo.is_local(xl)) {
    return -2;
  }
  const Long index = geo.index_from_coordinate(xl);
  const Long idx = fsel.f_local_idx.get_elem(index);
  return idx;
}

std::string make_selected_field_header(const Geometry& geo,
                                       const Int multiplicity,
                                       const Int sizeof_M, const crc32_t crc32)
{
  const Coordinate total_site = geo.total_site();
  std::ostringstream out;
  // const std::string todo = "NOT yet implemented";
  out << "BEGIN_SELECTED_FIELD_HEADER" << std::endl;
  out << "selected_field_version = 1.0" << std::endl;
  out << "total_site[0] = " << total_site[0] << std::endl;
  out << "total_site[1] = " << total_site[1] << std::endl;
  out << "total_site[2] = " << total_site[2] << std::endl;
  out << "total_site[3] = " << total_site[3] << std::endl;
  out << "multiplicity = " << multiplicity << std::endl;
  out << "sizeof(M) = " << sizeof_M << std::endl;
  out << ssprintf("selected_field_crc32 = %08X", crc32) << std::endl;
  out << "END_HEADER" << std::endl;
  return out.str();
}

Long read_selected_geo_info(Coordinate& total_site, int& multiplicity,
                            int& sizeof_M, crc32_t& crc,
                            const std::string& path)
{
  TIMER("read_selected_geo_info");
  Long pos = 0;
  if (get_id_node() == 0) {
    QFile fp = qfopen(path, "r");
    if (not fp.null()) {
      const std::string header = "BEGIN_SELECTED_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          std::vector<std::string> infos;
          infos.push_back(header);
          while (infos.back() != "END_HEADER\n" && infos.back() != "") {
            infos.push_back(qgetline(fp));
          }
          for (Int m = 0; m < 4; ++m) {
            reads(total_site[m],
                  info_get_prop(infos, ssprintf("total_site[%d] = ", m)));
          }
          reads(multiplicity, info_get_prop(infos, "multiplicity = "));
          reads(sizeof_M, info_get_prop(infos, "sizeof(M) = "));
          crc = read_crc32(info_get_prop(infos, "selected_field_crc32 = "));
        }
      }
      pos = qftell(fp);
    }
    qfclose(fp);
  }
  bcast(get_data_one_elem(pos));
  bcast(get_data_one_elem(total_site));
  bcast(get_data_one_elem(multiplicity));
  bcast(get_data_one_elem(sizeof_M));
  bcast(get_data_one_elem(crc));
  return pos;
}

bool is_selected_field(const std::string& path)
{
  TIMER("is_selected_field");
  Long nfile = 0;
  if (get_id_node() == 0) {
    QFile fp = qfopen(path, "r");
    if (not fp.null()) {
      const std::string header = "BEGIN_SELECTED_FIELD_HEADER\n";
      std::vector<char> check_line(header.size(), 0);
      if (1 == qfread(check_line.data(), header.size(), 1, fp)) {
        if (std::string(check_line.data(), check_line.size()) == header) {
          nfile = 1;
        }
      }
    }
    qfclose(fp);
  }
  bcast(get_data(nfile));
  return nfile > 0;
}

void set_sqrt_field(SelectedField<RealD>& f, const SelectedField<RealD>& f1)
{
  TIMER("set_sqrt_field(f,f1)");
  const Geometry& geo = f1.geo();
  const Int multiplicity = f1.multiplicity;
  f.init(geo, f1.n_elems, multiplicity);
  qacc_for(idx, f.n_elems, {
    const Vector<RealD> f1v = f1.get_elems_const(idx);
    Vector<RealD> fv = f.get_elems(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      fv[m] = std::sqrt(f1v[m]);
    }
  });
}

}  // namespace qlat
