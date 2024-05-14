#include <qlat/selected-shuffle.h>

namespace qlat
{  //

void SelectedShufflePlan::init()
{
  points_dist_type_send = PointsDistType::Local;
  points_dist_type_recv = PointsDistType::Random;
  SelectedPoints<Long>& spi_s = shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = shuffle_idx_points_recv;
  spi_s.init();
  spi_r.init();
  total_send_count = 0;
  total_recv_count = 0;
  sdispls.clear();
  rdispls.clear();
  sendcounts.clear();
  recvcounts.clear();
}

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sp_id_node_send_to, const Long n_points,
    const RngState& rs)
// assume `points_dist_type_send` is `PointsDistType::Local`.
{
  TIMER("set_selected_shuffle_id_node_send_to(spist,n_points,rs)");
  const Int num_node = get_num_node();
  sp_id_node_send_to.init(n_points, 1, PointsDistType::Local);
  RngState rsl = rs.split(get_id_node());
  qthread_for(idx, n_points, {
    RngState rsi = rsl.newtype(idx);
    const Int id_node_send_to = rand_gen(rsi) % num_node;
    qassert(0 <= id_node_send_to);
    qassert(id_node_send_to < num_node);
    sp_id_node_send_to.get_elem(idx) = id_node_send_to;
  });
}

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sp_id_node_send_to, const PointsSelection& psel,
    const RngState& rs)
{
  TIMER("set_selected_shuffle_id_node_send_to(spist,psel,rs)");
  const Long n_points = psel.size();
  const Int num_node = get_num_node();
  sp_id_node_send_to.init(n_points, 1, psel.points_dist_type);
  RngState rsl = rs.split(get_id_node());
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Long gindex = index_from_coordinate(xg, psel.total_site);
    RngState rsi = rsl.newtype(gindex);
    const Int id_node_send_to = rand_gen(rsi) % num_node;
    qassert(0 <= id_node_send_to);
    qassert(id_node_send_to < num_node);
    sp_id_node_send_to.get_elem(idx) = id_node_send_to;
  });
}

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const SelectedPoints<Int>& sp_id_node_send_to)
// Collective operation.
// partially set SelectedShufflePlan
//
// Missing elements:
// ssp.points_dist_type_recv
// ssp.shuffle_idx_points_recv
//
// To be set in `set_selected_shuffle_plan(ssp,psel,rs)`.
{
  TIMER("set_selected_shuffle_plan(ssp,spinst)");
  ssp.init();
  const Int num_node = get_num_node();
  const Long n_points = sp_id_node_send_to.n_points;
  qassert(sp_id_node_send_to.initialized == true);
  qassert(sp_id_node_send_to.points_dist_type != PointsDistType::Global);
  qassert(sp_id_node_send_to.multiplicity == 1);
  ssp.points_dist_type_send = sp_id_node_send_to.points_dist_type;
  SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  spi_s.init(n_points, 1, ssp.points_dist_type_send);
  set_zero(spi_s);
  ssp.total_send_count = n_points;
  ssp.sdispls.resize(num_node);
  ssp.rdispls.resize(num_node);
  ssp.sendcounts.resize(num_node);
  ssp.recvcounts.resize(num_node);
  set_zero(ssp.sdispls);
  set_zero(ssp.rdispls);
  set_zero(ssp.sendcounts);
  set_zero(ssp.recvcounts);
  qfor(idx, n_points, {
    const Int id_node_send_to = sp_id_node_send_to.get_elem(idx);
    ssp.sendcounts[id_node_send_to] += 1;
  });
  Int sdispl = 0;
  qfor(id_node, num_node, {
    ssp.sdispls[id_node] = sdispl;
    sdispl += ssp.sendcounts[id_node];
  });
  qassert(ssp.total_send_count == sdispl);
  vector<Int> c_idx_vec;
  c_idx_vec = ssp.sdispls;
  qfor(idx, n_points, {
    const Int id_node_send_to = sp_id_node_send_to.get_elem(idx);
    spi_s.get_elem(idx) = c_idx_vec[id_node_send_to];
    c_idx_vec[id_node_send_to] += 1;
  });
  qfor(id_node, num_node - 1,
       { qassert(c_idx_vec[id_node] == ssp.sdispls[id_node + 1]); });
  MPI_Alltoall(ssp.sendcounts.data(), sizeof(Int), MPI_BYTE,
               ssp.recvcounts.data(), sizeof(Int), MPI_BYTE, get_comm());
  Int rdispl = 0;
  qfor(id_node, num_node, {
    ssp.rdispls[id_node] = rdispl;
    rdispl += ssp.recvcounts[id_node];
  });
  ssp.total_recv_count = rdispl;
}

void set_selected_shuffle_plan(SelectedShufflePlan& ssp, const Long n_points,
                               const RngState& rs)
{
  TIMER("set_selected_shuffle_plan(ssp,n_points,rs)");
  SelectedPoints<Int> sp_id_node_send_to;
  set_selected_shuffle_id_node_send_to(sp_id_node_send_to, n_points, rs);
  set_selected_shuffle_plan(ssp, sp_id_node_send_to);
}

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const PointsSelection& psel, const RngState& rs)
// Call `set_selected_shuffle_plan(ssp,spinst)` to set ssp.
// In addition, set the missing:
// `ssp.points_dist_type_recv`
// `ssp.shuffle_idx_points_recv`
{
  TIMER("set_selected_shuffle_plan(ssp,psel,rs)");
  const Long n_points = psel.size();
  SelectedPoints<Int> sp_id_node_send_to;
  set_selected_shuffle_id_node_send_to(sp_id_node_send_to, psel, rs);
  set_selected_shuffle_plan(ssp, sp_id_node_send_to);
  ssp.points_dist_type_recv = PointsDistType::Random;
  SelectedPoints<Long> sp_idx0;
  sp_idx0.init(n_points, 1, ssp.points_dist_type_send);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Long gindex = index_from_coordinate(xg, psel.total_site);
    sp_idx0.get_elem(idx) = gindex;
  });
  SelectedPoints<Long> sp_idx;
  shuffle_selected_points(sp_idx, sp_idx0, ssp);
  qassert(sp_idx.n_points == ssp.total_recv_count);
  std::vector<std::pair<Long, Long>> idx_pair_vec(sp_idx.n_points);
  qthread_for(idx, sp_idx.n_points, {
    idx_pair_vec[idx].first = sp_idx.get_elem(idx);
    idx_pair_vec[idx].second = idx;
  });
  std::sort(idx_pair_vec.begin(), idx_pair_vec.end());
  SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  spi_r.init(sp_idx.n_points, 1, ssp.points_dist_type_recv);
  qthread_for(idx, sp_idx.n_points, {
    const Long idx_src = idx_pair_vec[idx].second;
    const Long idx_tgt = idx;
    spi_r.get_elem(idx_src) = idx_tgt;
  });
}

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp)
// const Long n_points = ssp.total_recv_count;
// const Int multiplicity = sp0.multiplicity;
// SelectedPoints<M> sp;
// sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
// SelectedPoints<Char> spc(sp.view_as_char());
// SelectedPoints<Char> spc0(sp0.view_as_char());
{
  TIMER("shuffle_selected_points_char(spc,spc0,ssp)");
  const SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  const SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  qassert(spi_s.initialized);
  qassert(spi_s.multiplicity == 1);
  qassert(spi_s.n_points == ssp.total_send_count);
  qassert(spc0.initialized == true);
  qassert(spc0.points_dist_type == ssp.points_dist_type_send);
  qassert(spc0.n_points == ssp.total_send_count);
  const Int multiplicity = spc0.multiplicity;
  qassert(spc.initialized == true);
  qassert(spc.points_dist_type == ssp.points_dist_type_recv);
  qassert(spc.n_points == ssp.total_recv_count);
  qassert(spc.multiplicity == multiplicity);
  SelectedPoints<Char> sp;
  if (spi_r.initialized) {
    qassert(spi_r.multiplicity == 1);
    qassert(spi_r.n_points == ssp.total_recv_count);
    sp.init(spc.n_points, multiplicity, ssp.points_dist_type_recv);
  } else {
    sp.set_view(spc);
  }
  SelectedPoints<Char> sp0;
  sp0.init(spc0.n_points, multiplicity, ssp.points_dist_type_send);
  qthread_for(idx, spc0.n_points, {
    const Vector<Char> v = spc0.get_elems_const(idx);
    const Long idx1 = spi_s.get_elem(idx);
    qassert(0 <= idx1);
    qassert(idx1 < sp0.n_points);
    Vector<Char> v1 = sp0.get_elems(idx1);
    assign(v1, v);
  });
  MPI_Datatype mpi_dtype;
  const Int mpi_ret = MPI_Type_contiguous(multiplicity, MPI_BYTE, &mpi_dtype);
  qassert(mpi_ret == 0);
  MPI_Type_commit(&mpi_dtype);
  MPI_Alltoallv(sp0.points.data(), ssp.sendcounts.data(), ssp.sdispls.data(),
                mpi_dtype, sp.points.data(), ssp.recvcounts.data(),
                ssp.rdispls.data(), mpi_dtype, get_comm());
  MPI_Type_free(&mpi_dtype);
  if (spi_r.initialized) {
    qthread_for(idx, sp.n_points, {
      const Vector<Char> v = sp.get_elems_const(idx);
      const Long idx1 = spi_r.get_elem(idx);
      qassert(0 <= idx1);
      qassert(idx1 < spc.n_points);
      Vector<Char> v1 = spc.get_elems(idx1);
      assign(v1, v);
    });
  }
}

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection(sp,psel,psel0,ssp)");
  qassert(psel0.points_dist_type == ssp.points_dist_type_send);
  const Long n_points = ssp.total_recv_count;
  psel.init(psel0.total_site, n_points);
  psel.points_dist_type = ssp.points_dist_type_recv;
  SelectedPoints<Char> pselc(psel.view_sp().view_as_char());
  const SelectedPoints<Char> pselc0(psel0.view_sp().view_as_char());
  shuffle_selected_points_char(pselc, pselc0, ssp);
}

void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0,
                             const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_field_selection(psel,fsel0,ssp)");
  PointsSelection psel0;
  set_psel_from_fsel(psel0, fsel0);
  shuffle_points_selection(psel, psel0, ssp);
}

}  // namespace qlat
