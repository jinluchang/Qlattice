#include <qlat/selected-shuffle.h>

namespace qlat
{  //

void SelectedShufflePlan::init()
{
  points_dist_type_send = PointsDistType::Local;
  points_dist_type_recv = PointsDistType::Random;
  num_selected_points_send = 0;
  num_selected_points_recv = 0;
  n_points_selected_points_send.clear();
  n_points_selected_points_recv.clear();
  n_points_selected_points_send.set_mem_type(MemType::Cpu);
  n_points_selected_points_recv.set_mem_type(MemType::Cpu);
  SelectedPoints<Long>& spi_s = shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = shuffle_idx_points_recv;
  SelectedPoints<Long>& spi_l = shuffle_idx_points_local;
  spi_s.init();
  spi_r.init();
  spi_l.init();
  spi_s.set_mem_type(MemType::Cpu);
  spi_r.set_mem_type(MemType::Cpu);
  spi_l.set_mem_type(MemType::Cpu);
  total_count_send = 0;
  total_count_recv = 0;
  total_count_local = 0;
  sendcounts.clear();
  recvcounts.clear();
  sdispls.clear();
  rdispls.clear();
  sendcounts.set_mem_type(MemType::Cpu);
  recvcounts.set_mem_type(MemType::Cpu);
  sdispls.set_mem_type(MemType::Cpu);
  rdispls.set_mem_type(MemType::Cpu);
}

void shuffle_selected_points_char(
    std::vector<SelectedPoints<Char>>& spc_vec,
    const std::vector<SelectedPoints<Char>>& spc0_vec,
    const SelectedShufflePlan& ssp)
// spc0_vec is the source to be send.
// spc_vec is the dest to be filled by received data.
// The size of the `std::vector`s and the size of `SelectedPoints`s (both send
// and recv) should be set (matching `ssp`) before calling this function.
// All `multiplicity` should be the same.
{
  TIMER_FLOPS("shuffle_selected_points_char(spc_vec,spc0_vec,ssp)");
  qassert(ssp.num_selected_points_send == (Long)spc0_vec.size());
  qassert(ssp.num_selected_points_recv == (Long)spc_vec.size());
  qassert(ssp.num_selected_points_send > 0);
  const Int multiplicity = spc0_vec[0].multiplicity;
  for (Int i = 0; i < ssp.num_selected_points_send; ++i) {
    qassert(ssp.n_points_selected_points_send[i] == spc0_vec[i].n_points);
    qassert(spc0_vec[i].initialized == true);
    qassert(spc0_vec[i].points_dist_type == ssp.points_dist_type_send);
    qassert(spc0_vec[i].multiplicity == multiplicity);
  }
  for (Int i = 0; i < ssp.num_selected_points_recv; ++i) {
    qassert(ssp.n_points_selected_points_recv[i] == spc_vec[i].n_points);
    qassert(spc_vec[i].initialized == true);
    qassert(spc_vec[i].points_dist_type == ssp.points_dist_type_recv);
    qassert(spc_vec[i].multiplicity == multiplicity);
  }
  const SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  const SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  const SelectedPoints<Long>& spi_l = ssp.shuffle_idx_points_local;
  qassert(spi_s.initialized);
  qassert(spi_s.multiplicity == 3);
  qassert(spi_s.n_points == ssp.total_count_send);
  qassert(spi_r.initialized);
  qassert(spi_r.multiplicity == 3);
  qassert(spi_r.n_points == ssp.total_count_recv);
  qassert(spi_l.initialized);
  qassert(spi_l.multiplicity == 4);
  qassert(spi_l.n_points == ssp.total_count_local);
  qthread_for(idx, spi_l.n_points, {
    const Vector<Long> v = spi_l.get_elems_const(idx);
    const Long idx_selected_points_send = v[0];
    const Long idx_within_field_send = v[1];
    const Long idx_selected_points_recv = v[2];
    const Long idx_within_field_recv = v[3];
    qassert(0 <= idx_selected_points_send and
            idx_selected_points_send < ssp.num_selected_points_send);
    qassert(0 <= idx_selected_points_recv and
            idx_selected_points_recv < ssp.num_selected_points_recv);
    const Vector<Char> v_val =
        spc0_vec[idx_selected_points_send].get_elems_const(
            idx_within_field_send);
    Vector<Char> v1_val =
        spc_vec[idx_selected_points_recv].get_elems(idx_within_field_recv);
    assign(v1_val, v_val);
  });
  // Initialized `sp` to be the target of the shuffle before final shuffle.
  SelectedPoints<Char> sp;
  sp.set_mem_type(MemType::Comm);
  sp.init(ssp.total_count_recv, multiplicity, PointsDistType::Local);
  // Copy spc0 to sp0. Reordered to be ready to send.
  SelectedPoints<Char> sp0;
  sp0.set_mem_type(MemType::Comm);
  sp0.init(ssp.total_count_send, multiplicity, PointsDistType::Local);
  qthread_for(idx, spi_s.n_points, {
    const Vector<Long> v = spi_s.get_elems_const(idx);
    const Long idx_selected_points_send = v[0];
    const Long idx_within_send_field = v[1];
    const Long idx_send_buffer = v[2];
    qassert(0 <= idx_selected_points_send and
            idx_selected_points_send < ssp.num_selected_points_send);
    const Vector<Char> v_val =
        spc0_vec[idx_selected_points_send].get_elems_const(
            idx_within_send_field);
    Vector<Char> v1_val = sp0.get_elems(idx_send_buffer);
    assign(v1_val, v_val);
  });
  // Perform shuffle with `MPI_Alltoallv`.
  {
    TIMER_FLOPS("shuffle_selected_points_char(spc,spc0,ssp)-mpi");
    const MpiDataType& mpi_dtype = get_mpi_data_type_contiguous(multiplicity);
    {
      TIMER_FLOPS("shuffle_selected_points_char(spc,spc0,ssp)-MPI_Alltoallv");
      mpi_alltoallv(sp0.points.data(), ssp.sendcounts.data(),
                    ssp.sdispls.data(), mpi_dtype.mpi_dtype, sp.points.data(),
                    ssp.recvcounts.data(), ssp.rdispls.data(),
                    mpi_dtype.mpi_dtype, get_comm());
      timer.flops += (ssp.total_count_send + ssp.total_count_recv) / 2;
    }
    sp0.init();
    timer.flops += (ssp.total_count_send + ssp.total_count_recv) / 2;
  }
  // perform final reordering.
  qthread_for(idx, spi_r.n_points, {
    const Vector<Long> v = spi_r.get_elems_const(idx);
    const Long idx_selected_points_recv = v[0];
    const Long idx_within_field_recv = v[1];
    const Long idx_buffer_recv = v[2];
    qassert(0 <= idx_selected_points_recv and
            idx_selected_points_recv < ssp.num_selected_points_recv);
    const Vector<Char> v_val = sp.get_elems_const(idx_buffer_recv);
    Vector<Char> v1_val =
        spc_vec[idx_selected_points_recv].get_elems(idx_within_field_recv);
    assign(v1_val, v_val);
  });
  timer.flops +=
      (ssp.total_count_send + ssp.total_count_recv) / 2 + ssp.total_count_local;
}

// ------------------------------


// ------------------------------

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

static void set_selected_shuffle_plan_no_reorder(
    SelectedShufflePlan& ssp, const SelectedPoints<Int>& sp_id_node_send_to)
// Collective operation.
// partially set SelectedShufflePlan
//
// ssp.points_dist_type_recv = PointsDistType::Random.
// ssp.shuffle_idx_points_recv is trivially set (no reorder at all).
//
// To be set in `set_selected_shuffle_plan(ssp,psel,rs)`.
{
  TIMER("set_selected_shuffle_plan_no_reorder(ssp,spinst)");
  ssp.init();
  const Int num_node = get_num_node();
  const Int id_node_local = get_id_node();
  qassert(sp_id_node_send_to.initialized == true);
  qassert(sp_id_node_send_to.points_dist_type != PointsDistType::Global);
  qassert(sp_id_node_send_to.multiplicity == 1);
  ssp.num_selected_points_send = 1;
  ssp.n_points_selected_points_send.resize(ssp.num_selected_points_send);
  set_zero(ssp.n_points_selected_points_send);
  ssp.points_dist_type_send = sp_id_node_send_to.points_dist_type;
  SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  SelectedPoints<Long>& spi_l = ssp.shuffle_idx_points_local;
  ssp.sdispls.resize(num_node);
  ssp.rdispls.resize(num_node);
  ssp.sendcounts.set_mem_type(MemType::Comm);
  ssp.recvcounts.set_mem_type(MemType::Comm);
  ssp.sendcounts.resize(num_node);
  ssp.recvcounts.resize(num_node);
  set_zero(ssp.sdispls);
  set_zero(ssp.rdispls);
  set_zero(ssp.sendcounts);
  set_zero(ssp.recvcounts);
  ssp.n_points_selected_points_send[0] = sp_id_node_send_to.n_points;
  qfor(idx, ssp.n_points_selected_points_send[0], {
    const Int id_node_send_to = sp_id_node_send_to.get_elem(idx);
    if (id_node_send_to == id_node_local) {
      ssp.total_count_local += 1;
    } else {
      ssp.total_count_send += 1;
      ssp.sendcounts[id_node_send_to] += 1;
    }
  });
  Long sdispl = 0;
  qfor(id_node, num_node, {
    ssp.sdispls[id_node] = sdispl;
    sdispl += ssp.sendcounts[id_node];
  });
  qassert(ssp.total_count_send == sdispl);
  Long idx_local = 0;
  Long idx_send = 0;
  Long c_idx_local = 0;
  vector<Long> c_idx_vec(MemType::Cpu);
  c_idx_vec = ssp.sdispls;
  spi_l.init(ssp.total_count_local, 4, PointsDistType::Local);
  spi_s.init(ssp.total_count_send, 3, PointsDistType::Local);
  set_zero(spi_l);
  set_zero(spi_s);
  qfor(idx, ssp.n_points_selected_points_send[0], {
    const Int id_node_send_to = sp_id_node_send_to.get_elem(idx);
    if (id_node_send_to == id_node_local) {
      Vector<Long> v = spi_l.get_elems(idx_local);
      v[0] = 0;
      v[1] = idx;
      v[2] = 0;
      v[3] = c_idx_local;
      c_idx_local += 1;
      idx_local += 1;
    } else {
      Vector<Long> v = spi_s.get_elems(idx_send);
      v[0] = 0;
      v[1] = idx;
      v[2] = c_idx_vec[id_node_send_to];
      c_idx_vec[id_node_send_to] += 1;
      idx_send += 1;
    }
  });
  qassert(idx_local == ssp.total_count_local);
  qassert(idx_send == ssp.total_count_send);
  qassert(c_idx_local == ssp.total_count_local);
  qfor(id_node, num_node - 1,
       { qassert(c_idx_vec[id_node] == ssp.sdispls[id_node + 1]); });
  qassert(c_idx_vec[num_node - 1] == ssp.total_count_send);
  MPI_Alltoall(ssp.sendcounts.data(), sizeof(Long), MPI_BYTE,
               ssp.recvcounts.data(), sizeof(Long), MPI_BYTE, get_comm());
  Long rdispl = 0;
  qfor(id_node, num_node, {
    ssp.rdispls[id_node] = rdispl;
    rdispl += ssp.recvcounts[id_node];
  });
  ssp.total_count_recv = rdispl;
  ssp.points_dist_type_recv = PointsDistType::Random;
  spi_r.init(ssp.total_count_recv, 3, PointsDistType::Local);
  qthread_for(idx, spi_r.n_points, {
    Vector<Long> v = spi_r.get_elems(idx);
    v[0] = 0;
    v[1] = ssp.total_count_local + idx;
    v[2] = idx;
  });
  ssp.sendcounts.set_mem_type(MemType::Cpu);
  ssp.recvcounts.set_mem_type(MemType::Cpu);
  ssp.num_selected_points_recv = 1;
  ssp.n_points_selected_points_recv.resize(ssp.num_selected_points_recv);
  set_zero(ssp.n_points_selected_points_recv);
  ssp.n_points_selected_points_recv[0] =
      ssp.total_count_local + ssp.total_count_recv;
}

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const PointsSelection& psel, const RngState& rs)
// Collective operation.
// make shuffle plan
// psel.points_dist_type == PointsDistType::Local
// ssp.points_dist_type_recv = PointsDistType::Random
// Sort the shuffled points by order of the gindex of points.
//
// Call `set_selected_shuffle_plan_no_reorder(ssp,spinst)` to set ssp.
// In addition, set the missing:
// `ssp.points_dist_type_recv`
// `ssp.shuffle_idx_points_recv`
{
  TIMER("set_selected_shuffle_plan(ssp,psel,rs)");
  const Long n_points = psel.size();
  SelectedPoints<Int> sp_id_node_send_to;
  set_selected_shuffle_id_node_send_to(sp_id_node_send_to, psel, rs);
  set_selected_shuffle_plan_no_reorder(ssp, sp_id_node_send_to);
  ssp.points_dist_type_recv = PointsDistType::Random;
  SelectedPoints<Long> sp_idx0;
  sp_idx0.set_mem_type(MemType::Cpu);
  sp_idx0.init(n_points, 1, ssp.points_dist_type_send);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Long gindex = index_from_coordinate(xg, psel.total_site);
    sp_idx0.get_elem(idx) = gindex;
  });
  SelectedPoints<Long> sp_idx;
  sp_idx.set_mem_type(MemType::Cpu);
  shuffle_selected_points(sp_idx, sp_idx0, ssp);
  std::vector<std::pair<Long, Long>> idx_pair_vec(sp_idx.n_points);
  // Set the `pair.first` to be the rank that determine the target order.
  // Let `pair.second` remembers the initial order.
  qthread_for(idx, sp_idx.n_points, {
    idx_pair_vec[idx].first = sp_idx.get_elem(idx);
    idx_pair_vec[idx].second = idx;
  });
  // Sort according to rank
  std::sort(idx_pair_vec.begin(), idx_pair_vec.end());
  // Set the pair.second to be the target idx
  qthread_for(idx, sp_idx.n_points, {
    idx_pair_vec[idx].first = idx_pair_vec[idx].second;
    idx_pair_vec[idx].second = idx;
  });
  // Reorder back to the original order
  std::sort(idx_pair_vec.begin(), idx_pair_vec.end());
  SelectedPoints<Long>& spi_l = ssp.shuffle_idx_points_local;
  qthread_for(idx, spi_l.n_points, {
    Vector<Long> v = spi_l.get_elems(idx);
    const Long idx_init = v[3];
    qassert(idx_init == idx_pair_vec[idx_init].first);
    const Long idx_tgt = idx_pair_vec[idx_init].second;
    v[3] = idx_tgt;
  });
  SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  qthread_for(idx, spi_r.n_points, {
    Vector<Long> v = spi_r.get_elems(idx);
    const Long idx_init = v[1];
    qassert(idx_init == idx_pair_vec[idx_init].first);
    const Long idx_tgt = idx_pair_vec[idx_init].second;
    v[1] = idx_tgt;
  });
}

// ------------------------------

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp)
// const Long n_points = ssp.total_count_recv;
// const Int multiplicity = sp0.multiplicity;
// SelectedPoints<M> sp;
// sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
// SelectedPoints<Char> spc(sp.view_as_char());
// SelectedPoints<Char> spc0(sp0.view_as_char());
{
  TIMER_FLOPS("shuffle_selected_points_char(spc,spc0,ssp)");
  qassert(ssp.num_selected_points_send == 1);
  qassert(ssp.num_selected_points_recv == 1);
  std::vector<SelectedPoints<Char>> spc_vec(1);
  std::vector<SelectedPoints<Char>> spc0_vec(1);
  spc_vec[0].set_view(spc);
  spc0_vec[0].set_view(spc0);
  shuffle_selected_points_char(spc_vec, spc0_vec, ssp);
  timer.flops +=
      (ssp.total_count_send + ssp.total_count_recv) / 2 + ssp.total_count_local;
}

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection(sp,psel,psel0,ssp)");
  qassert(psel0.points_dist_type == ssp.points_dist_type_send);
  qassert(ssp.num_selected_points_send == 1);
  qassert(ssp.num_selected_points_recv == 1);
  const Long n_points = ssp.n_points_selected_points_recv[0];
  psel.init(psel0.total_site, n_points);
  psel.points_dist_type = ssp.points_dist_type_recv;
  SelectedPoints<Char> pselc(psel.view_sp().view_as_char());
  const SelectedPoints<Char> pselc0(psel0.view_sp().view_as_char());
  shuffle_selected_points_char(pselc, pselc0, ssp);
}

// ------------------------------

void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0,
                             const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_field_selection(psel,fsel0,ssp)");
  PointsSelection psel0;
  set_psel_from_fsel(psel0, fsel0);
  shuffle_points_selection(psel, psel0, ssp);
}

}  // namespace qlat
