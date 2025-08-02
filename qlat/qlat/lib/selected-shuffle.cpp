#include <qlat/selected-shuffle.h>

#include <array>

namespace qlat
{  //

void SelectedShufflePlan::init()
{
  points_dist_type_send = PointsDistType::Local;
  points_dist_type_recv = PointsDistType::Random;
  total_site = Coordinate();
  size_node_send.clear();
  coor_node_send.clear();
  size_node_recv.clear();
  coor_node_recv.clear();
  num_selected_points_send = 0;
  num_selected_points_recv = 0;
  n_points_selected_points_send.clear();
  n_points_selected_points_recv.clear();
  n_points_selected_points_send.set_mem_type(MemType::Comm);
  n_points_selected_points_recv.set_mem_type(MemType::Comm);
  SelectedPoints<Long>& spi_s = shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = shuffle_idx_points_recv;
  SelectedPoints<Long>& spi_l = shuffle_idx_points_local;
  spi_s.init();
  spi_r.init();
  spi_l.init();
  spi_s.set_mem_type(MemType::Comm);
  spi_r.set_mem_type(MemType::Comm);
  spi_l.set_mem_type(MemType::Comm);
  total_count_send = 0;
  total_count_recv = 0;
  total_count_local = 0;
  sendcounts.clear();
  recvcounts.clear();
  sdispls.clear();
  rdispls.clear();
  sendcounts.set_mem_type(MemType::Comm);
  recvcounts.set_mem_type(MemType::Comm);
  sdispls.set_mem_type(MemType::Comm);
  rdispls.set_mem_type(MemType::Comm);
}

// ------------------------------

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
  const Int num_selected_points_send = ssp.num_selected_points_send;
  const Int num_selected_points_recv = ssp.num_selected_points_recv;
  qassert(num_selected_points_send == (Long)spc0_vec.size());
  qassert(f_glb_sum(num_selected_points_send) > 0);
  const bool b_has_send = num_selected_points_send > 0;
  const Int multiplicity =
      f_bcast_any(b_has_send ? spc0_vec[0].multiplicity : 0, b_has_send);
  vector<SelectedPoints<Char>> view0_vec(num_selected_points_send,
                                         MemType::Comm);
  vector<SelectedPoints<Char>> view_vec(num_selected_points_recv,
                                        MemType::Comm);
  set_zero(view0_vec);
  set_zero(view_vec);
  for (Int i = 0; i < num_selected_points_send; ++i) {
    qassert(ssp.n_points_selected_points_send[i] == spc0_vec[i].n_points);
    qassert(spc0_vec[i].initialized == true);
    qassert(spc0_vec[i].points_dist_type == ssp.points_dist_type_send);
    qassert(spc0_vec[i].multiplicity == multiplicity);
    view0_vec[i].set_view(spc0_vec[i]);
  }
  spc_vec.resize(num_selected_points_recv);
  for (Int i = 0; i < num_selected_points_recv; ++i) {
    spc_vec[i].init(ssp.n_points_selected_points_recv[i], multiplicity,
                    ssp.points_dist_type_recv);
    view_vec[i].set_view(spc_vec[i]);
  }
  view_vec.set_mem_type(MemType::CommAcc);
  view0_vec.set_mem_type(MemType::CommAcc);
  SelectedPoints<Long> spi_s, spi_r, spi_l;
  spi_s.set_view(ssp.shuffle_idx_points_send);
  spi_r.set_view(ssp.shuffle_idx_points_recv);
  spi_l.set_view(ssp.shuffle_idx_points_local);
  qassert(spi_s.initialized);
  qassert(spi_s.multiplicity == 3);
  qassert(spi_s.n_points == ssp.total_count_send);
  qassert(spi_s.points.mem_type == MemType::CommAcc);
  qassert(spi_r.initialized);
  qassert(spi_r.multiplicity == 3);
  qassert(spi_r.n_points == ssp.total_count_recv);
  qassert(spi_r.points.mem_type == MemType::CommAcc);
  qassert(spi_l.initialized);
  qassert(spi_l.multiplicity == 4);
  qassert(spi_l.n_points == ssp.total_count_local);
  qassert(spi_l.points.mem_type == MemType::CommAcc);
  {
    TIMER_FLOPS("shuffle_selected_points_char-local");
    timer.flops += spi_l.n_points * multiplicity;
    qmem_for(idx, spi_l.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_l.get_elems_const(idx);
      const Long idx_selected_points_send = v[0];
      const Long idx_within_field_send = v[1];
      const Long idx_selected_points_recv = v[2];
      const Long idx_within_field_recv = v[3];
      qassert(0 <= idx_selected_points_send and
              idx_selected_points_send < num_selected_points_send);
      qassert(0 <= idx_selected_points_recv and
              idx_selected_points_recv < num_selected_points_recv);
      const Vector<Char> v_val =
          view0_vec[idx_selected_points_send].get_elems_const(
              idx_within_field_send);
      Vector<Char> v1_val =
          view_vec[idx_selected_points_recv].get_elems(idx_within_field_recv);
      assign(v1_val, v_val);
    });
  }
  // Initialized `sp` to be the target of the shuffle before final shuffle.
  SelectedPoints<Char> sp;
  sp.set_mem_type(MemType::CommAcc);
  sp.init(ssp.total_count_recv, multiplicity, PointsDistType::Local);
  // Copy spc0 to sp0. Reordered to be ready to send.
  SelectedPoints<Char> sp0;
  sp0.set_mem_type(MemType::CommAcc);
  sp0.init(ssp.total_count_send, multiplicity, PointsDistType::Local);
  {
    TIMER_FLOPS("shuffle_selected_points_char-send");
    timer.flops += spi_s.n_points * multiplicity;
    qmem_for(idx, spi_s.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_s.get_elems_const(idx);
      const Long idx_selected_points_send = v[0];
      const Long idx_within_field_send = v[1];
      const Long idx_buffer_send = v[2];
      qassert(0 <= idx_selected_points_send and
              idx_selected_points_send < num_selected_points_send);
      const Vector<Char> v_val =
          view0_vec[idx_selected_points_send].get_elems_const(
              idx_within_field_send);
      Vector<Char> v1_val = sp0.get_elems(idx_buffer_send);
      assign(v1_val, v_val);
    });
  }
  // Perform shuffle with `mpi_alltoallv`.
  {
    TIMER_FLOPS("shuffle_selected_points_char-mpi");
    const MpiDataType& mpi_dtype = get_mpi_data_type_contiguous(multiplicity);
    {
      TIMER_FLOPS("shuffle_selected_points_char-mpi_alltoallv");
      mpi_alltoallv(sp0.points.data(), ssp.sendcounts.data(),
                    ssp.sdispls.data(), mpi_dtype.mpi_dtype, sp.points.data(),
                    ssp.recvcounts.data(), ssp.rdispls.data(),
                    mpi_dtype.mpi_dtype, get_comm());
      timer.flops +=
          (ssp.total_count_send + ssp.total_count_recv) / 2 * multiplicity;
    }
    sp0.init();
    timer.flops +=
        (ssp.total_count_send + ssp.total_count_recv) / 2 * multiplicity;
  }
  {
    TIMER_FLOPS("shuffle_selected_points_char-recv");
    timer.flops += spi_r.n_points * multiplicity;
    // perform final reordering.
    qmem_for(idx, spi_r.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_r.get_elems_const(idx);
      const Long idx_selected_points_recv = v[0];
      const Long idx_within_field_recv = v[1];
      const Long idx_buffer_recv = v[2];
      qassert(0 <= idx_selected_points_recv and
              idx_selected_points_recv < num_selected_points_recv);
      const Vector<Char> v_val = sp.get_elems_const(idx_buffer_recv);
      Vector<Char> v1_val =
          view_vec[idx_selected_points_recv].get_elems(idx_within_field_recv);
      assign(v1_val, v_val);
    });
  }
  timer.flops += ((ssp.total_count_send + ssp.total_count_recv) / 2 +
                  ssp.total_count_local) *
                 multiplicity;
}

void shuffle_selected_points_back_char(
    std::vector<SelectedPoints<Char>>& spc_vec,
    const std::vector<SelectedPoints<Char>>& spc0_vec,
    const SelectedShufflePlan& ssp)
// spc0_vec is the source to be send.
// spc_vec is the dest to be filled by received data.
// The size of the `std::vector`s and the size of `SelectedPoints`s (both send
// and recv) should be set (matching `ssp`) before calling this function.
// All `multiplicity` should be the same.
// The shuffle direction is opposite to `shuffle_selected_points_char` (or the
// direction suggested by `ssp`).
{
  TIMER_FLOPS("shuffle_selected_points_back_char(spc_vec,spc0_vec,ssp)");
  const Int num_selected_points_send = ssp.num_selected_points_send;
  const Int num_selected_points_recv = ssp.num_selected_points_recv;
  qassert(num_selected_points_recv == (Long)spc0_vec.size());
  qassert(f_glb_sum(num_selected_points_recv) > 0);
  const bool b_has_send = num_selected_points_recv > 0;
  const Int multiplicity =
      f_bcast_any(b_has_send ? spc0_vec[0].multiplicity : 0, b_has_send);
  vector<SelectedPoints<Char>> view0_vec(num_selected_points_recv,
                                         MemType::Comm);
  vector<SelectedPoints<Char>> view_vec(num_selected_points_send,
                                        MemType::Comm);
  set_zero(view0_vec);
  set_zero(view_vec);
  for (Int i = 0; i < num_selected_points_recv; ++i) {
    qassert(ssp.n_points_selected_points_recv[i] == spc0_vec[i].n_points);
    qassert(spc0_vec[i].initialized == true);
    qassert(spc0_vec[i].points_dist_type == ssp.points_dist_type_recv);
    qassert(spc0_vec[i].multiplicity == multiplicity);
    view0_vec[i].set_view(spc0_vec[i]);
  }
  spc_vec.resize(num_selected_points_send);
  for (Int i = 0; i < num_selected_points_send; ++i) {
    spc_vec[i].init(ssp.n_points_selected_points_send[i], multiplicity,
                    ssp.points_dist_type_send);
    view_vec[i].set_view(spc_vec[i]);
  }
  view_vec.set_mem_type(MemType::CommAcc);
  view0_vec.set_mem_type(MemType::CommAcc);
  SelectedPoints<Long> spi_s, spi_r, spi_l;
  spi_s.set_view(ssp.shuffle_idx_points_send);
  spi_r.set_view(ssp.shuffle_idx_points_recv);
  spi_l.set_view(ssp.shuffle_idx_points_local);
  qassert(spi_s.initialized);
  qassert(spi_s.multiplicity == 3);
  qassert(spi_s.n_points == ssp.total_count_send);
  qassert(spi_r.initialized);
  qassert(spi_r.multiplicity == 3);
  qassert(spi_r.n_points == ssp.total_count_recv);
  qassert(spi_l.initialized);
  qassert(spi_l.multiplicity == 4);
  qassert(spi_l.n_points == ssp.total_count_local);
  {
    TIMER_FLOPS("shuffle_selected_points_back_char-local");
    timer.flops += spi_l.n_points * multiplicity;
    qmem_for(idx, spi_l.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_l.get_elems_const(idx);
      const Long idx_selected_points_send = v[0];
      const Long idx_within_field_send = v[1];
      const Long idx_selected_points_recv = v[2];
      const Long idx_within_field_recv = v[3];
      qassert(0 <= idx_selected_points_send and
              idx_selected_points_send < num_selected_points_send);
      qassert(0 <= idx_selected_points_recv and
              idx_selected_points_recv < num_selected_points_recv);
      const Vector<Char> v_val =
          view0_vec[idx_selected_points_recv].get_elems_const(
              idx_within_field_recv);
      Vector<Char> v1_val =
          view_vec[idx_selected_points_send].get_elems(idx_within_field_send);
      assign(v1_val, v_val);
    });
  }
  // Initialized `sp` to be the target of the shuffle before final shuffle.
  SelectedPoints<Char> sp;
  sp.set_mem_type(MemType::CommAcc);
  sp.init(ssp.total_count_send, multiplicity, PointsDistType::Local);
  // Copy spc0 to sp0. Reordered to be ready to send.
  SelectedPoints<Char> sp0;
  sp0.set_mem_type(MemType::CommAcc);
  sp0.init(ssp.total_count_recv, multiplicity, PointsDistType::Local);
  {
    TIMER_FLOPS("shuffle_selected_points_back_char-send");
    timer.flops += spi_r.n_points * multiplicity;
    qmem_for(idx, spi_r.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_r.get_elems_const(idx);
      const Long idx_selected_points_recv = v[0];
      const Long idx_within_field_recv = v[1];
      const Long idx_buffer_recv = v[2];
      qassert(0 <= idx_selected_points_recv and
              idx_selected_points_recv < num_selected_points_recv);
      const Vector<Char> v_val =
          view0_vec[idx_selected_points_recv].get_elems_const(
              idx_within_field_recv);
      Vector<Char> v1_val = sp0.get_elems(idx_buffer_recv);
      assign(v1_val, v_val);
    });
  }
  // Perform shuffle with `mpi_alltoallv`.
  {
    TIMER_FLOPS("shuffle_selected_points_back_char-mpi");
    const MpiDataType& mpi_dtype = get_mpi_data_type_contiguous(multiplicity);
    {
      TIMER_FLOPS("shuffle_selected_points_back_char-mpi_alltoallv");
      mpi_alltoallv(sp0.points.data(), ssp.recvcounts.data(),
                    ssp.rdispls.data(), mpi_dtype.mpi_dtype, sp.points.data(),
                    ssp.sendcounts.data(), ssp.sdispls.data(),
                    mpi_dtype.mpi_dtype, get_comm());
      timer.flops +=
          (ssp.total_count_send + ssp.total_count_recv) / 2 * multiplicity;
    }
    sp0.init();
    timer.flops +=
        (ssp.total_count_send + ssp.total_count_recv) / 2 * multiplicity;
  }
  {
    TIMER_FLOPS("shuffle_selected_points_back_char-recv");
    timer.flops += spi_s.n_points * multiplicity;
    // perform final reordering.
    qmem_for(idx, spi_s.n_points, MemType::CommAcc, {
      const Vector<Long> v = spi_s.get_elems_const(idx);
      const Long idx_selected_points_send = v[0];
      const Long idx_within_field_send = v[1];
      const Long idx_buffer_send = v[2];
      qassert(0 <= idx_selected_points_send and
              idx_selected_points_send < num_selected_points_send);
      const Vector<Char> v_val = sp.get_elems_const(idx_buffer_send);
      Vector<Char> v1_val =
          view_vec[idx_selected_points_send].get_elems(idx_within_field_send);
      assign(v1_val, v_val);
    });
  }
  timer.flops += ((ssp.total_count_send + ssp.total_count_recv) / 2 +
                  ssp.total_count_local) *
                 multiplicity;
}

// ------------------------------

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp)
{
  TIMER_FLOPS("shuffle_selected_points_char(spc,spc0,ssp)");
  qassert(ssp.num_selected_points_send == 1);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(1);
  spc0_vec[0].set_view(spc0);
  shuffle_selected_points_char(spc_vec, spc0_vec, ssp);
  qassert(spc_vec.size() == 1);
  qswap_cast(spc, spc_vec[0]);
  timer.flops += ((ssp.total_count_send + ssp.total_count_recv) / 2 +
                  ssp.total_count_local) *
                 spc0.multiplicity;
}

void shuffle_selected_points_back_char(SelectedPoints<Char>& spc,
                                       const SelectedPoints<Char>& spc0,
                                       const SelectedShufflePlan& ssp)
{
  TIMER_FLOPS("shuffle_selected_points_back_char(spc,spc0,ssp)");
  qassert(ssp.num_selected_points_send == 1);
  qassert(ssp.num_selected_points_recv == 1);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(1);
  spc0_vec[0].set_view(spc0);
  shuffle_selected_points_back_char(spc_vec, spc0_vec, ssp);
  qassert(spc_vec.size() == 1);
  qswap_cast(spc, spc_vec[0]);
  timer.flops += ((ssp.total_count_send + ssp.total_count_recv) / 2 +
                  ssp.total_count_local) *
                 spc0.multiplicity;
}

void shuffle_points_selection(std::vector<PointsSelection>& psel_vec,
                              const std::vector<PointsSelection>& psel0_vec,
                              const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection(psel_vec,psel0_vec,ssp)");
  qassert(ssp.num_selected_points_send == (Long)psel0_vec.size());
  qassert(f_glb_sum((Long)psel0_vec.size()) > 0);
  const Coordinate total_site =
      f_bcast_any(psel0_vec.size() > 0 ? psel0_vec[0].total_site : Coordinate(),
                  psel0_vec.size() > 0);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(psel0_vec.size());
  for (Int i = 0; i < (Int)psel0_vec.size(); ++i) {
    qassert(psel0_vec[i].size() == ssp.n_points_selected_points_send[i]);
    qassert(psel0_vec[i].points_dist_type == ssp.points_dist_type_send);
    qassert(psel0_vec[i].total_site == total_site);
    spc0_vec[i].set_view_cast(psel0_vec[i].view_sp());
  }
  shuffle_selected_points_char(spc_vec, spc0_vec, ssp);
  psel_vec.resize(spc_vec.size());
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    Coordinate total_site2 = total_site;
    qswap_cast(psel_vec[i], spc_vec[i], total_site2);
  }
}

void shuffle_points_selection_back(
    std::vector<PointsSelection>& psel_vec,
    const std::vector<PointsSelection>& psel0_vec,
    const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection_back(psel_vec,psel0_vec,ssp)");
  qassert(ssp.num_selected_points_recv == (Long)psel0_vec.size());
  qassert(f_glb_sum((Long)psel0_vec.size()) > 0);
  const Coordinate total_site =
      f_bcast_any(psel0_vec.size() > 0 ? psel0_vec[0].total_site : Coordinate(),
                  psel0_vec.size() > 0);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(psel0_vec.size());
  for (Int i = 0; i < (Int)psel0_vec.size(); ++i) {
    qassert(psel0_vec[i].size() == ssp.n_points_selected_points_recv[i]);
    qassert(psel0_vec[i].points_dist_type == ssp.points_dist_type_recv);
    qassert(psel0_vec[i].total_site == total_site);
    spc0_vec[i].set_view_cast(psel0_vec[i].view_sp());
  }
  shuffle_selected_points_back_char(spc_vec, spc0_vec, ssp);
  psel_vec.resize(spc_vec.size());
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    Coordinate total_site2 = total_site;
    qswap_cast(psel_vec[i], spc_vec[i], total_site2);
  }
}

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection(psel,psel0,ssp)");
  qassert(psel0.points_dist_type == ssp.points_dist_type_send);
  qassert(ssp.num_selected_points_send == 1);
  qassert(ssp.num_selected_points_recv == 1);
  const Long n_points = ssp.n_points_selected_points_recv[0];
  psel.init(psel0.total_site, n_points);
  psel.points_dist_type = ssp.points_dist_type_recv;
  SelectedPoints<Char> pselc;
  const SelectedPoints<Char> pselc0(psel0.view_sp().view_as_char());
  shuffle_selected_points_char(pselc, pselc0, ssp);
  Coordinate total_site2 = psel0.total_site;
  qswap_cast(psel, pselc, total_site2);
}

void shuffle_points_selection_back(PointsSelection& psel,
                                   const PointsSelection& psel0,
                                   const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_points_selection_back(psel,psel0,ssp)");
  qassert(psel0.points_dist_type == ssp.points_dist_type_recv);
  qassert(ssp.num_selected_points_send == 1);
  qassert(ssp.num_selected_points_recv == 1);
  SelectedPoints<Char> pselc;
  const SelectedPoints<Char> pselc0(psel0.view_sp().view_as_char());
  shuffle_selected_points_back_char(pselc, pselc0, ssp);
  Coordinate total_site2 = psel0.total_site;
  qswap_cast(psel, pselc, total_site2);
}

// ------------------------------

static void set_selected_shuffle_plan_no_reorder(
    SelectedShufflePlan& ssp, const SelectedPoints<Long>& sp_instruction,
    const vector<Long>& n_points_selected_points_send,
    const PointsDistType points_dist_type_send,
    const PointsDistType points_dist_type_recv)
// Collective operation.
// partially set SelectedShufflePlan
//
// v = sp_instruction.get_elems(idx)
// v[0] = idx_selected_points_send
// v[1] = idx_within_field_send
// v[2] = id_node_send_to
// v[3] = idx_selected_points_recv
// v[4] = rank_within_field_recv
{
  TIMER("set_selected_shuffle_plan_no_reorder(ssp,sp_inst,vec,pdts,pdtr)");
  ssp.init();
  const Int num_node = get_num_node();
  const Int id_node_local = get_id_node();
  qassert(sp_instruction.initialized == true);
  qassert(sp_instruction.points_dist_type == PointsDistType::Local);
  qassert(sp_instruction.multiplicity == 5);
  ssp.points_dist_type_send = points_dist_type_send;
  ssp.points_dist_type_recv = points_dist_type_recv;
  ssp.num_selected_points_send = n_points_selected_points_send.size();
  ssp.num_selected_points_recv = 1;
  ssp.n_points_selected_points_send = n_points_selected_points_send;
  ssp.n_points_selected_points_recv.resize(1);
  set_zero(ssp.n_points_selected_points_recv);
  SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  SelectedPoints<Long>& spi_l = ssp.shuffle_idx_points_local;
  ssp.sdispls.resize(num_node);
  ssp.rdispls.resize(num_node);
  ssp.sendcounts.resize(num_node);
  ssp.recvcounts.resize(num_node);
  set_zero(ssp.sdispls);
  set_zero(ssp.rdispls);
  set_zero(ssp.sendcounts);
  set_zero(ssp.recvcounts);
  qfor(idx, sp_instruction.n_points, {
    const Vector<Long> v = sp_instruction.get_elems_const(idx);
    // const Int idx_selected_points_send = v[0];
    // const Long idx_within_send_field = v[1];
    const Int id_node_send_to = v[2];
    // const Int idx_selected_points_recv = v[3];
    // const Long rank_within_field_recv = v[4];
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
  vector<Long> c_idx_vec(MemType::Comm);
  c_idx_vec = ssp.sdispls;
  spi_l.init(ssp.total_count_local, 4, PointsDistType::Local);
  spi_s.init(ssp.total_count_send, 3, PointsDistType::Local);
  set_zero(spi_l);
  set_zero(spi_s);
  qfor(idx, sp_instruction.n_points, {
    const Vector<Long> v = sp_instruction.get_elems_const(idx);
    const Int idx_selected_points_send = v[0];
    const Long idx_within_send_field = v[1];
    const Int id_node_send_to = v[2];
    // const Int idx_selected_points_recv = v[3];
    // const Long rank_within_field_recv = v[4];
    if (id_node_send_to == id_node_local) {
      Vector<Long> v = spi_l.get_elems(idx_local);
      v[0] = idx_selected_points_send;
      v[1] = idx_within_send_field;
      v[2] = 0;
      v[3] = c_idx_local;
      c_idx_local += 1;
      idx_local += 1;
    } else {
      Vector<Long> v = spi_s.get_elems(idx_send);
      v[0] = idx_selected_points_send;
      v[1] = idx_within_send_field;
      v[2] = c_idx_vec[id_node_send_to];
      c_idx_vec[id_node_send_to] += 1;
      idx_send += 1;
    }
  });
  spi_s.set_mem_type(MemType::CommAcc);
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
  spi_r.init(ssp.total_count_recv, 3, PointsDistType::Local);
  qthread_for(idx, spi_r.n_points, {
    Vector<Long> v = spi_r.get_elems(idx);
    v[0] = 0;
    v[1] = ssp.total_count_local + idx;
    v[2] = idx;
  });
  ssp.n_points_selected_points_recv[0] =
      ssp.total_count_local + ssp.total_count_recv;
}

void set_selected_shuffle_plan(
    SelectedShufflePlan& ssp, const SelectedPoints<Long>& sp_instruction,
    const vector<Long>& n_points_selected_points_send,
    const PointsDistType points_dist_type_send,
    const PointsDistType points_dist_type_recv)
// Collective operation.
// make shuffle plan
//
// Internally call `set_selected_shuffle_plan_no_reorder`.
//
// v = sp_instruction.get_elems(idx)
// v[0] = idx_selected_points_send
// v[1] = idx_within_field_send
// v[2] = id_node_send_to
// v[3] = idx_selected_points_recv
// v[4] = rank_within_field_recv
//
// After this function, still need to set:
// total_site
// size_node_send
// coor_node_send
// size_node_recv
// coor_node_recv
{
  TIMER("set_selected_shuffle_plan(ssp,sp_inst,vec,pdts,pdtr)");
  set_selected_shuffle_plan_no_reorder(
      ssp, sp_instruction, n_points_selected_points_send, points_dist_type_send,
      points_dist_type_recv);
  //
  SelectedPoints<Long>& spi_s = ssp.shuffle_idx_points_send;
  SelectedPoints<Long>& spi_r = ssp.shuffle_idx_points_recv;
  SelectedPoints<Long>& spi_l = ssp.shuffle_idx_points_local;
  spi_s.set_mem_type(MemType::CommAcc);
  spi_r.set_mem_type(MemType::CommAcc);
  spi_l.set_mem_type(MemType::CommAcc);
  // Make a sampling send data
  std::vector<SelectedPoints<Long>> spc0_vec(ssp.num_selected_points_send);
  for (Int i = 0; i < ssp.num_selected_points_send; ++i) {
    spc0_vec[i].set_mem_type(MemType::Comm);
    spc0_vec[i].init(ssp.n_points_selected_points_send[i], 2,
                     ssp.points_dist_type_send);
    set_zero(spc0_vec[i]);
  }
  qthread_for(idx, sp_instruction.n_points, {
    const Vector<Long> v = sp_instruction.get_elems_const(idx);
    const Int idx_selected_points_send = v[0];
    const Long idx_within_send_field = v[1];
    // const Int id_node_send_to = v[2];
    const Int idx_selected_points_recv = v[3];
    const Long rank_within_field_recv = v[4];
    Vector<Long> v_s =
        spc0_vec[idx_selected_points_send].get_elems(idx_within_send_field);
    v_s[0] = idx_selected_points_recv;
    v_s[1] = rank_within_field_recv;
  });
  std::vector<SelectedPoints<Long>> spc_vec(1);
  spc_vec[0].set_mem_type(MemType::CommAcc);
  shuffle_selected_points(spc_vec, spc0_vec, ssp);
  qassert(spc_vec.size() == 1);
  SelectedPoints<Long>& spc = spc_vec[0];
  std::vector<Long> n_points_selected_points_recv(ssp.num_selected_points_recv,
                                                  0);
  std::vector<std::vector<std::array<Long, 2>>> idx_arr_vec_vec(
      ssp.num_selected_points_recv);
  qfor(idx, spc.n_points, {
    const Vector<Long> v = spc.get_elems_const(idx);
    const Int idx_selected_points_recv = v[0];
    const Long rank_within_field_recv = v[1];
    if (idx_selected_points_recv >= ssp.num_selected_points_recv) {
      ssp.num_selected_points_recv = idx_selected_points_recv + 1;
      n_points_selected_points_recv.resize(ssp.num_selected_points_recv);
      idx_arr_vec_vec.resize(ssp.num_selected_points_recv);
    }
    n_points_selected_points_recv[idx_selected_points_recv] += 1;
    std::array<Long, 2> arr;
    arr[0] = rank_within_field_recv;
    arr[1] = idx;
    idx_arr_vec_vec[idx_selected_points_recv].push_back(arr);
  });
  ssp.n_points_selected_points_recv = n_points_selected_points_recv;
  qassert(ssp.n_points_selected_points_recv.size() ==
          ssp.num_selected_points_recv);
  Long n_points = 0;
  qfor(i, ssp.num_selected_points_recv, {
    qassert((Long)idx_arr_vec_vec[i].size() ==
            ssp.n_points_selected_points_recv[i]);
    n_points += ssp.n_points_selected_points_recv[i];
  });
  qassert(n_points == spc.n_points);
  qthread_for(i, ssp.num_selected_points_recv, {
    std::sort(idx_arr_vec_vec[i].begin(), idx_arr_vec_vec[i].end());
  });
  SelectedPoints<Long> sp_idx;
  sp_idx.set_mem_type(MemType::Comm);
  sp_idx.init(spc.n_points, 2, PointsDistType::Local);
  set_zero(sp_idx);
  qfor(i, ssp.num_selected_points_recv, {
    const Int idx_selected_points_recv = i;
    const std::vector<std::array<Long, 2>>& idx_arr_vec = idx_arr_vec_vec[i];
    qthread_for(j, (Long)idx_arr_vec.size(), {
      const Long idx_within_field_recv = j;
      const std::array<Long, 2>& arr = idx_arr_vec[j];
      // const Long rank_within_field_recv = arr[0];
      const Long idx = arr[1];
      Vector<Long> v = sp_idx.get_elems(idx);
      v[0] = idx_selected_points_recv;
      v[1] = idx_within_field_recv;
    });
  });
  sp_idx.set_mem_type(MemType::CommAcc);
  qmem_for(idx, spi_l.n_points, MemType::CommAcc, {
    Vector<Long> v = spi_l.get_elems(idx);
    qassert(v[2] == 0);
    const Long idx_init = v[3];
    const Vector<Long> v_idx = sp_idx.get_elems_const(idx_init);
    const Int idx_selected_points_recv = v_idx[0];
    const Long idx_within_field_recv = v_idx[1];
    v[2] = idx_selected_points_recv;
    v[3] = idx_within_field_recv;
  });
  qmem_for(idx, spi_r.n_points, MemType::CommAcc, {
    Vector<Long> v = spi_r.get_elems(idx);
    qassert(v[0] == 0);
    const Long idx_init = v[1];
    const Vector<Long> v_idx = sp_idx.get_elems_const(idx_init);
    const Int idx_selected_points_recv = v_idx[0];
    const Long idx_within_field_recv = v_idx[1];
    v[0] = idx_selected_points_recv;
    v[1] = idx_within_field_recv;
  });
}

// ------------------------------

void set_selected_shuffle_instruction_r_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send,
    const std::vector<PointsSelection>& psel_vec, const RngState& rs)
// Shuffle the data randomly based on `gindex` and `rs`.
//
// n_points_selected_points_send.size() == psel_vec.size()
// n_points_selected_points_send[idx_selected_points_send] ==
// psel_vec[idx_selected_points_send].size()
//
// v = sp_instruction.get_elems(idx)
// v[0] = idx_selected_points_send
// v[1] = idx_within_field_send
// v[2] = id_node_send_to
// v[3] = idx_selected_points_recv
// v[4] = rank_within_field_recv
{
  TIMER(
      "set_selected_shuffle_instruction_r_from_l(sp_inst,vec,pdt,psel_vec,rs)");
  qassert(f_glb_sum((Long)psel_vec.size()) > 0);
  const Int num_node = get_num_node();
  n_points_selected_points_send.clear();
  n_points_selected_points_send.set_mem_type(MemType::Comm);
  n_points_selected_points_send.resize(psel_vec.size());
  points_dist_type_send = static_cast<PointsDistType>(f_bcast_any(
      psel_vec.size() > 0 ? static_cast<Int>(psel_vec[0].points_dist_type) : 0,
      psel_vec.size() > 0));
  Long n_points = 0;
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    qassert(points_dist_type_send == psel_vec[i].points_dist_type);
    n_points_selected_points_send[i] = psel_vec[i].size();
    n_points += psel_vec[i].size();
  }
  sp_instruction.init();
  sp_instruction.set_mem_type(MemType::Cpu);
  sp_instruction.init(n_points, 5, PointsDistType::Local);
  const RngState rs_shuffle = rs.split("shuffle_r_from_l");
  Long n_points_processed = 0;
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    const PointsSelection& psel = psel_vec[i];
    const Int idx_selected_points_send = i;
    const Int idx_selected_points_recv = i;
    qthread_for(idx, psel.size(), {
      const Long idx_within_send_field = idx;
      const Coordinate& xg = psel[idx];
      const Long gindex = index_from_coordinate(xg, psel.total_site);
      RngState rsi = rs_shuffle.newtype(gindex);
      const Int id_node_send_to = rand_gen(rsi) % num_node;
      const Long rank_within_field_recv = gindex;
      qassert(0 <= id_node_send_to);
      qassert(id_node_send_to < num_node);
      Vector<Long> v = sp_instruction.get_elems(n_points_processed + idx);
      v[0] = idx_selected_points_send;
      v[1] = idx_within_send_field;
      v[2] = id_node_send_to;
      v[3] = idx_selected_points_recv;
      v[4] = rank_within_field_recv;
    });
    n_points_processed += psel.size();
  }
  qassert(n_points_processed == n_points);
}

void set_selected_shuffle_plan_r_from_l(
    SelectedShufflePlan& ssp, const std::vector<PointsSelection>& psel_vec,
    const std::vector<Geometry>& geo_vec, const RngState& rs)
// Collective operation.
// make shuffle plan
// ssp.points_dist_type_recv = PointsDistType::Random
// psel_vec[i].points_dist_type == PointsDistType::Local (or other types)
// Sort the shuffled points by order of the gindex of points.
{
  TIMER("set_selected_shuffle_plan_r_from_l(ssp,psel_vec,geo_vec,rs)");
  SelectedPoints<Long> sp_instruction;
  vector<Long> n_points_selected_points_send;
  PointsDistType points_dist_type_send;
  set_selected_shuffle_instruction_r_from_l(
      sp_instruction, n_points_selected_points_send, points_dist_type_send,
      psel_vec, rs);
  const PointsDistType points_dist_type_recv = PointsDistType::Random;
  set_selected_shuffle_plan(ssp, sp_instruction, n_points_selected_points_send,
                            points_dist_type_send, points_dist_type_recv);
  ssp.total_site =
      f_bcast_any(psel_vec.size() > 0 ? psel_vec[0].total_site : Coordinate(),
                  psel_vec.size() > 0);
  ssp.size_node_send.resize(geo_vec.size());
  ssp.coor_node_send.resize(geo_vec.size());
  for (Int i = 0; i < (Int)geo_vec.size(); ++i) {
    const Geometry& geo = geo_vec[i];
    qassert(ssp.total_site == geo.total_site());
    ssp.size_node_send[i] = geo.geon.size_node;
    ssp.coor_node_send[i] = geo.geon.coor_node;
  }
}

void set_selected_shuffle_plan_r_from_l(SelectedShufflePlan& ssp,
                                        const PointsSelection& psel,
                                        const Geometry& geo, const RngState& rs)
// Collective operation.
// make shuffle plan
// ssp.points_dist_type_recv = PointsDistType::Random
// psel.points_dist_type == PointsDistType::Local (or other types)
// Sort the shuffled points by order of the gindex of points.
{
  TIMER("set_selected_shuffle_plan_r_from_l(ssp,psel,geo,rs)");
  std::vector<PointsSelection> psel_vec(1);
  std::vector<Geometry> geo_vec(1);
  psel_vec[0] = psel;
  geo_vec[0] = geo;
  set_selected_shuffle_plan_r_from_l(ssp, psel_vec, geo_vec, rs);
}

// ------------------------------

Long id_node_from_t_slice_id_field(const Int t_slice, const Int t_size,
                                   const Int id_field, const Int num_field,
                                   const Int num_node)
{
  const Int n_t_slice_per_node = (t_size * num_field - 1) / num_node + 1;
  const Int index = t_slice * num_field + id_field;
  const Int id_node = index / n_t_slice_per_node;
  qassert(0 <= id_node);
  qassert(id_node < num_node);
  return id_node;
}

Long idx_sp_from_t_slice_id_field(const Int t_slice, const Int t_size,
                                  const Int id_field, const Int num_field,
                                  const Int num_node)
{
  const Int n_t_slice_per_node = (t_size * num_field - 1) / num_node + 1;
  const Int index = t_slice * num_field + id_field;
  const Int idx_sp = index % n_t_slice_per_node;
  qassert(0 <= idx_sp);
  qassert(idx_sp < n_t_slice_per_node);
  return idx_sp;
}

void set_selected_shuffle_instruction_t_slice_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send,
    const std::vector<PointsSelection>& psel_vec)
// n_points_selected_points_send.size() == psel_vec.size()
// n_points_selected_points_send[idx_selected_points_send] ==
// psel_vec[idx_selected_points_send].size()
//
// v = sp_instruction.get_elems(idx)
// v[0] = idx_selected_points_send
// v[1] = idx_within_field_send
// v[2] = id_node_send_to
// v[3] = idx_selected_points_recv
// v[4] = rank_within_field_recv
{
  TIMER(
      "set_selected_shuffle_instruction_t_slice_from_l(sp_inst,vec,pdt,psel_"
      "vec)");
  qassert(f_glb_sum((Long)psel_vec.size()) > 0);
  const Int num_field = psel_vec.size();
  const Int num_node = get_num_node();
  n_points_selected_points_send.clear();
  n_points_selected_points_send.set_mem_type(MemType::Comm);
  n_points_selected_points_send.resize(psel_vec.size());
  points_dist_type_send = static_cast<PointsDistType>(f_bcast_any(
      psel_vec.size() > 0 ? static_cast<Int>(psel_vec[0].points_dist_type) : 0,
      psel_vec.size() > 0));
  const Coordinate total_site = f_bcast_any(
      psel_vec[0].size() > 0 ? psel_vec[0].total_site : Coordinate(),
      psel_vec[0].size() > 0);
  const Int t_size = total_site[3];
  Long n_points = 0;
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    qassert(points_dist_type_send == psel_vec[i].points_dist_type);
    qassert(total_site == psel_vec[i].total_site);
    n_points_selected_points_send[i] = psel_vec[i].size();
    n_points += psel_vec[i].size();
  }
  sp_instruction.init();
  sp_instruction.set_mem_type(MemType::Cpu);
  sp_instruction.init(n_points, 5, PointsDistType::Local);
  Long n_points_processed = 0;
  for (Int i = 0; i < (Int)psel_vec.size(); ++i) {
    const Int id_field = i;
    const PointsSelection& psel = psel_vec[i];
    const Int idx_selected_points_send = i;
    qthread_for(idx, psel.size(), {
      const Long idx_within_send_field = idx;
      const Coordinate& xg = psel[idx];
      const Int t_slice = xg[3];
      const Long gindex = index_from_coordinate(xg, psel.total_site);
      const Int id_node_send_to = id_node_from_t_slice_id_field(
          t_slice, t_size, id_field, num_field, num_node);
      const Int idx_selected_points_recv = idx_sp_from_t_slice_id_field(
          t_slice, t_size, id_field, num_field, num_node);
      const Long rank_within_field_recv = gindex;
      qassert(0 <= id_node_send_to);
      qassert(id_node_send_to < num_node);
      Vector<Long> v = sp_instruction.get_elems(n_points_processed + idx);
      v[0] = idx_selected_points_send;
      v[1] = idx_within_send_field;
      v[2] = id_node_send_to;
      v[3] = idx_selected_points_recv;
      v[4] = rank_within_field_recv;
    });
    n_points_processed += psel.size();
  }
  qassert(n_points_processed == n_points);
}

void set_selected_shuffle_plan_t_slice_from_l(
    SelectedShufflePlan& ssp, const std::vector<PointsSelection>& psel_vec,
    const std::vector<Geometry>& geo_vec)
// Collective operation.
// make shuffle plan
// ssp.points_dist_type_recv = PointsDistType::Local
// psel_vec[i].points_dist_type == PointsDistType::Local
// Sort the shuffled points by order of the gindex of points.
{
  TIMER("set_selected_shuffle_plan_t_slice_from_l(ssp,psel_vec,geo_vec)");
  SelectedPoints<Long> sp_instruction;
  vector<Long> n_points_selected_points_send;
  PointsDistType points_dist_type_send;
  set_selected_shuffle_instruction_t_slice_from_l(
      sp_instruction, n_points_selected_points_send, points_dist_type_send,
      psel_vec);
  PointsDistType points_dist_type_recv = PointsDistType::Local;
  if (points_dist_type_send == PointsDistType::Full) {
    points_dist_type_recv = PointsDistType::Full;
  }
  set_selected_shuffle_plan(ssp, sp_instruction, n_points_selected_points_send,
                            points_dist_type_send, points_dist_type_recv);
  ssp.total_site =
      f_bcast_any(psel_vec.size() > 0 ? psel_vec[0].total_site : Coordinate(),
                  psel_vec.size() > 0);
  ssp.size_node_send.resize(geo_vec.size());
  ssp.coor_node_send.resize(geo_vec.size());
  for (Int i = 0; i < (Int)geo_vec.size(); ++i) {
    const Geometry& geo = geo_vec[i];
    qassert(ssp.total_site == geo.total_site());
    ssp.size_node_send[i] = geo.geon.size_node;
    ssp.coor_node_send[i] = geo.geon.coor_node;
  }
  const Int num_field = psel_vec.size();
  const Int id_node = get_id_node();
  const Int num_node = get_num_node();
  const Int t_size = ssp.total_site[3];
  const Coordinate size_node = Coordinate(1, 1, 1, t_size);
  for (Int t_slice = 0; t_slice < t_size; ++t_slice) {
    const Coordinate coor_node = Coordinate(0, 0, 0, t_slice);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      if (id_node == id_node_from_t_slice_id_field(t_slice, t_size, id_field,
                                                   num_field, num_node)) {
        const Int idx_sp = idx_sp_from_t_slice_id_field(
            t_slice, t_size, id_field, num_field, num_node);
        const Int num_sp = std::max(idx_sp + 1, (Int)ssp.size_node_recv.size());
        ssp.size_node_recv.resize(num_sp);
        ssp.coor_node_recv.resize(num_sp);
        ssp.size_node_recv[idx_sp] = size_node;
        ssp.coor_node_recv[idx_sp] = coor_node;
      }
    }
  }
}

// ------------------------------

void set_selected_shuffle_instruction_dist_t_slice_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send, const PointsSelection& psel,
    const Int num_field)
// n_points_selected_points_send.size() == 1
// n_points_selected_points_send[0] == psel.size()
//
// v = sp_instruction.get_elems(idx)
// v[0] = idx_selected_points_send = 0
// v[1] = idx_within_field_send
// v[2] = id_node_send_to
// v[3] = idx_selected_points_recv
// v[4] = rank_within_field_recv
{
  TIMER(
      "set_selected_shuffle_instruction_dist_t_slice_from_l(sp_inst,vec,pdt,"
      "psel,num_field)");
  n_points_selected_points_send.clear();
  n_points_selected_points_send.set_mem_type(MemType::Comm);
  n_points_selected_points_send.resize(1);
  n_points_selected_points_send[0] = psel.size();
  points_dist_type_send = psel.points_dist_type;
  const Coordinate total_site = psel.total_site;
  const Int t_size = total_site[3];
  const Int num_node = get_num_node();
  // id_node_vec_vec[t_slice] =>
  // vector of (id_node_send_to, idx_selected_points_recv,) that should contain
  // this t_slice
  std::vector<std::vector<std::array<Int, 2>>> id_node_vec_vec(t_size);
  // t_slice_vec_vec[id_node] =>
  // vector of t_slice that should be contained in this id_node
  std::vector<std::vector<Int>> t_slice_vec_vec(num_node);
  for (Int t_slice = 0; t_slice < t_size; ++t_slice) {
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      const Int id_node = id_node_from_t_slice_id_field(
          t_slice, t_size, id_field, num_field, num_node);
      std::vector<Int>& t_slice_vec = t_slice_vec_vec[id_node];
      bool has_t_slice = false;
      for (Int i = 0; i < (Int)t_slice_vec.size(); ++i) {
        if (t_slice_vec[i] == t_slice) {
          has_t_slice = true;
          break;
        }
      }
      if (not has_t_slice) {
        std::vector<std::array<Int, 2>>& id_node_vec = id_node_vec_vec[t_slice];
        std::array<Int, 2> id_node_arr;
        id_node_arr[0] = id_node;
        id_node_arr[1] = t_slice_vec.size();
        t_slice_vec.push_back(t_slice);
        id_node_vec.push_back(id_node_arr);
      }
    }
  }
  const Int idx_selected_points_send = 0;
  std::vector<std::array<Long, 5>> instruction_vec;
  qfor(idx, psel.size(), {
    const Long idx_within_send_field = idx;
    const Coordinate& xg = psel[idx];
    const Int t_slice = xg[3];
    const Long gindex = index_from_coordinate(xg, psel.total_site);
    const std::vector<std::array<Int, 2>>& id_node_vec =
        id_node_vec_vec[t_slice];
    for (Int i = 0; i < (Int)id_node_vec.size(); ++i) {
      const Int id_node_send_to = id_node_vec[i][0];
      const Int idx_selected_points_recv = id_node_vec[i][1];
      const Long rank_within_field_recv = gindex;
      qassert(0 <= id_node_send_to);
      qassert(id_node_send_to < num_node);
      std::array<Long, 5> v;
      v[0] = idx_selected_points_send;
      v[1] = idx_within_send_field;
      v[2] = id_node_send_to;
      v[3] = idx_selected_points_recv;
      v[4] = rank_within_field_recv;
      instruction_vec.push_back(v);
    }
  });
  sp_instruction.init();
  sp_instruction.set_mem_type(MemType::Cpu);
  sp_instruction.init(instruction_vec.size(), 5, PointsDistType::Local);
  qthread_for(idx, sp_instruction.n_points, {
    Vector<Long> v = sp_instruction.get_elems(idx);
    assign(v, instruction_vec[idx]);
  });
}

void set_selected_shuffle_plan_dist_t_slice_from_l(SelectedShufflePlan& ssp,
                                                   const PointsSelection& psel,
                                                   const Geometry& geo,
                                                   const Int num_field)
// Collective operation.
// make shuffle plan
// ssp.points_dist_type_recv = PointsDistType::Local
// psel_vec[i].points_dist_type == PointsDistType::Local
// Sort the shuffled points by order of the gindex of points.
{
  TIMER("set_selected_shuffle_plan_dist_t_slice_from_l(ssp,psel,geo,num_field)")
  SelectedPoints<Long> sp_instruction;
  vector<Long> n_points_selected_points_send;
  PointsDistType points_dist_type_send;
  set_selected_shuffle_instruction_dist_t_slice_from_l(
      sp_instruction, n_points_selected_points_send, points_dist_type_send,
      psel, num_field);
  PointsDistType points_dist_type_recv = PointsDistType::Local;
  if (points_dist_type_send == PointsDistType::Full) {
    points_dist_type_recv = PointsDistType::Full;
  }
  set_selected_shuffle_plan(ssp, sp_instruction, n_points_selected_points_send,
                            points_dist_type_send, points_dist_type_recv);
  ssp.total_site = psel.total_site;
  ssp.size_node_send.resize(1);
  ssp.coor_node_send.resize(1);
  ssp.size_node_send[0] = geo.geon.size_node;
  ssp.coor_node_send[0] = geo.geon.coor_node;
  const Int id_node = get_id_node();
  const Int num_node = get_num_node();
  const Int t_size = ssp.total_site[3];
  const Coordinate size_node = Coordinate(1, 1, 1, t_size);
  std::vector<bool> b_vec(t_size, false);
  Int num_sp = 0;
  for (Int t_slice = 0; t_slice < t_size; ++t_slice) {
    const Coordinate coor_node = Coordinate(0, 0, 0, t_slice);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      if (id_node == id_node_from_t_slice_id_field(t_slice, t_size, id_field,
                                                   num_field, num_node)) {
        if (b_vec[t_slice]) {
          continue;
        }
        b_vec[t_slice] = true;
        num_sp += 1;
        ssp.size_node_recv.resize(num_sp);
        ssp.coor_node_recv.resize(num_sp);
        ssp.size_node_recv[num_sp - 1] = size_node;
        ssp.coor_node_recv[num_sp - 1] = coor_node;
      }
    }
  }
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
