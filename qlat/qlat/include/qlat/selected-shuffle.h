#pragma once

#include <qlat/field.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct SelectedShufflePlan {
  PointsDistType points_dist_type_send;
  PointsDistType points_dist_type_recv;
  Long num_selected_points_send;
  Long num_selected_points_recv;
  // n_points_selected_points_send.size() == num_selected_points_send
  vector<Long> n_points_selected_points_send;
  // n_points_selected_points_recv.size() == num_selected_points_recv
  vector<Long> n_points_selected_points_recv;
  // Prepare send buffer from selected points according to this idx field. (before sending)
  // multiplicity = 3 (idx_selected_points_send, idx_within_field_send, idx_buffer_send)
  SelectedPoints<Long> shuffle_idx_points_send;
  // Shuffle recv buffer to fill selected points according to this idx field. (after receiving)
  // multiplicity = 3 (idx_selected_points_recv, idx_within_field_recv, idx_buffer_recv)
  SelectedPoints<Long> shuffle_idx_points_recv;
  // Local field according to this idx field.after receiving
  // multiplicity = 4 (idx_selected_points_send, idx_within_field_send, idx_selected_points_recv, idx_within_field_recv,)
  SelectedPoints<Long> shuffle_idx_points_local;
  // shuffle_idx_points_send.n_points == total_count_send
  Long total_count_send;
  // shuffle_idx_points_recv.n_points == total_count_recv
  Long total_count_recv;
  // shuffle_idx_points_local.n_points == total_count_local
  Long total_count_local;
  // Used in mpi_alltoallv
  vector<Long> sendcounts;
  vector<Long> recvcounts;
  vector<Long> sdispls;
  vector<Long> rdispls;
  //
  void init();
};

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const PointsSelection& psel, const RngState& rs);

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp);

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_points(SelectedPoints<M>& sp,
                             const SelectedPoints<M>& sp0,
                             const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_points(sp,sp0,ssp)");
  const Long n_points = ssp.total_count_recv + ssp.total_count_local;
  const Int multiplicity = sp0.multiplicity;
  sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
  SelectedPoints<Char> spc(sp.view_as_char());
  const SelectedPoints<Char> spc0(sp0.view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
}

// -------------------

void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0,
                             const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, const SelectedField<M>& sf0,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,sf0,ssp)");
  const Long n_points = ssp.total_count_recv + ssp.total_count_local;
  const Int multiplicity = sf0.multiplicity;
  sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
  SelectedPoints<Char> spc(sp.view_as_char());
  const SelectedPoints<Char> spc0(sf0.view_sp().view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
}

// -------------------

}  // namespace qlat
