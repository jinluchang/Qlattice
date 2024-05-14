#pragma once

#include <qlat/field.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct SelectedShufflePlan {
  PointsDistType points_dist_type_send;
  PointsDistType points_dist_type_recv;
  SelectedPoints<Long> shuffle_idx_points_send;
  // Reorder field according to this idx field.before sending
  SelectedPoints<Long> shuffle_idx_points_recv;
  // Reorder field according to this idx field.after receiving
  Long total_send_count;
  Long total_recv_count;
  vector<Int> sendcounts;
  vector<Int> recvcounts;
  vector<Int> sdispls;
  vector<Int> rdispls;
  //
  void init();
};

// void set_selected_shuffle_id_node_send_to(
//     SelectedPoints<Int>& sp_id_node_send_to, const Long n_points,
//     const RngState& rs);

// void set_selected_shuffle_id_node_send_to(
//     SelectedPoints<Int>& sp_id_node_send_to, const PointsSelection& psel,
//     const RngState& rs);

// void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
//                                const SelectedPoints<Int>& sp_id_node_send_to);

// void set_selected_shuffle_plan(SelectedShufflePlan& ssp, const Long n_points,
//                                const RngState& rs);

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const PointsSelection& psel, const RngState& rs);

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp);

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp);

void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0,
                             const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_points(SelectedPoints<M>& sp,
                             const SelectedPoints<M>& sp0,
                             const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_points(sp,sp0,ssp)");
  const Long n_points = ssp.total_recv_count;
  const Int multiplicity = sp0.multiplicity;
  sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
  SelectedPoints<Char> spc(sp.view_as_char());
  const SelectedPoints<Char> spc0(sp0.view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
}

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, const SelectedField<M>& sf0,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,sf0,ssp)");
  const Long n_points = ssp.total_recv_count;
  const Int multiplicity = sf0.multiplicity;
  sp.init(n_points, multiplicity, ssp.points_dist_type_recv);
  SelectedPoints<Char> spc(sp.view_as_char());
  const SelectedPoints<Char> spc0(sf0.view_sp().view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
}

}  // namespace qlat
