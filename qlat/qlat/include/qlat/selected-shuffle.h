#pragma once

#include <qlat/field.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct SelectedShufflePlan {
  PointsDistType points_dist_type_send;
  PointsDistType points_dist_type_recv;
  // When the coordinate is not relevant, it is set to be Coordinate().
  Coordinate total_site;
  // size_node_send and coor_node_send are relevant when points_dist_type_send is one of `Full` or `Local`.
  vector<Coordinate> size_node_send;
  vector<Coordinate> coor_node_send;
  // size_node_recv and coor_node_recv are relevant when points_dist_type_recv is one of `Full` or `Local`.
  vector<Coordinate> size_node_recv;
  vector<Coordinate> coor_node_recv;
  // num_selected_points_send and num_selected_points_recv do not need to be the same across nodes.
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

// -------------------

void shuffle_selected_points_char(
    std::vector<SelectedPoints<Char>>& spc_vec,
    const std::vector<SelectedPoints<Char>>& spc0_vec,
    const SelectedShufflePlan& ssp);

void shuffle_selected_points_back_char(
    std::vector<SelectedPoints<Char>>& spc_vec,
    const std::vector<SelectedPoints<Char>>& spc0_vec,
    const SelectedShufflePlan& ssp);

// -------------------

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp);

void shuffle_selected_points_back_char(SelectedPoints<Char>& spc,
                                       const SelectedPoints<Char>& spc0,
                                       const SelectedShufflePlan& ssp);

void shuffle_points_selection(std::vector<PointsSelection>& psel_vec,
                              const std::vector<PointsSelection>& psel0_vec,
                              const SelectedShufflePlan& ssp);

void shuffle_points_selection_back(
    std::vector<PointsSelection>& psel_vec,
    const std::vector<PointsSelection>& psel0_vec,
    const SelectedShufflePlan& ssp);

void shuffle_points_selection(PointsSelection& psel,
                              const PointsSelection& psel0,
                              const SelectedShufflePlan& ssp);

void shuffle_points_selection_back(PointsSelection& psel,
                                   const PointsSelection& psel0,
                                   const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_points(std::vector<SelectedPoints<M>>& sp_vec,
                             const std::vector<SelectedPoints<M>>& sp0_vec,
                             const SelectedShufflePlan& ssp)
// completely reset sp_vec
{
  TIMER("shuffle_selected_points(sp_vec,sp0_vec,ssp)");
  Qassert(ssp.num_selected_points_send == (Long)sp0_vec.size());
  Qassert(f_glb_sum((Long)sp0_vec.size()) > 0);
  const Int multiplicity = f_bcast_any(
      sp0_vec.size() > 0 ? sp0_vec[0].multiplicity : 0, sp0_vec.size() > 0);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(sp0_vec.size());
  for (Int i = 0; i < (Int)sp0_vec.size(); ++i) {
    Qassert(sp0_vec[i].n_points == ssp.n_points_selected_points_send[i]);
    Qassert(sp0_vec[i].multiplicity == multiplicity);
    Qassert(sp0_vec[i].points_dist_type == ssp.points_dist_type_send);
    spc0_vec[i].set_view_cast(sp0_vec[i]);
  }
  shuffle_selected_points_char(spc_vec, spc0_vec, ssp);
  sp_vec.resize(spc_vec.size());
  for (Int i = 0; i < (Int)sp_vec.size(); ++i) {
    qswap_cast(sp_vec[i], spc_vec[i]);
  }
}

template <class M>
void shuffle_selected_points_back(std::vector<SelectedPoints<M>>& sp_vec,
                                  const std::vector<SelectedPoints<M>>& sp0_vec,
                                  const SelectedShufflePlan& ssp)
// completely reset sp_vec
{
  TIMER("shuffle_selected_points_back(sp_vec,sp0_vec,ssp)");
  Qassert(ssp.num_selected_points_recv == (Long)sp0_vec.size());
  Qassert(f_glb_sum((Long)sp0_vec.size()) > 0);
  const Int multiplicity = f_bcast_any(
      sp0_vec.size() > 0 ? sp0_vec[0].multiplicity : 0, sp0_vec.size() > 0);
  std::vector<SelectedPoints<Char>> spc_vec;
  std::vector<SelectedPoints<Char>> spc0_vec(sp0_vec.size());
  for (Int i = 0; i < (Int)sp0_vec.size(); ++i) {
    Qassert(sp0_vec[i].n_points == ssp.n_points_selected_points_recv[i]);
    Qassert(sp0_vec[i].multiplicity == multiplicity);
    Qassert(sp0_vec[i].points_dist_type == ssp.points_dist_type_recv);
    spc0_vec[i].set_view_cast(sp0_vec[i]);
  }
  shuffle_selected_points_back_char(spc_vec, spc0_vec, ssp);
  sp_vec.resize(spc_vec.size());
  for (Int i = 0; i < (Int)sp_vec.size(); ++i) {
    qswap_cast(sp_vec[i], spc_vec[i]);
  }
}

template <class M>
void shuffle_selected_points(SelectedPoints<M>& sp,
                             const SelectedPoints<M>& sp0,
                             const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_points(sp,sp0,ssp)");
  Qassert(ssp.num_selected_points_send == 1);
  Qassert(ssp.num_selected_points_recv == 1);
  SelectedPoints<Char> spc;
  SelectedPoints<Char> spc0;
  spc0.set_view_cast(sp0);
  shuffle_selected_points_char(spc, spc0, ssp);
  qswap_cast(sp, spc);
}

template <class M>
void shuffle_selected_points_back(SelectedPoints<M>& sp,
                                  const SelectedPoints<M>& sp0,
                                  const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_points_back(sp,sp0,ssp)");
  Qassert(ssp.num_selected_points_send == 1);
  Qassert(ssp.num_selected_points_recv == 1);
  SelectedPoints<Char> spc;
  SelectedPoints<Char> spc0;
  spc0.set_view_cast(sp0);
  shuffle_selected_points_back_char(spc, spc0, ssp);
  qswap_cast(sp, spc);
}

// -------------------

void set_selected_shuffle_plan(
    SelectedShufflePlan& ssp, const SelectedPoints<Long>& sp_instruction,
    const vector<Long>& n_points_selected_points_send,
    const PointsDistType points_dist_type_send,
    const PointsDistType points_dist_type_recv);

// -------------------

void set_selected_shuffle_instruction_r_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send,
    const std::vector<PointsSelection>& psel_vec, const RngState& rs);

void set_selected_shuffle_plan_r_from_l(
    SelectedShufflePlan& ssp, const std::vector<PointsSelection>& psel_vec,
    const std::vector<Geometry>& geo_vec, const RngState& rs);

void set_selected_shuffle_plan_r_from_l(SelectedShufflePlan& ssp,
                                        const PointsSelection& psel,
                                        const Geometry& geo,
                                        const RngState& rs);

// -------------------

Long id_node_from_t_slice_id_field(const Int t_slice, const Int t_size,
                                   const Int id_field, const Int num_field,
                                   const Int num_node);

Long idx_sp_from_t_slice_id_field(const Int t_slice, const Int t_size,
                                  const Int id_field, const Int num_field,
                                  const Int num_node);

void set_selected_shuffle_instruction_t_slice_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send,
    const std::vector<PointsSelection>& psel_vec);

void set_selected_shuffle_plan_t_slice_from_l(
    SelectedShufflePlan& ssp, const std::vector<PointsSelection>& psel_vec,
    const std::vector<Geometry>& geo_vec);

// -------------------

void set_selected_shuffle_instruction_dist_t_slice_from_l(
    SelectedPoints<Long>& sp_instruction,
    vector<Long>& n_points_selected_points_send,
    PointsDistType& points_dist_type_send, const PointsSelection& psel,
    const Int num_field);

void set_selected_shuffle_plan_dist_t_slice_from_l(SelectedShufflePlan& ssp,
                                                   const PointsSelection& psel,
                                                   const Geometry& geo,
                                                   const Int num_field);

// -------------------

void shuffle_field_selection(PointsSelection& psel, const FieldSelection& fsel0,
                             const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, const SelectedField<M>& sf0,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,sf0,ssp)");
  Qassert(ssp.num_selected_points_send == 1);
  Qassert(ssp.num_selected_points_recv == 1);
  SelectedPoints<Char> spc;
  const SelectedPoints<Char> spc0(sf0.view_sp().view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
  qswap_cast(sp, spc);
}

// -------------------

}  // namespace qlat
