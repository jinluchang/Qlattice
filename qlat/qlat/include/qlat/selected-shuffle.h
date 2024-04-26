#pragma once

#include <qlat/field.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct SelectedShufflePlan {
  SelectedPoints<Long>
      send_shuffle_idx_points;  // Reorder field according to this idx field.
  SelectedPoints<Long>
      recv_shuffle_idx_points;  // Reorder field according to this idx field.
  Long total_send_count;
  Long total_recv_count;
  vector<Int> sendcounts;
  vector<Int> recvcounts;
  vector<Int> sdispls;
  vector<Int> rdispls;
  //
  void init();
};

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sp_id_node_send_to, const Long n_points,
    const RngState& rs);

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const SelectedPoints<Int>& sp_id_node_send_to);

void set_selected_shuffle_plan(SelectedShufflePlan& ssp, const Long n_points,
                               const RngState& rs);

void shuffle_selected_points_char(SelectedPoints<Char>& spc,
                                  const SelectedPoints<Char>& spc0,
                                  const SelectedShufflePlan& ssp);

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,sf,ssp)");
  const Long n_points = ssp.total_recv_count;
  const Int multiplicity = sf.geo().multiplicity;
  sp.init(n_points, multiplicity);
  sp.distributed = true;
  SelectedPoints<Char> spc(sp.view_as_char());
  const SelectedPoints<Char> spc0(sf.view_sp().view_as_char());
  shuffle_selected_points_char(spc, spc0, ssp);
}

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, PointsSelection& psel,
                            const SelectedField<M>& sf,
                            const FieldSelection& fsel,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,psel,sf,fsel,ssp)");
  const Long n_elems = fsel.n_elems;
  qassert(n_elems == ssp.total_send_count);
  shuffle_selected_field(sp, sf, ssp);
  psel.init(ssp.total_recv_count);
  psel.distributed = true;
  PointsSelection psel0;
  set_psel_from_fsel(psel0, fsel);
  SelectedPoints<Char> pselc(psel.view_sp().view_as_char());
  const SelectedPoints<Char> pselc0(psel0.view_sp().view_as_char());
  shuffle_selected_points_char(pselc, pselc0, ssp);
}

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, PointsSelection& psel,
                            const SelectedField<M>& sf,
                            const FieldSelection& fsel, const RngState& rs)
{
  TIMER("shuffle_selected_field(sp,psel,sf,fsel,rs)");
  const Long n_elems = fsel.n_elems;
  SelectedShufflePlan ssp;
  set_selected_shuffle_plan(ssp, n_elems, rs);
  shuffle_selected_field(sp, psel, sf, fsel, ssp);
}

}  // namespace qlat
