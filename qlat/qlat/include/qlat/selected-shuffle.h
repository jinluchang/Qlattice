#pragma once

#include <qlat/field.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

struct SelectedShufflePlan {
  SelectedPoints<Long>
      local_shuffle_idx_field;  // Reorder field according to this idx field.
  Long total_send_count;
  Long total_recv_count;
  vector<int> sendcounts;
  vector<int> recvcounts;
  vector<int> sdispls;
  vector<int> rdispls;
  //
  void init();
};

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sf_id_node_send_to, const Long n_points,
    const RngState& rs);

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const SelectedPoints<Int>& sf_id_node_send_to);

void set_selected_shuffle_plan(SelectedShufflePlan& ssp, const Long n_points,
                               const RngState& rs);

void shuffle_selected_field_char(SelectedPoints<char>& spc,
                                 const SelectedField<char>& sfc,
                                 const SelectedShufflePlan& ssp);

void set_points_selection_from_selected_points(
    PointsSelection& psel, const SelectedPoints<Coordinate>& spx);

void set_selected_field_from_field_selection(SelectedField<Coordinate>& sfx,
                                             const FieldSelection& fsel);

template <class M>
void shuffle_selected_field(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                            const SelectedShufflePlan& ssp)
{
  TIMER("shuffle_selected_field(sp,sf,ssp)");
  const Long n_points = ssp.total_recv_count;
  const Int multiplicity = sf.geo().multiplicity;
  sp.init(n_points, multiplicity);
  sp.distributed = true;
  SelectedPoints<char> spc(sp.template view_as<char>());
  const SelectedField<char> sfc(sf.template view_as<char>());
  shuffle_selected_field_char(spc, sfc, ssp);
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
  SelectedField<Coordinate> sfx;
  set_selected_field_from_field_selection(sfx, fsel);
  SelectedPoints<Coordinate> spx;
  shuffle_selected_field(spx, sfx, ssp);
  set_points_selection_from_selected_points(psel, spx);
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
