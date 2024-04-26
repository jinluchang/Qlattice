#include <qlat/selected-shuffle.h>

namespace qlat
{  //

void SelectedShufflePlan::init()
{
  SelectedPoints<Long>& sfi = local_shuffle_idx_field;
  sfi.init();
  total_send_count = 0;
  total_recv_count = 0;
  sdispls.clear();
  rdispls.clear();
  sendcounts.clear();
  recvcounts.clear();
}

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sf_id_node_send_to, const Long n_points,
    const RngState& rs)
{
  TIMER("set_selected_shuffle_id_node_send_to(sf_id_node_send_to,n_points,rs)");
  const Int num_node = get_num_node();
  sf_id_node_send_to.init(n_points, 1, true);
  RngState rsl = rs.split(get_id_node());
  qthread_for(idx, n_points, {
    RngState rsi = rsl.newtype(idx);
    const Int id_node_send_to = rand_gen(rsi) % num_node;
    qassert(0 <= id_node_send_to);
    qassert(id_node_send_to < num_node);
    sf_id_node_send_to.get_elem(idx) = id_node_send_to;
  });
}

void set_selected_shuffle_id_node_send_to(
    SelectedPoints<Int>& sf_id_node_send_to,
    const SelectedPoints<Coordinate>& sfx, const Coordinate& total_site,
    const RngState& rs)
{
  TIMER("set_selected_shuffle_id_node_send_to(sf_id_node_send_to,sfx,total_site,rs)");
  const Long n_points = sfx.n_points;
  qassert(sfx.multiplicity == 1);
  const Int num_node = get_num_node();
  sf_id_node_send_to.init(n_points, 1, true);
  RngState rsl = rs.split(get_id_node());
  qthread_for(idx, n_points, {
    const Coordinate& xg = sfx.get_elem(idx);
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rsi = rsl.newtype(gindex);
    const Int id_node_send_to = rand_gen(rsi) % num_node;
    qassert(0 <= id_node_send_to);
    qassert(id_node_send_to < num_node);
    sf_id_node_send_to.get_elem(idx) = id_node_send_to;
  });
}

void set_selected_shuffle_plan(SelectedShufflePlan& ssp,
                               const SelectedPoints<Int>& sf_id_node_send_to)
// Collective operation.
{
  TIMER("set_selected_shuffle_plan(ssp,sf_id_node_send_to)");
  ssp.init();
  const Int num_node = get_num_node();
  const Long n_points = sf_id_node_send_to.n_points;
  qassert(sf_id_node_send_to.initialized == true);
  qassert(sf_id_node_send_to.distributed == true);
  qassert(sf_id_node_send_to.multiplicity == 1);
  SelectedPoints<Long>& sfi = ssp.local_shuffle_idx_field;
  sfi.init(n_points, 1, true);
  set_zero(sfi);
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
    const Int id_node_send_to = sf_id_node_send_to.get_elem(idx);
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
    const Int id_node_send_to = sf_id_node_send_to.get_elem(idx);
    sfi.get_elem(idx) = c_idx_vec[id_node_send_to];
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
  SelectedPoints<Int> sf_id_node_send_to;
  set_selected_shuffle_id_node_send_to(sf_id_node_send_to, n_points, rs);
  set_selected_shuffle_plan(ssp, sf_id_node_send_to);
}

void shuffle_selected_field_char(SelectedPoints<char>& spc,
                                 const SelectedField<char>& sfc,
                                 const SelectedShufflePlan& ssp)
// const Long n_points = ssp.total_recv_count;
// const Int multiplicity = sf.geo().multiplicity;
// SelectedPoints<M> sp;
// sp.init(n_points, multiplicity);
// sp.distributed = true;
// SelectedPoints<char> spc(sp.template view_as<char>());
// SelectedField<char> sfc(sf.template view_as<char>());
{
  TIMER("shuffle_selected_field_char(spc,sfc,ssp)");
  const int multiplicity = sfc.geo().multiplicity;
  qassert(spc.initialized == true);
  qassert(spc.distributed == true);
  qassert(spc.n_points == ssp.total_recv_count);
  qassert(spc.multiplicity == multiplicity);
  const SelectedPoints<Long>& sfi = ssp.local_shuffle_idx_field;
  qassert(sfi.n_points == sfc.n_elems);
  qassert(sfi.multiplicity == 1);
  qassert(sfc.geo().is_only_local);
  SelectedField<char> sf1;
  sf1.init(sfc.geo(), sfc.n_elems, multiplicity);
  qthread_for(idx, sfc.n_elems, {
    const Vector<char> v = sfc.get_elems_const(idx);
    const Long idx1 = sfi.get_elem(idx);
    qassert(0 <= idx1);
    qassert(idx1 < sf1.n_elems);
    Vector<char> v1 = sf1.get_elems(idx1);
    assign(v1, v);
  });
  MPI_Datatype mpi_dtype;
  int mpi_ret = MPI_Type_contiguous(multiplicity, MPI_BYTE, &mpi_dtype);
  qassert(mpi_ret == 0);
  MPI_Type_commit(&mpi_dtype);
  MPI_Alltoallv(sfc.field.data(), ssp.sendcounts.data(), ssp.sdispls.data(),
                mpi_dtype, spc.points.data(), ssp.recvcounts.data(),
                ssp.rdispls.data(), mpi_dtype, get_comm());
  MPI_Type_free(&mpi_dtype);
}

void set_points_selection_from_selected_points(
    PointsSelection& psel, const SelectedPoints<Coordinate>& spx)
{
  TIMER("set_points_selection_from_selected_points(psel,spx)");
  qassert(spx.initialized);
  qassert(spx.multiplicity == 1);
  qassert(spx.n_points == spx.points.size());
  psel.initialized = spx.initialized;
  psel.distributed = spx.distributed;
  psel.xgs = spx.points;
}

void set_selected_points_from_points_selection(SelectedPoints<Coordinate>& sfx,
                                               const PointsSelection& psel)
{
  TIMER("set_selected_field_from_points_selection(sfx,psel)");
  sfx.init(psel, 1);
  qassert(sfx.n_points == psel.size());
  qthread_for(idx, psel.size(), {
    const Coordinate& xg = psel[idx];
    sfx.get_elem(idx) = xg;
  });
}

void set_selected_field_from_field_selection(SelectedField<Coordinate>& sfx,
                                             const FieldSelection& fsel)
{
  TIMER("set_selected_field_from_field_selection(sfx,fsel)");
  qassert(fsel.f_rank.initialized);
  const Geometry& geo = fsel.f_rank.geo();
  sfx.init(fsel, 1);
  qassert(sfx.n_elems == fsel.n_elems);
  qthread_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    sfx.get_elem(idx) = xg;
  });
}

}  // namespace qlat
