// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/field.h>

#include <map>
#include <set>
#include <vector>

namespace qlat
{  //

struct API CommMarks : Field<int8_t> {
};

typedef void (*SetMarksField)(CommMarks& marks, const Geometry& geo,
                              const Int multiplicity, const std::string& tag);

template <class Vec>
qacc void set_marks_field_path(CommMarks& marks, const Coordinate xl,
                               const Vec& path)
{
  const Geometry& geo = marks.geo();
  Coordinate xl1 = xl;
  for (int i = 0; i < (int)path.size(); ++i) {
    const int dir = path[i];
    qassert(-DIMN <= dir && dir < DIMN);
    if (0 <= dir) {
      if (not geo.is_local(xl1)) {
        marks.get_elem(xl1, dir) = 1;
      }
      xl1[dir] += 1;
    } else {
      xl1[-dir - 1] -= 1;
      if (not geo.is_local(xl1)) {
        marks.get_elem(xl1, -dir - 1) = 1;
      }
    }
  }
}

void set_marks_field_all(CommMarks& marks, const Geometry& geo, const Int multiplicity,
                         const std::string& tag);

void set_marks_field_1(CommMarks& marks, const Geometry& geo, const Int multiplicity,
                       const std::string& tag);

void set_marks_field_gf_hamilton(CommMarks& marks, const Geometry& geo, const Int multiplicity,
                                 const std::string& tag);

void set_marks_field_gm_force(CommMarks& marks, const Geometry& geo, const Int multiplicity,
                              const std::string& tag);

struct API CommPackInfo {
  Long offset;
  Long buffer_idx;
  Long size;
};

struct API CommMsgInfo {
  Int id_node;
  Long buffer_idx;
  Long size;
};

struct API CommPlan {
  Long total_send_size;  // send buffer size
  std::vector<CommMsgInfo> send_msg_infos;
  std::vector<CommPackInfo> send_pack_infos;
  Long total_recv_size;  // recv buffer size
  std::vector<CommMsgInfo> recv_msg_infos;
  std::vector<CommPackInfo> recv_pack_infos;
};

struct API CommPlanKey {
  std::string key;
  SetMarksField set_marks_field;
  std::string tag;
  box<Geometry> geo;
  Int multiplicity;
};

const int lattice_size_multiplier = 3;
// g_offset calculated assume lattice_size_multiplier*total_site
// (This is a hack, please fix me)

void g_offset_id_node_from_offset(Long& g_offset, int& id_node,
                                  const Long offset, const Geometry& geo, const Int multiplicity);

Long offset_send_from_g_offset(const Long g_offset, const Geometry& geo, const Int multiplicity);

Long offset_recv_from_g_offset(const Long g_offset, const Geometry& geo, const Int multiplicity);

CommPlan make_comm_plan(const CommMarks& marks);

CommPlan make_comm_plan(const CommPlanKey& cpk);

API inline Cache<std::string, CommPlan>& get_comm_plan_cache()
{
  static Cache<std::string, CommPlan> cache("CommPlanCache", 32);
  return cache;
}

const CommPlan& get_comm_plan(const CommPlanKey& cpk);

const CommPlan& get_comm_plan(const SetMarksField& set_marks_field,
                              const std::string& tag, const Geometry& geo,
                              const Int multiplicity);

template <class M>
void refresh_expanded(Field<M>& f, const CommPlan& plan)
{
  const Long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded");
  timer.flops += total_bytes / 2;
  vector<M> send_buffer(plan.total_send_size, MemType::Comm);
  vector<M> recv_buffer(plan.total_recv_size, MemType::Comm);
#pragma omp parallel for
  for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i) {
    const CommPackInfo& cpi = plan.send_pack_infos[i];
    memcpy(&send_buffer[cpi.buffer_idx], &f.get_elem_offset(cpi.offset),
           cpi.size * sizeof(M));
  }
  {
    SYNC_NODE();
    TIMER_FLOPS("refresh_expanded-comm");
    timer.flops +=
        (plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;
    std::vector<MPI_Request> reqs;
    {
      TIMER("refresh_expanded-comm-init");
      const int mpi_tag = 10;
      for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = plan.recv_msg_infos[i];
        mpi_irecv(&recv_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), reqs);
      }
      for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = plan.send_msg_infos[i];
        mpi_isend(&send_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), reqs);
      }
    }
    mpi_waitall(reqs);
    SYNC_NODE();
  }
#pragma omp parallel for
  for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i) {
    const CommPackInfo& cpi = plan.recv_pack_infos[i];
    memcpy(&f.get_elem_offset(cpi.offset), &recv_buffer[cpi.buffer_idx],
           cpi.size * sizeof(M));
  }
}

template <class M>
void refresh_expanded(
    Field<M>& f, const SetMarksField& set_marks_field = set_marks_field_all,
    const std::string& tag = "")
{
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan =
      get_comm_plan(set_marks_field, tag, f.geo(), f.multiplicity);
  QLAT_DIAGNOSTIC_POP;
  refresh_expanded(f, plan);
}

template <class M>
void refresh_expanded_1(Field<M>& f)
{
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan =
      get_comm_plan(set_marks_field_1, "", f.geo(), f.multiplicity);
  QLAT_DIAGNOSTIC_POP;
  refresh_expanded(f, plan);
}

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD_EXPAND
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                              \
                                                                    \
  QLAT_EXTERN template void refresh_expanded(Field<TYPENAME>& f,    \
                                             const CommPlan& plan); \
                                                                    \
  QLAT_EXTERN template void refresh_expanded(                       \
      Field<TYPENAME>& f, const SetMarksField& set_marks_field,     \
      const std::string& tag);                                      \
                                                                    \
  QLAT_EXTERN template void refresh_expanded_1(Field<TYPENAME>& f)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
