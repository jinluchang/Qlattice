// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/setup.h>
#include <qlat/field-utils.h>
#include <qlat/mpi.h>


#include <map>
#include <set>
#include <vector>

namespace qlat
{  //

struct API CommMarks : Field<int8_t> {
};

typedef void (*SetMarksField)(CommMarks& marks, const Geometry& geo,
                              const std::string& tag);

inline void set_marks_field_all(CommMarks& marks, const Geometry& geo,
                                const std::string& tag)
// tag is not used
{
  TIMER_VERBOSE("set_marks_field_all");
  (void)tag;
  marks.init();
  marks.init(geo);
  set_zero(marks);
#pragma omp parallel for
  for (long offset = 0; offset < geo.local_volume_expanded() * geo.multiplicity;
       ++offset) {
    const Coordinate xl = geo.coordinate_from_offset(offset);
    if (not geo.is_local(xl)) {
      marks.get_elem_offset(offset) = 1;
    }
  }
}

inline void set_marks_field_1(CommMarks& marks, const Geometry& geo,
                              const std::string& tag)
// tag is not used
{
  TIMER_VERBOSE("set_marks_field_1");
  (void)tag;
  marks.init();
  marks.init(geo);
  set_zero(marks);
  Geometry geo_full = geo;
  geo_full.eo = 0;
#pragma omp parallel for
  for (long index = 0; index < geo_full.local_volume(); ++index) {
    const Coordinate xl = geo_full.coordinate_from_index(index);
    for (int dir = -4; dir < 4; ++dir) {
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      if (geo.is_on_node(xl1) and !geo.is_local(xl1)) {
        Vector<int8_t> v = marks.get_elems(xl1);
        for (int m = 0; m < geo.multiplicity; ++m) {
          v[m] = 1;
        }
      }
    }
  }
}

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

struct API CommPackInfo {
  long offset;
  long buffer_idx;
  long size;
};

struct API CommMsgInfo {
  int id_node;
  long buffer_idx;
  long size;
};

struct API CommPlan {
  long total_send_size;  // send buffer size
  std::vector<CommMsgInfo> send_msg_infos;
  std::vector<CommPackInfo> send_pack_infos;
  long total_recv_size;  // recv buffer size
  std::vector<CommMsgInfo> recv_msg_infos;
  std::vector<CommPackInfo> recv_pack_infos;
};

struct API CommPlanKey {
  std::string key;
  SetMarksField set_marks_field;
  std::string tag;
  box<Geometry> geo;
};

const int lattice_size_multiplier = 3;
// g_offset calculated assume lattice_size_multiplier*total_site
// (This is a hack, please fix me)

inline void g_offset_id_node_from_offset(long& g_offset, int& id_node,
                                         const long offset, const Geometry& geo)
// offset is expanded
// g_offset calculated assume lattice_size_multiplier*total_site
{
  const Coordinate total_site = geo.total_site();
  const Coordinate xl = geo.coordinate_from_offset(offset);
  qassert(not geo.is_local(xl));
  qassert(geo.is_on_node(xl));
  const Coordinate xg =
      mod(geo.coordinate_g_from_l(xl), lattice_size_multiplier * total_site);
  const Coordinate coor_node = mod(xg, total_site) / geo.node_site;
  id_node = index_from_coordinate(coor_node, geo.geon.size_node);
  g_offset = index_from_coordinate(xg, lattice_size_multiplier * total_site) *
                 geo.multiplicity +
             offset % geo.multiplicity;
}

inline long offset_send_from_g_offset(const long g_offset, const Geometry& geo)
// g_offset calculated assume lattice_size_multiplier*total_site
// return offset is local
{
  const Coordinate total_site = geo.total_site();
  const Coordinate xg = coordinate_from_index(
      g_offset / geo.multiplicity, lattice_size_multiplier * total_site);
  Coordinate xl = geo.coordinate_l_from_g(xg);
  for (int mu = 0; mu < DIMN; ++mu) {
    while (xl[mu] >= geo.node_site[mu]) {
      xl[mu] -= total_site[mu];
    }
    while (xl[mu] < 0) {
      xl[mu] += total_site[mu];
    }
  }
  qassert(geo.is_local(xl));
  return geo.offset_from_coordinate(xl) + g_offset % geo.multiplicity;
}

inline long offset_recv_from_g_offset(const long g_offset, const Geometry& geo)
// g_offset calculated assume lattice_size_multiplier*total_site
// return offset is expanded
{
  const Coordinate total_site = geo.total_site();
  const Coordinate xg = coordinate_from_index(
      g_offset / geo.multiplicity, lattice_size_multiplier * total_site);
  Coordinate xl = geo.coordinate_l_from_g(xg);
  for (int mu = 0; mu < DIMN; ++mu) {
    while (xl[mu] >= geo.node_site[mu] + geo.expansion_right[mu]) {
      xl[mu] -= lattice_size_multiplier * total_site[mu];
    }
    while (xl[mu] < -geo.expansion_left[mu]) {
      xl[mu] += lattice_size_multiplier * total_site[mu];
    }
  }
  qassert(not geo.is_local(xl));
  qassert(geo.is_on_node(xl));
  return geo.offset_from_coordinate(xl) + g_offset % geo.multiplicity;
}

inline CommPlan make_comm_plan(const CommMarks& marks)
{
  TIMER_VERBOSE("make_comm_plan");
  const Geometry& geo = marks.geo();
  CommPlan ret;
  ret.total_send_size = 0;
  ret.total_recv_size = 0;
  //
  std::map<int, std::vector<long> >
      src_id_node_g_offsets;  // src node id ; vector of g_offset
  for (long offset = 0; offset < geo.local_volume_expanded() * geo.multiplicity;
       ++offset) {
    const int8_t r = marks.get_elem_offset(offset);
    if (r != 0) {
      int id_node;
      long g_offset;
      g_offset_id_node_from_offset(g_offset, id_node, offset, geo);
      if (id_node != get_id_node()) {
        qassert(0 <= id_node and id_node < get_num_node());
        src_id_node_g_offsets[id_node].push_back(g_offset);
      }
    }
  }
  //
  vector<long> src_id_node_count(
      get_num_node(), 0);  // number of total send pkgs for each node
  {
    long count = 0;
    for (std::map<int, std::vector<long> >::const_iterator it =
             src_id_node_g_offsets.begin();
         it != src_id_node_g_offsets.end(); ++it) {
      src_id_node_count[it->first] += 1;
      CommMsgInfo cmi;
      cmi.id_node = it->first;
      cmi.buffer_idx = count;
      cmi.size = it->second.size();
      ret.recv_msg_infos.push_back(cmi);
      count += cmi.size;
    }
    ret.total_recv_size = count;
    // ret.total_recv_size finish
    // ret.recv_msg_infos finish
  }
  glb_sum(get_data(src_id_node_count));
  ret.send_msg_infos.resize(src_id_node_count[get_id_node()]);
  //
  std::map<int, std::vector<long> >
      dst_id_node_g_offsets;  // dst node id ; vector of g_offset
  {
    std::vector<MPI_Request> reqs;
    {
      const int mpi_tag = 8;
      std::vector<CommMsgInfo> send_send_msg_infos(
          src_id_node_g_offsets.size());
      for (int i = 0; i < (int)ret.send_msg_infos.size(); ++i) {
        CommMsgInfo& cmi = ret.send_msg_infos[i];
        mpi_irecv(&cmi, sizeof(CommMsgInfo), MPI_BYTE, MPI_ANY_SOURCE, mpi_tag,
                 get_comm(), reqs);
      }
      int k = 0;
      for (std::map<int, std::vector<long> >::const_iterator it =
               src_id_node_g_offsets.begin();
           it != src_id_node_g_offsets.end(); ++it) {
        CommMsgInfo& cmi = send_send_msg_infos[k];
        cmi.id_node = get_id_node();
        cmi.buffer_idx = 0;
        cmi.size = it->second.size();
        mpi_isend(&cmi, sizeof(CommMsgInfo), MPI_BYTE, it->first, mpi_tag,
                  get_comm(), reqs);
        k += 1;
      }
      mpi_waitall(reqs);
      for (int i = 0; i < (int)ret.send_msg_infos.size(); ++i) {
        CommMsgInfo& cmi = ret.send_msg_infos[i];
        dst_id_node_g_offsets[cmi.id_node].resize(cmi.size);
      }
    }
    {
      const int mpi_tag = 9;
      int k = 0;
      long count = 0;
      for (std::map<int, std::vector<long> >::iterator it =
               dst_id_node_g_offsets.begin();
           it != dst_id_node_g_offsets.end(); ++it) {
        CommMsgInfo& cmi = ret.send_msg_infos[k];
        cmi.id_node = it->first;
        cmi.buffer_idx = count;
        cmi.size = it->second.size();
        count += cmi.size;
        mpi_irecv(it->second.data(), it->second.size(), MPI_LONG, it->first,
                  mpi_tag, get_comm(), reqs);
        k += 1;
      }
      ret.total_send_size = count;
      // ret.total_send_size finish
      k = 0;
      for (std::map<int, std::vector<long> >::const_iterator it =
               src_id_node_g_offsets.begin();
           it != src_id_node_g_offsets.end(); ++it) {
        mpi_isend((void*)it->second.data(), it->second.size(), MPI_LONG,
                  it->first, mpi_tag, get_comm(), reqs);
        k += 1;
      }
      // ret.send_msg_infos finish
      mpi_waitall(reqs);
    }
  }
  {
    long current_buffer_idx = 0;
    int k = 0;
    for (std::map<int, std::vector<long> >::const_iterator it =
             src_id_node_g_offsets.begin();
         it != src_id_node_g_offsets.end(); ++it) {
      const int src_id_node = it->first;
      const std::vector<long>& g_offsets = it->second;
      qassert(src_id_node == ret.recv_msg_infos[k].id_node);
      qassert(current_buffer_idx == ret.recv_msg_infos[k].buffer_idx);
      qassert((long)g_offsets.size() == ret.recv_msg_infos[k].size);
      long current_offset = -1;
      for (long i = 0; i < (long)g_offsets.size(); ++i) {
        const long g_offset = g_offsets[i];
        const long offset =
            offset_recv_from_g_offset(g_offset, geo);  // offset is expanded
        if (offset != current_offset) {
          CommPackInfo cpi;
          cpi.offset = offset;
          cpi.buffer_idx = current_buffer_idx;
          cpi.size = 1;
          ret.recv_pack_infos.push_back(cpi);
          current_offset = offset + 1;
          current_buffer_idx += 1;
        } else {
          CommPackInfo& cpi = ret.recv_pack_infos.back();
          cpi.size += 1;
          current_offset = offset + 1;
          current_buffer_idx += 1;
        }
      }
      k += 1;
    }
  }
  {
    long current_buffer_idx = 0;
    int k = 0;
    for (std::map<int, std::vector<long> >::const_iterator it =
             dst_id_node_g_offsets.begin();
         it != dst_id_node_g_offsets.end(); ++it) {
      const int dst_id_node = it->first;
      const std::vector<long>& g_offsets = it->second;
      qassert(dst_id_node == ret.send_msg_infos[k].id_node);
      qassert(current_buffer_idx == ret.send_msg_infos[k].buffer_idx);
      qassert((long)g_offsets.size() == ret.send_msg_infos[k].size);
      long current_offset = -1;
      for (long i = 0; i < (long)g_offsets.size(); ++i) {
        const long g_offset = g_offsets[i];
        const long offset =
            offset_send_from_g_offset(g_offset, geo);  // offset is local
        if (offset != current_offset) {
          CommPackInfo cpi;
          cpi.offset = offset;
          cpi.buffer_idx = current_buffer_idx;
          cpi.size = 1;
          ret.send_pack_infos.push_back(cpi);
          current_offset = offset + 1;
          current_buffer_idx += 1;
        } else {
          CommPackInfo& cpi = ret.send_pack_infos.back();
          cpi.size += 1;
          current_offset = offset + 1;
          current_buffer_idx += 1;
        }
      }
      k += 1;
    }
  }
  return ret;
}

inline CommPlan make_comm_plan(const CommPlanKey& cpk)
{
  CommMarks marks;
  cpk.set_marks_field(marks, cpk.geo(), cpk.tag);
  return make_comm_plan(marks);
}

API inline Cache<std::string, CommPlan>& get_comm_plan_cache()
{
  static Cache<std::string, CommPlan> cache("CommPlanCache", 32);
  return cache;
}

inline const CommPlan& get_comm_plan(const CommPlanKey& cpk)
{
  if (!get_comm_plan_cache().has(cpk.key)) {
    get_comm_plan_cache()[cpk.key] = make_comm_plan(cpk);
  }
  return get_comm_plan_cache()[cpk.key];
}

inline const CommPlan& get_comm_plan(const SetMarksField& set_marks_field,
                                     const std::string& tag,
                                     const Geometry& geo)
{
  CommPlanKey cpk;
  std::ostringstream out;
  out << (void*)set_marks_field << "," << tag << "," << geo.multiplicity << ","
      << geo.eo << "," << show(geo.expansion_left) << ","
      << show(geo.expansion_right) << "," << show(geo.total_site());
  cpk.key = out.str();
  cpk.set_marks_field = set_marks_field;
  cpk.tag = tag;
  clear(cpk.geo);
  cpk.geo.set(geo);
  return get_comm_plan(cpk);
}

template <class M>
void refresh_expanded(Field<M>& f, const CommPlan& plan)
{
  const long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded");
  timer.flops += total_bytes / 2;
  vector<M> send_buffer(plan.total_send_size);
  vector<M> recv_buffer(plan.total_recv_size);
#pragma omp parallel for
  for (long i = 0; i < (long)plan.send_pack_infos.size(); ++i) {
    const CommPackInfo& cpi = plan.send_pack_infos[i];
    memcpy(&send_buffer[cpi.buffer_idx], &f.get_elem_offset(cpi.offset),
           cpi.size * sizeof(M));
  }
  {
    sync_node();
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
    sync_node();
  }
#pragma omp parallel for
  for (long i = 0; i < (long)plan.recv_pack_infos.size(); ++i) {
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
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, f.geo());
  refresh_expanded(f, plan);
}

template <class M>
void refresh_expanded_1(Field<M>& f)
{
  const CommPlan& plan = get_comm_plan(set_marks_field_1, "", f.geo());
  refresh_expanded(f, plan);
}

inline void set_marks_field_gf_hamilton(CommMarks& marks, const Geometry& geo,
                                        const std::string& tag)
{
  TIMER_VERBOSE("set_marks_field_gf_hamilton");
  qassert(geo.multiplicity == 4);
  marks.init();
  marks.init(geo);
  set_zero(marks);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        set_marks_field_path(marks, xl,
                             make_array<int>(mu, nu, -mu - 1, -nu - 1));
        if (tag == "plaq+rect") {
          set_marks_field_path(
              marks, xl,
              make_array<int>(mu, mu, nu, -mu - 1, -mu - 1, -nu - 1));
          set_marks_field_path(
              marks, xl,
              make_array<int>(nu, nu, mu, -nu - 1, -nu - 1, -mu - 1));
        }
      }
    }
  }
}

inline void set_marks_field_gm_force(CommMarks& marks, const Geometry& geo,
                                     const std::string& tag)
{
  TIMER_VERBOSE("set_marks_field_gm_force");
  qassert(geo.multiplicity == 4);
  marks.init();
  marks.init(geo);
  set_zero(marks);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = -4; nu < 4; ++nu) {
        if (nu == mu or -nu - 1 == mu) {
          continue;
        }
        set_marks_field_path(marks, xl, make_array<int>(nu, mu, -nu - 1));
        if (tag == "plaq+rect") {
          set_marks_field_path(marks, xl,
                               make_array<int>(nu, nu, mu, -nu - 1, -nu - 1));
          set_marks_field_path(marks, xl,
                               make_array<int>(nu, mu, mu, -nu - 1, -mu - 1));
          set_marks_field_path(marks, xl,
                               make_array<int>(-mu - 1, nu, mu, mu, -nu - 1));
        }
      }
    }
  }
}

// inline void set_marks_field_m2(CommMarks& marks, const Geometry& geo,
//                                const std::string& tag)
// // tag is not used
// {
//   TIMER_VERBOSE("set_marks_field_m2");
//   marks.init();
//   marks.init(geo);
// #pragma omp parallel for
//   for (long record = 0; record < geo.local_volume_expanded(); ++record) {
//     const Coordinate xl = geo.coordinateFromRecord(record);
//     if (xl[0] < 1 - geo.expansion_left[0] or
//         xl[0] >= geo.node_site[0] + geo.expansion_right[0] -
//                      1  // 1 for checkerboarding. This function is now
//                      useless
//                         // anyway.
//         or xl[1] < 2 - geo.expansion_left[1] or
//         xl[1] >= geo.node_site[1] + geo.expansion_right[1] - 2 or
//         xl[2] < 2 - geo.expansion_left[2] or
//         xl[2] >= geo.node_site[2] + geo.expansion_right[2] - 2 or
//         xl[3] < 2 - geo.expansion_left[3] or
//         xl[3] >= geo.node_site[3] + geo.expansion_right[3] - 2) {
//       Vector<int8_t> v = marks.get_elems(xl);
//       for (int m = 0; m < geo.multiplicity; ++m) {
//         v[m] = 1;
//       }
//     }
//   }
// }

// template <class M>
// void refresh_expanded_m2(Field<M>& f)
// {
//   const CommPlan& plan = get_comm_plan(set_marks_field_m2, "", f.geo());
//   refresh_expanded(f, plan);
// }

// template <class M>
// void refresh_expanded_(Field<M>& field_comm)
// {
//   TIMER("refresh_expanded");
//   // tested for expansion = 2 case.
//   std::map<Coordinate, std::vector<M> > send_map;
//   std::map<Coordinate, int> send_map_consume;
//
//   Coordinate pos; // coordinate position of a site relative to this node
//   Coordinate local_pos; // coordinate position of a site relative to its home
//   node Coordinate node_pos; // home node coordinate of a site in node space
//
//   // populate send_map with the data that we need to send to other nodes
//   long record_size = field_comm.geo().local_volume_expanded();
//   for(long record = 0; record < record_size; record++){
//     pos = field_comm.geo().coordinateFromRecord(record);
//     if(field_comm.geo().is_local(pos)) continue;
//     for(int mu = 0; mu < DIMN; mu++){
//       local_pos[mu] = pos[mu] % field_comm.geo().node_site[mu];
//       node_pos[mu] = pos[mu] / field_comm.geo().node_site[mu];
//       if(local_pos[mu] < 0){
//         local_pos[mu] += field_comm.geo().node_site[mu];
//         node_pos[mu]--;
//       }
//     }
//     std::vector<M> &vec = send_map[node_pos];
//     for(int mu = 0; mu < field_comm.geo().multiplicity; mu++)
//       vec.push_back(field_comm.get_elems_const(local_pos)[mu]);
//   }
//
//   std::vector<M> recv_vec;
//   // will store data received from other nodes
//   // Iterate over all the nodes to which we need to send data.
//   // We ultimately copy the received data into the corresponding
//   // value of sendmap.
//   typename std::map<Coordinate, std::vector<M> >::iterator it;
//
//   // pure communication
//
//   for(it = send_map.begin(); it != send_map.end(); it++){
//     node_pos = it->first;
//     std::vector<M> &send_vec = it->second;
//     long size = send_vec.size();
//     size_t size_bytes = size * sizeof(M);
//     recv_vec.resize(std::max((long)2500, size));
//
//     M *send = send_vec.data();
//     M *recv = recv_vec.data();
//
//     Coordinate coor_this, coort, coorf;
//     int id_this, idt, idf;
//     // assuming periodic boundary condition. maybe need some fixing?
//     id_this = get_id_node();
//     coor_this = qlat::coordinate_from_index(id_this,
//         field_comm.geo().geon.size_node);
//     coort = coor_this - node_pos;
//     regularize_coordinate(coort, field_comm.geo().geon.size_node);
//     coorf = coor_this + node_pos;
//     regularize_coordinate(coorf, field_comm.geo().geon.size_node);
//
//     idt = qlat::index_from_coordinate(coort,
//     field_comm.geo().geon.size_node); idf =
//     qlat::index_from_coordinate(coorf, field_comm.geo().geon.size_node);
//
//     MPI_Request req;
//     MPI_Isend((void*)send, size_bytes, MPI_BYTE, idt, 0, get_comm(), &req);
//     const int ret = MPI_Recv((void*)recv, size_bytes, MPI_BYTE,
//         idf, 0, get_comm(), MPI_STATUS_IGNORE);
//     MPI_Wait(&req, MPI_STATUS_IGNORE);
//     qassert(!ret);
//
//     memcpy(send, recv, size_bytes);
//
//     send_map_consume[node_pos] = 0;
//   }
//   // Now send_map[node_pos] is the vector of data recieved from the node
//   // pointed to by key.
//   for(long record = 0; record < record_size; record++){
//     pos = field_comm.geo().coordinateFromRecord(record);
//     if(field_comm.geo().is_local(pos)) continue;
//     for(int mu = 0; mu < DIMN; mu++){
//       local_pos[mu] = pos[mu] % field_comm.geo().node_site[mu];
//       node_pos[mu] = pos[mu] / field_comm.geo().node_site[mu];
//       if(local_pos[mu] < 0){
//         local_pos[mu] += field_comm.geo().node_site[mu];
//         node_pos[mu]--;
//       }
//     }
//     // send_map_consume[key] keeps track of our progress in consuming the
//     // received data in sendmap[key], so that we know which offset of
//     // send_map[node_pos] corresponds to which site.
//     int consume = send_map_consume[node_pos];
//     std::vector<M> &vec = send_map[node_pos];
//     for(int mu = 0; mu < field_comm.geo().multiplicity; mu++){
//       field_comm.get_elems(pos)[mu] = vec[consume];
//       consume++;
//     }
//     send_map_consume[node_pos] = consume;
//   }
// }

}  // namespace qlat
