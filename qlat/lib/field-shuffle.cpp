#define QLAT_INSTANTIATE_FIELD_SHUFFLE

#include <qlat/field-shuffle.h>

namespace qlat
{  //

std::vector<GeometryNode> make_dist_io_geons(const Coordinate& new_size_node)
{
  TIMER("make_dist_io_geons");
  const int new_num_node = product(new_size_node);
  const int num_node = get_num_node();
  std::vector<GeometryNode> ret;
  const int min_size_chunk = new_num_node / num_node;
  const int remain = new_num_node % num_node;
  const int id_node_in_shuffle =
      get_num_node() == new_num_node
          ? get_id_node()
          : get_id_node_in_shuffle(get_id_node(), new_num_node, num_node);
  const int size_chunk =
      id_node_in_shuffle < remain ? min_size_chunk + 1 : min_size_chunk;
  const int chunk_start =
      id_node_in_shuffle * min_size_chunk +
      (id_node_in_shuffle < remain ? id_node_in_shuffle : remain);
  const int chunk_end = std::min(new_num_node, chunk_start + size_chunk);
  for (int new_id_node = chunk_start; new_id_node < chunk_end; ++new_id_node) {
    GeometryNode geon;
    geon.initialized = true;
    geon.num_node = new_num_node;
    geon.id_node = new_id_node;
    geon.size_node = new_size_node;
    geon.coor_node = coordinate_from_index(new_id_node, new_size_node);
    ret.push_back(geon);
  }
  return ret;
}

std::vector<Geometry> make_dist_io_geos(const Coordinate& total_site,
                                        const int multiplicity,
                                        const Coordinate& new_size_node)
{
  TIMER("make_dist_io_geos");
  const std::vector<GeometryNode> geons = make_dist_io_geons(new_size_node);
  std::vector<Geometry> ret;
  const Coordinate new_node_site = total_site / new_size_node;
  for (int i = 0; i < (int)geons.size(); ++i) {
    Geometry geo_recv;
    geo_recv.init(geons[i], new_node_site, multiplicity);
    ret.push_back(geo_recv);
  }
  return ret;
}

ShufflePlan make_shuffle_plan(std::vector<FieldSelection>& fsels,
                              const FieldSelection& fsel,
                              const Coordinate& new_size_node)
{
  TIMER("make_shuffle_plan");
  if (new_size_node == fsel.f_rank.geo().geon.size_node) {
    fsels.clear();
    fsels.resize(1);
    fsels[0] = fsel;
    ShufflePlan sp;
    sp.is_no_shuffle = true;
    sp.new_size_node = new_size_node;
    sp.scp.global_comm_size = fsel.f_rank.geo().total_volume();
    return sp;
  }
  return make_shuffle_plan_generic(fsels, fsel, new_size_node, identity<long>);
}

ShufflePlan make_shuffle_plan(const ShufflePlanKey& spk)
{
  FieldSelection fsel;
  set_field_selection(fsel, spk.total_site);
  std::vector<FieldSelection> fsels;
  return make_shuffle_plan(fsels, fsel, spk.new_size_node);
}

ShufflePlan make_shuffle_plan_fft(const Coordinate& total_site, const int dir)
{
  TIMER_VERBOSE("make_shuffle_plan_fft");
  Geometry geo;
  geo.init(total_site, 1);
  // vol_perp_dir -> vpd
  const long vol_perp_dir =
      product(geo.node_site) / geo.node_site[dir];  // local volume perp to dir
  const Coordinate size_node = geo.geon.size_node;
  const Coordinate coor_node = geo.geon.coor_node;
  const long num_node_dir = size_node[dir];
  const long id_node_dir = coor_node[dir];
  long vpd_start;
  long vpd_size;  // share of the ``vol_perp_dir'' on this node
  split_work(vpd_start, vpd_size, vol_perp_dir, num_node_dir, id_node_dir);
  std::vector<Geometry> geos_recv(total_site[dir]);
  const Coordinate new_size_node(1, 1, total_site[dir], geo.geon.num_node);
  for (int i = 0; i < total_site[dir]; ++i) {
    const Coordinate new_coor_node(0, 0, i, geo.geon.id_node);
    const GeometryNode geon(index_from_coordinate(new_coor_node, new_size_node),
                            new_size_node);
    geos_recv[i].init(geon, Coordinate(vpd_size, 1, 1, 1), 1);
  }
  // to be shared
  ShufflePlan ret;
  // geo_send
  ret.geo_send = geo;
  // geos_recv
  ret.geos_recv = geos_recv;
  // total_send_size
  ret.scp.total_send_size = geo.local_volume();
  // send_id_node_size
  // send_new_id_node_size
  std::map<int, long> send_id_node_size;
  std::map<int, long> send_new_id_node_size;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    // custom
    const int id_node =
        get_id_node_fft(xl, geo.node_site, coor_node, size_node, dir);
    // custom
    const Coordinate new_coor_node(1, 1, xg[dir], id_node);
    const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
    //
    send_id_node_size[id_node] += 1;
    send_new_id_node_size[new_id_node] += 1;
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_id_node_size.begin();
         it != send_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.scp.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.scp.total_send_size);
  }
  // send_new_id_node_idx
  std::map<int, long> send_new_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_new_id_node_size.begin();
         it != send_new_id_node_size.end(); ++it) {
      const int new_id_node = it->first;
      const long node_size = it->second;
      send_new_id_node_idx[new_id_node] = count;
      count += node_size;
    }
    qassert(count == geo.local_volume());
    qassert(count == ret.scp.total_send_size);
  }
  // send_pack_infos
  {
    long last_buffer_idx = -1;
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      // custom
      const int id_node =
          get_id_node_fft(xl, geo.node_site, coor_node, size_node, dir);
      // custom
      const Coordinate new_coor_node(1, 1, xg[dir], id_node);
      const int new_id_node =
          index_from_coordinate(new_coor_node, new_size_node);
      //
      const long buffer_idx = send_new_id_node_idx[new_id_node];
      if (buffer_idx == last_buffer_idx and
          ret.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        ret.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        ret.send_pack_infos.push_back(pi);
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  ret.scp.total_recv_size = 0;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    ret.scp.total_recv_size += geo_recv.local_volume();
  }
  // recv_id_node_size
  std::map<int, long> recv_id_node_size;
  for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
    const Geometry& geo_recv = ret.geos_recv[i];
    for (long index = 0; index < geo_recv.local_volume(); ++index) {
      const Coordinate xl = geo_recv.coordinate_from_index(index);
      const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
      // custom
      Coordinate coor_node_orig = coordinate_from_index(xg[3], size_node);
      coor_node_orig[dir] = xg[2] / geo.node_site[dir];
      const long id_node = index_from_coordinate(coor_node_orig, size_node);
      //
      recv_id_node_size[id_node] += 1;
    }
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        ret.scp.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    if (count != ret.scp.total_recv_size) {
      qwarn(fname + ssprintf(": count = %ld", count));
      qwarn(fname + ssprintf(": ret.scp.total_recv_size = %ld",
                             ret.scp.total_recv_size));
    }
    qassert(count == ret.scp.total_recv_size);
  }
  // recv_id_node_idx
  std::map<int, long> recv_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      recv_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == ret.scp.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < ret.geos_recv.size(); ++i) {
      const Geometry& geo_recv = ret.geos_recv[i];
      long last_buffer_idx = -1;
      for (long index = 0; index < geo_recv.local_volume(); ++index) {
        const Coordinate xl = geo_recv.coordinate_from_index(index);
        const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
        // custom
        Coordinate coor_node_orig = coordinate_from_index(xg[3], size_node);
        coor_node_orig[dir] = xg[2] / geo.node_site[dir];
        const long id_node = index_from_coordinate(coor_node_orig, size_node);
        //
        const long buffer_idx = recv_id_node_idx[id_node];
        if (buffer_idx == last_buffer_idx and
            ret.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
          ret.recv_pack_infos.back().size += 1;
        } else {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          ret.recv_pack_infos.push_back(pi);
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  long num_send_packs = ret.send_pack_infos.size();
  long num_recv_packs = ret.recv_pack_infos.size();
  long num_send_msgs = ret.scp.send_msg_infos.size();
  long num_recv_msgs = ret.scp.recv_msg_infos.size();
  displayln_info(0,
                 fname + ssprintf(": num_send_packs = %10ld", num_send_packs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_packs = %10ld", num_recv_packs));
  displayln_info(0,
                 fname + ssprintf(": num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_msgs  = %10ld", num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(
      0, fname + ssprintf(": total num_send_packs = %10ld", num_send_packs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_packs = %10ld", num_recv_packs));
  displayln_info(
      0, fname + ssprintf(": total num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_msgs  = %10ld", num_recv_msgs));
  ret.scp.global_comm_size = ret.scp.total_send_size;
  glb_sum(ret.scp.global_comm_size);
  displayln_info(0, fname + ssprintf(": global_comm_size = %10ld",
                                     ret.scp.global_comm_size));
  return ret;
}

ShufflePlan make_shuffle_plan_shift(FieldSelection& fsel_shift,
                                    const FieldSelection& fsel,
                                    const Coordinate& shift,
                                    const bool is_reflect)
{
  TIMER_VERBOSE("make_shuffle_plan_shift");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate& new_size_node = geo.geon.size_node;
  ShuffleShiftGIndexMap func;
  func.total_site = geo.total_site();
  func.shift = shift;
  func.is_reflect = is_reflect;
  std::vector<FieldSelection> fsels;
  const ShufflePlan sp =
      make_shuffle_plan_generic(fsels, fsel, new_size_node, func);
  qassert(fsels.size() == 1);
  fsel_shift = fsels[0];
  return sp;
}

// old code

#if 0

ShufflePlan make_shuffle_plan_nofsel_nofunc(const ShufflePlanKey& spk)
// assume the order is preserved in transfer.
// obsolete
{
  TIMER_VERBOSE("make_shuffle_plan_nofsel_nofunc");
  const Coordinate& total_site = spk.total_site;
  const Coordinate& new_size_node = spk.new_size_node;
  const Coordinate new_node_site = total_site / new_size_node;
  qassert(new_size_node * new_node_site == total_site);
  const int new_num_node = product(new_size_node);
  Geometry geo;
  geo.init(total_site, 1);
  const int num_node = geo.geon.num_node;
  std::vector<Geometry> geos_recv =
      make_dist_io_geos(geo.total_site(), geo.multiplicity, new_size_node);
  // to be shared
  ShufflePlan sp;
  // geo_send
  sp.geo_send = geo;
  // geos_recv
  sp.geos_recv = geos_recv;
  // total_send_size
  sp.scp.total_send_size = geo.local_volume();
  // send_id_node_size
  // send_new_id_node_size
  std::map<int, long> send_id_node_size;
  std::map<int, long> send_new_id_node_size;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    // custom
    const Coordinate new_coor_node = xg / new_node_site;
    const int new_id_node = index_from_coordinate(new_coor_node, new_size_node);
    // custom
    const int id_node_in_shuffle = get_id_node_in_shuffle_from_new_id_node(
        new_id_node, new_num_node, num_node);
    //
    send_id_node_size[id_node_in_shuffle] += 1;
    send_new_id_node_size[new_id_node] += 1;
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_id_node_size.begin();
         it != send_id_node_size.end(); ++it) {
      const int id_node_in_shuffle = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = get_id_node_from_id_node_in_shuffle(
            id_node_in_shuffle, new_num_node, num_node);
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        sp.scp.send_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == geo.local_volume());
    qassert(count == sp.scp.total_send_size);
  }
  // send_new_id_node_idx
  std::map<int, long> send_new_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = send_new_id_node_size.begin();
         it != send_new_id_node_size.end(); ++it) {
      const int new_id_node = it->first;
      const long node_size = it->second;
      send_new_id_node_idx[new_id_node] = count;
      count += node_size;
    }
    qassert(count == geo.local_volume());
    qassert(count == sp.scp.total_send_size);
  }
  // send_pack_infos
  {
    long last_buffer_idx = -1;
    for (long index = 0; index < geo.local_volume(); ++index) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      // custom
      const Coordinate new_coor_node = xg / new_node_site;
      const int new_id_node =
          index_from_coordinate(new_coor_node, new_size_node);
      //
      const long buffer_idx = send_new_id_node_idx[new_id_node];
      if (buffer_idx == last_buffer_idx and
          sp.send_pack_infos.back().size < get_shuffle_max_pack_size()) {
        sp.send_pack_infos.back().size += 1;
      } else {
        last_buffer_idx = buffer_idx;
        ShufflePlanSendPackInfo pi;
        pi.field_idx = index;
        pi.buffer_idx = buffer_idx;
        pi.size = 1;
        sp.send_pack_infos.push_back(pi);
      }
      send_new_id_node_idx[new_id_node] += 1;
      last_buffer_idx += 1;
    }
  }
  // total_recv_size
  sp.scp.total_recv_size = 0;
  for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
    const Geometry& geos_recv = sp.geos_recv[i];
    sp.scp.total_recv_size += geos_recv.local_volume();
  }
  // recv_id_node_size
  std::map<int, long> recv_id_node_size;
  for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
    const Geometry& geo_recv = sp.geos_recv[i];
    for (long index = 0; index < geo_recv.local_volume(); ++index) {
      const Coordinate xl = geo_recv.coordinate_from_index(index);
      const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
      const Coordinate coor_node = xg / geo.node_site;
      // custom
      const long id_node = index_from_coordinate(coor_node, geo.geon.size_node);
      //
      recv_id_node_size[id_node] += 1;
    }
  }
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      long node_size_remain = node_size;
      while (node_size_remain > 0) {
        ShufflePlanMsgInfo mi;
        mi.id_node = id_node;
        mi.idx = count;
        mi.size = std::min(node_size_remain, get_shuffle_max_msg_size());
        sp.scp.recv_msg_infos.push_back(mi);
        node_size_remain -= mi.size;
        count += mi.size;
      }
    }
    qassert(count == sp.scp.total_recv_size);
  }
  // recv_id_node_idx
  std::map<int, long> recv_id_node_idx;
  {
    long count = 0;
    for (std::map<int, long>::const_iterator it = recv_id_node_size.begin();
         it != recv_id_node_size.end(); ++it) {
      const int id_node = it->first;
      const long node_size = it->second;
      recv_id_node_idx[id_node] = count;
      count += node_size;
    }
    qassert(count == sp.scp.total_recv_size);
  }
  // recv_pack_infos
  {
    for (size_t i = 0; i < sp.geos_recv.size(); ++i) {
      const Geometry& geo_recv = sp.geos_recv[i];
      long last_buffer_idx = -1;
      for (long index = 0; index < geo_recv.local_volume(); ++index) {
        const Coordinate xl = geo_recv.coordinate_from_index(index);
        const Coordinate xg = geo_recv.coordinate_g_from_l(xl);
        // custom
        const Coordinate coor_node = xg / geo.node_site;
        const int id_node =
            index_from_coordinate(coor_node, geo.geon.size_node);
        //
        const long buffer_idx = recv_id_node_idx[id_node];
        if (buffer_idx == last_buffer_idx and
            sp.recv_pack_infos.back().size < get_shuffle_max_pack_size()) {
          sp.recv_pack_infos.back().size += 1;
        } else {
          last_buffer_idx = buffer_idx;
          ShufflePlanRecvPackInfo pi;
          pi.local_geos_idx = i;
          pi.field_idx = index;
          pi.buffer_idx = buffer_idx;
          pi.size = 1;
          sp.recv_pack_infos.push_back(pi);
        }
        recv_id_node_idx[id_node] += 1;
        last_buffer_idx += 1;
      }
    }
  }
  long num_send_packs = sp.send_pack_infos.size();
  long num_recv_packs = sp.recv_pack_infos.size();
  long num_send_msgs = sp.scp.send_msg_infos.size();
  long num_recv_msgs = sp.scp.recv_msg_infos.size();
  displayln_info(0,
                 fname + ssprintf(": num_send_packs = %10ld", num_send_packs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_packs = %10ld", num_recv_packs));
  displayln_info(0,
                 fname + ssprintf(": num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(0,
                 fname + ssprintf(": num_recv_msgs  = %10ld", num_recv_msgs));
  glb_sum(num_send_packs);
  glb_sum(num_recv_packs);
  glb_sum(num_send_msgs);
  glb_sum(num_recv_msgs);
  displayln_info(
      0, fname + ssprintf(": total num_send_packs = %10ld", num_send_packs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_packs = %10ld", num_recv_packs));
  displayln_info(
      0, fname + ssprintf(": total num_send_msgs  = %10ld", num_send_msgs));
  displayln_info(
      0, fname + ssprintf(": total num_recv_msgs  = %10ld", num_recv_msgs));
  sp.scp.global_comm_size = sp.scp.total_send_size;
  glb_sum(sp.scp.global_comm_size);
  displayln_info(0, fname + ssprintf(": global_comm_size = %10ld",
                                     sp.scp.global_comm_size));
  return sp;
}

#endif

}  // namespace qlat
