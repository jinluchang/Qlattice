// utils_field_expand.cpp
// Gen Wang
// Apr. 2024

#include <qlat/vector_utils/utils_field_expand.h>

namespace qlat
{

void setup_expand(const Geometry& geo, const Int multiplicity,
                  qlat::vector<Long>& pack_send, qlat::vector<Long>& pack_recv,
                  const SetMarksField& set_marks_field, const std::string& tag)
{
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long Nsend = plan.total_send_size;
  const Long Nrecv = plan.total_recv_size;
  // printf("setup %8d %8d \n", int(Nsend * 2), int(Nrecv * 2));
  pack_send.set_mem_type(MemType::Cpu);
  pack_recv.set_mem_type(MemType::Cpu);
  pack_send.resize(Nsend * 2);
  pack_recv.resize(Nrecv * 2);
  //
  {
    TIMER("refresh setup index");
    Long cur = 0;
    for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i) {
      const CommPackInfo& cpi = plan.send_pack_infos[i];
#pragma omp parallel for
      for (Int off = 0; off < cpi.size; off++) {
        pack_send[(cur + off) * 2 + 0] = cpi.buffer_idx + off;
        pack_send[(cur + off) * 2 + 1] = cpi.offset + off;
      }
      cur += cpi.size;
    }
    //
    cur = 0;
    for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i) {
      const CommPackInfo& cpi = plan.recv_pack_infos[i];
#pragma omp parallel for
      for (Int off = 0; off < cpi.size; off++) {
        pack_recv[(cur + off) * 2 + 0] = cpi.offset + off;
        pack_recv[(cur + off) * 2 + 1] = cpi.buffer_idx + off;
      }
      cur += cpi.size;
    }
  }
  // put the results to Acc
  pack_send.set_mem_type(MemType::Acc);
  pack_recv.set_mem_type(MemType::Acc);
}

bool operator<(const expand_index_Key& x, const expand_index_Key& y)
{
  Int sr = x.tag.compare(y.tag);
  if (sr < 0) {
    return false;
  }
  if (sr > 0) {
    return true;
  }
  sr = Compare_geo_less(x.geo, y.geo);
  if (sr < 0) {
    return false;
  }
  if (sr > 0) {
    return true;
  }
  return x.multiplicity < y.multiplicity;
}

void set_marks_field_dir(CommMarks& marks, const Geometry& geo,
                         const Int multiplicity, const std::string& tag)
// tag is partialy used
{
  TIMER_VERBOSE("set_marks_field_dir");
  Int set_tag = -10000;
  bool has_dir_mask = false;
  Int dir_skip[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  //
  const auto enable_dir = [&](const Int dir) {
    Qassert(dirmt <= dir and dir <= dirt);
    dir_skip[dir - dirmt] = 0;
  };
  const auto enable_pos_dir = [&](const Int mu) { enable_dir(mu); };
  const auto enable_neg_dir = [&](const Int mu) { enable_dir(-mu - 1); };
  const auto enable_pm_dir = [&](const Int mu) {
    enable_pos_dir(mu);
    enable_neg_dir(mu);
  };
  const auto enable_single_dir = [&](const Int dir) {
    set_tag = dir;
    enable_dir(dir);
  };
  //
  if (tag == std::string("dirL")) {
    set_tag = dirL;
    enable_dir(dirmx);
    enable_dir(dirmy);
    enable_dir(dirmz);
    enable_dir(dirmt);
  }
  if (tag == std::string("dirR")) {
    set_tag = dirR;
    enable_dir(dirx);
    enable_dir(diry);
    enable_dir(dirz);
    enable_dir(dirt);
  }
  if (tag == std::string("dirx")) {
    enable_single_dir(dirx);
  }
  if (tag == std::string("diry")) {
    enable_single_dir(diry);
  }
  if (tag == std::string("dirz")) {
    enable_single_dir(dirz);
  }
  if (tag == std::string("dirt")) {
    enable_single_dir(dirt);
  }
  if (tag == std::string("dirmx")) {
    enable_single_dir(dirmx);
  }
  if (tag == std::string("dirmy")) {
    enable_single_dir(dirmy);
  }
  if (tag == std::string("dirmz")) {
    enable_single_dir(dirmz);
  }
  if (tag == std::string("dirmt")) {
    enable_single_dir(dirmt);
  }
  if (tag == std::string("dirX")) {
    has_dir_mask = true;
    enable_pm_dir(0);
  }
  if (tag == std::string("dirY")) {
    has_dir_mask = true;
    enable_pm_dir(1);
  }
  if (tag == std::string("dirZ")) {
    has_dir_mask = true;
    enable_pm_dir(2);
  }
  if (tag == std::string("dirT")) {
    has_dir_mask = true;
    enable_pm_dir(3);
  }
  if (tag == std::string("dirXYZT")) {
    has_dir_mask = true;
    enable_pm_dir(0);
    enable_pm_dir(1);
    enable_pm_dir(2);
    enable_pm_dir(3);
  }
  if (tag == std::string("dirxy")) {
    has_dir_mask = true;
    enable_pos_dir(0);
    enable_pos_dir(1);
  }
  if (tag == std::string("diryz")) {
    has_dir_mask = true;
    enable_pos_dir(1);
    enable_pos_dir(2);
  }
  if (tag == std::string("dirxz")) {
    has_dir_mask = true;
    enable_pos_dir(0);
    enable_pos_dir(2);
  }
  if (tag == std::string("dirxyz")) {
    has_dir_mask = true;
    enable_pos_dir(0);
    enable_pos_dir(1);
    enable_pos_dir(2);
  }
  if (tag == std::string("dirmxy")) {
    has_dir_mask = true;
    enable_neg_dir(0);
    enable_neg_dir(1);
  }
  if (tag == std::string("dirmyz")) {
    has_dir_mask = true;
    enable_neg_dir(1);
    enable_neg_dir(2);
  }
  if (tag == std::string("dirmxz")) {
    has_dir_mask = true;
    enable_neg_dir(0);
    enable_neg_dir(2);
  }
  if (tag == std::string("dirmxyz")) {
    has_dir_mask = true;
    enable_neg_dir(0);
    enable_neg_dir(1);
    enable_neg_dir(2);
  }
  if (tag == std::string("dirXY")) {
    has_dir_mask = true;
    enable_pm_dir(0);
    enable_pm_dir(1);
  }
  if (tag == std::string("dirYZ")) {
    has_dir_mask = true;
    enable_pm_dir(1);
    enable_pm_dir(2);
  }
  if (tag == std::string("dirXZ")) {
    has_dir_mask = true;
    enable_pm_dir(0);
    enable_pm_dir(2);
  }
  if (tag == std::string("dirXYZ")) {
    has_dir_mask = true;
    enable_pm_dir(0);
    enable_pm_dir(1);
    enable_pm_dir(2);
  }
  //
  const bool check_on_node = set_tag >= -3 - 1 or has_dir_mask;
  Qassert(check_on_node);
  marks.init();
  marks.init(geo, multiplicity);
  set_zero(marks);
  Geometry geo_full = geo;
  geo_full.eo = 0;
//
#pragma omp parallel for
  for (Long index = 0; index < geo_full.local_volume(); ++index) {
    const Coordinate xl = geo_full.coordinate_from_index(index);
    for (Int dir = -4; dir < 4; ++dir) {
      if (dir_skip[dir - dirmt] != 0) {
        continue;
      }
      //
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      //
      if (check_on_node) {
        Qassert(geo.is_on_node(xl1));
      }  // always need to be found on geo
      //
      if (geo.is_on_node(xl1) and !geo.is_local(xl1)) {
        Vector<int8_t> v = marks.get_elems(xl1);
        for (Int m = 0; m < multiplicity; ++m) {
          v[m] = 1;
        }
      }
    }
  }
}

}  // namespace qlat
