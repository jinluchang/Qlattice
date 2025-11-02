// utils_field_expand.cpp
// Gen Wang
// Apr. 2024

#include <qlat/vector_utils/utils_field_expand.h>

namespace qlat
{

void setup_expand(const Geometry& geo, const Int multiplicity, qlat::vector<Long >& pack_send, qlat::vector<Long >& pack_recv, const SetMarksField& set_marks_field, const std::string& tag){
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long Nsend = plan.total_send_size;
  const Long Nrecv = plan.total_recv_size;
  //printf("setup %8d %8d \n", int(Nsend * 2), int(Nrecv * 2));
  pack_send.resize( Nsend * 2 );
  pack_recv.resize( Nrecv * 2 );
  //
  {
    TIMER("refresh setup index");
    Long cur = 0;
    for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i){
      const CommPackInfo& cpi = plan.send_pack_infos[i];
      #pragma omp parallel for
      for(Int off=0;off<cpi.size;off++){
        pack_send[(cur+off)*2 + 0] = cpi.buffer_idx + off;
        pack_send[(cur+off)*2 + 1] = cpi.offset + off;
      }
      cur += cpi.size;
    }
    //
         cur = 0;
    for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i){
      const CommPackInfo& cpi = plan.recv_pack_infos[i];
      #pragma omp parallel for
      for(Int off=0;off<cpi.size;off++){
        pack_recv[(cur + off)*2 + 0] = cpi.offset + off;
        pack_recv[(cur + off)*2 + 1] = cpi.buffer_idx + off;
      }
      cur += cpi.size;
    }
  }
}

bool operator<(const expand_index_Key& x, const expand_index_Key& y){
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
  if(tag == std::string("dirL")){set_tag = 500;}
  if(tag == std::string("dirR")){set_tag = 501;}
  //
  if(tag == std::string("")){set_tag = -100;}
  if(tag == std::string("dirx")){ set_tag = 0;}
  if(tag == std::string("diry")){ set_tag = 1;}
  if(tag == std::string("dirz")){ set_tag = 2;}
  if(tag == std::string("dirt")){ set_tag = 3;}
  if(tag == std::string("dirmx")){set_tag = -0 - 1;}
  if(tag == std::string("dirmy")){set_tag = -1 - 1;}
  if(tag == std::string("dirmz")){set_tag = -2 - 1;}
  if(tag == std::string("dirmt")){set_tag = -3 - 1;}
  //
  Qassert(set_tag != -10000);
  marks.init();
  marks.init(geo, multiplicity);
  set_zero(marks);
  Geometry geo_full = geo;
  geo_full.eo = 0;
  //
  #pragma omp parallel for
  for (Long index = 0; index < geo_full.local_volume(); ++index) {
    const Coordinate xl = geo_full.coordinate_from_index(index);
    for (Int dir = -4; dir < 4; ++dir)
    {
      if((set_tag >= -3 - 1 and set_tag < 4) and dir != set_tag){continue ;}
      if(set_tag == 500 and dir >= 0){continue ;}//only do left
      if(set_tag == 501 and dir <  0){continue ;}//only do right
      //
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      //
      if(set_tag >= -3 - 1){
        Qassert(geo.is_on_node(xl1));
      }// always need to be found on geo
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

}

