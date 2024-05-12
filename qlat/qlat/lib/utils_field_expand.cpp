// utils_field_expand.cpp
// Gen Wang
// Apr. 2024

#include <qlat/vector_utils/utils_field_expand.h>

namespace qlat
{

void setup_expand(const Geometry& geo, const Int multiplicity, qlat::vector_acc<Long >& pack_send, qlat::vector_acc<Long >& pack_recv, const SetMarksField& set_marks_field, const std::string& tag)
{
  const CommPlan& plan = get_comm_plan(set_marks_field, tag, geo, multiplicity);
  const Long Nsend = plan.total_send_size;
  const Long Nrecv = plan.total_recv_size;
  /////printf("setup %8d %8d \n", int(Nsend * 2), int(Nrecv * 2));

  pack_send.resize( Nsend * 2 );
  pack_recv.resize( Nrecv * 2 );

  {
  TIMER("refresh setup index");
  Long cur = 0;
  for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i){
    const CommPackInfo& cpi = plan.send_pack_infos[i];
    #pragma omp parallel for
    for(int off=0;off<cpi.size;off++){
      pack_send[(cur+off)*2 + 0] = cpi.buffer_idx + off;
      pack_send[(cur+off)*2 + 1] = cpi.offset + off;
    }
    cur += cpi.size;
  }

       cur = 0;
  for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i){
    const CommPackInfo& cpi = plan.recv_pack_infos[i];
    #pragma omp parallel for
    for(int off=0;off<cpi.size;off++){
      pack_recv[(cur + off)*2 + 0] = cpi.offset + off;
      pack_recv[(cur + off)*2 + 1] = cpi.buffer_idx + off;
    }
    cur += cpi.size;
  }
  }
}

bool compare_geo(const Geometry& g0, const Geometry& g1)
{
  int equal = 1;
  if(g0.initialized           != g1.initialized ){ return 0; }
  if(g0.eo                    != g1.eo ){ return 0; }
  if(g0.is_only_local         != g1.is_only_local    ){ return 0; }

  if(g0.geon                  != g1.geon ){ return 0; }
  if(g0.node_site             != g1.node_site    ){ return 0; }
  if(g0.node_site_expanded    != g1.node_site_expanded    ){ return 0; }

  if(g0.expansion_left        != g1.expansion_left  ){ return 0; }
  if(g0.expansion_right       != g1.expansion_right ){ return 0; }

  //if(g0.total_site()    != g1.total_site()    ){ return 0; }

  return equal;
}

bool Compare_geo(const Geometry& g0, const Geometry& g1)
{
  return compare_geo(g0, g1);
}

bool compare_less(const Geometry& g0, const Geometry& g1)
{
  if(g0.total_site()    < g1.total_site()    ){  return true;}
  if(g1.total_site()    < g0.total_site()    ){  return false;}

  if(g0.expansion_left  < g1.expansion_left  ){  return true;}
  if(g1.expansion_left  < g0.expansion_left  ){  return false;}

  if(g0.expansion_right < g1.expansion_right ){  return true;}
  if(g1.expansion_right < g0.expansion_right ){  return false;}

  return false;
}

bool operator<(const expand_index_Key& x, const expand_index_Key& y)
{
  int sr = x.tag.compare(y.tag);
  if (sr < 0) {
    return false;
  }
  if (sr > 0) {
    return true;
  }
  sr = compare_less(x.geo, y.geo);
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
  int set_tag = -10000;
  if(tag == std::string("dirL")){set_tag = 500;}
  if(tag == std::string("dirR")){set_tag = 501;}

  if(tag == std::string("")){set_tag = -100;}
  if(tag == std::string("dirx")){ set_tag = 0;}
  if(tag == std::string("diry")){ set_tag = 1;}
  if(tag == std::string("dirz")){ set_tag = 2;}
  if(tag == std::string("dirt")){ set_tag = 3;}
  if(tag == std::string("dirmx")){set_tag = -0 - 1;}
  if(tag == std::string("dirmy")){set_tag = -1 - 1;}
  if(tag == std::string("dirmz")){set_tag = -2 - 1;}
  if(tag == std::string("dirmt")){set_tag = -3 - 1;}

  Qassert(set_tag != -10000);
  marks.init();
  marks.init(geo, multiplicity);
  set_zero(marks);
  Geometry geo_full = geo;
  geo_full.eo = 0;

  #pragma omp parallel for
  for (Long index = 0; index < geo_full.local_volume(); ++index) {
    const Coordinate xl = geo_full.coordinate_from_index(index);
    for (int dir = -4; dir < 4; ++dir) 
    {
      if((set_tag >= -3 - 1 and set_tag < 4) and dir != set_tag){continue ;}
      if(set_tag == 500 and dir >= 0){continue ;}//only do left
      if(set_tag == 500 and dir <  0){continue ;}//only do right

      const Coordinate xl1 = coordinate_shifts(xl, dir);

      if(set_tag >= -3 - 1){
        Qassert(geo.is_on_node(xl1));
      }// always need to be found on geo

      if (geo.is_on_node(xl1) and !geo.is_local(xl1)) {
        Vector<int8_t> v = marks.get_elems(xl1);
        for (int m = 0; m < multiplicity; ++m) {
          v[m] = 1;
        }
      }
    }
  }

}

//inline Cache<expand_index_Key, expand_index_buf >& get_expand_index_buf_cache()
//{
//  static Cache<expand_index_Key, expand_index_buf > cache("expand_index_Key", 64);
//  return cache;
//}
//
//inline expand_index_buf& get_expand_index_buf_plan(const expand_index_Key& ekey)
//{
//  if (!get_expand_index_buf_cache().has(ekey)) {
//    //Geometry geo(ekey.total_site, 1);
//    //Geometry geo_ext = geo_resize(geo, ekey.expansion_left, ekey.expansion_right);
//    get_expand_index_buf_cache()[ekey] = expand_index_buf(ekey.geo);
//  }
//  expand_index_buf& buf = get_expand_index_buf_cache()[ekey];
//  return buf;
//}

//inline expand_index_buf& get_expand_index_buf_plan(const Geometry& geo)
//{
//  expand_index_Key ekey(geo);
//  return get_expand_index_buf_plan(ekey);
//}

//template <class M>
//void refresh_expanded_GPU(Field<M>& f, int GPU = 1)
//{
//  const CommPlan& plan = get_comm_plan(set_marks_field_all, "", f.geo());
//  const Long total_bytes =
//      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
//  if (0 == total_bytes) {
//    return;
//  }
//  TIMER_FLOPS("refresh_expanded_GPU");
//  timer.flops += total_bytes / 2;
//
//  std::vector<MPI_Request> reqs_send;
//  std::vector<MPI_Request> reqs_recv;
//
//  qlat::vector_gpu<char >& sbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf0"), GPU);
//  qlat::vector_gpu<char >& rbuf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf1"), GPU);
//  //qlat::vector_gpu<char >& pack_buf = qlat::get_vector_gpu_plan<char >(0, std::string("general_buf2"), -1);
//  expand_index_buf& ebuf = get_expand_index_buf_plan(f.geo());
//  //expand_index_buf ebuf(f.geo());
//  //printf("send %8d %8d \n", int(ebuf.pack_send.size()), int( 2*plan.total_send_size));
//  //printf("recv %8d %8d \n", int(ebuf.pack_recv.size()), int( 2*plan.total_recv_size));
//  Qassert(ebuf.pack_send.size() == 2*plan.total_send_size and ebuf.pack_recv.size() == 2*plan.total_recv_size);
//  //expand_index_buf ebuf(f.geo());
//
//  const Long Nsend = plan.total_send_size;
//  const Long Nrecv = plan.total_recv_size;
//
//  sbuf.resizeL(Nsend * sizeof(M) / sizeof(char));
//  rbuf.resizeL(Nrecv * sizeof(M) / sizeof(char));
//  //pack_buf.resizeL( 2 * plan.total_send_size * sizeof(Long) / sizeof(char));
//  //pack_buf.resizeL( 2 * plan.total_recv_size * sizeof(Long) / sizeof(char));
//
//  M* sP = (M*) &sbuf[0];
//  M* rP = (M*) &rbuf[0];
//
//  Qassert(sizeof(M) % sizeof(double) == 0);
//  ////setup reciev
//  const int mpi_tag = 10;
//  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
//    const CommMsgInfo& cmi = plan.recv_msg_infos[i]; 
//    mpi_irecv(&rP[cmi.buffer_idx], cmi.size * sizeof(M)/sizeof(double), MPI_DOUBLE,
//              cmi.id_node, mpi_tag, get_comm(), reqs_recv);
//  }
//
//  //qlat::vector_acc<long > pack_infos;
//  Long* pack_send = (Long*) &ebuf.pack_send[0];
//  Long* pack_recv = (Long*) &ebuf.pack_recv[0];
//
//  //{
//  //TIMER("refresh setup index");
//  //////pack_infos.resize(2 * plan.total_send_size);
//  //Long cur = 0;
//  //for (Long i = 0; i < (Long)plan.send_pack_infos.size(); ++i){
//  //  const CommPackInfo& cpi = plan.send_pack_infos[i];
//  //  #pragma omp parallel for
//  //  for(int off=0;off<cpi.size;off++){
//  //    pack_infos[(cur+off)*2 + 0] = cpi.buffer_idx + off;
//  //    pack_infos[(cur+off)*2 + 1] = cpi.offset + off;
//  //    //send_pack_infos[i*3 + 2] = cur;
//  //  }
//  //  cur += cpi.size;
//  //}
//  //}
//
//  qGPU_for(isp, Nsend, GPU, {
//    Long ri = pack_send[isp* 2 + 0];
//    Long si = pack_send[isp* 2 + 1];
//    sP[ri] = f.get_elem_offset(si);
//  });
//
//  { 
//    //TIMER("refresh_expanded-comm-init");
//    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
//      const CommMsgInfo& cmi = plan.send_msg_infos[i];
//      mpi_isend(&sP[cmi.buffer_idx], cmi.size * sizeof(M)/sizeof(double), MPI_DOUBLE,
//                cmi.id_node, mpi_tag, get_comm(), reqs_send);
//    }
//  }
//
//  //{
//  //TIMER("refresh setup index");
//  /////const Long Nsize = 0;
//  /////pack_infos.resize(2 * plan.total_recv_size);
//  //Long cur = 0;
//  //for (Long i = 0; i < (Long)plan.recv_pack_infos.size(); ++i){
//  //  const CommPackInfo& cpi = plan.recv_pack_infos[i];
//  //  #pragma omp parallel for
//  //  for(int off=0;off<cpi.size;off++){
//  //    pack_infos[(cur + off)*2 + 0] = cpi.offset + off;
//  //    pack_infos[(cur + off)*2 + 1] = cpi.buffer_idx + off;
//  //  }
//  //  cur += cpi.size;
//  //}
//  //}
//
//  mpi_waitall(reqs_recv);////receive done and write
//  qGPU_for(isp, Nrecv, GPU, {
//    const Long ri = pack_recv[isp* 2 + 0];
//    const Long si = pack_recv[isp* 2 + 1];
//    f.get_elem_offset(ri) = rP[si];
//  });
//
//  mpi_waitall(reqs_send);
//  //safe_free_vector_gpu_plan<char >(std::string("general_buf0"), GPU);
//  //safe_free_vector_gpu_plan<char >(std::string("general_buf1"), GPU);
//}


}

