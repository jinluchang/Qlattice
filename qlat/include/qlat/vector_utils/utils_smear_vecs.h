// utils_smear_vecs.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_SMEAR_VECS_H
#define UTILS_SMEAR_VECS_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
#include <qlat/qcd-smear.h>
#include "utils_Vec_redistribute.h"
#include "utils_shift_vecs.h"
#include "check_fun.h"
 
namespace qlat{

///////////psrc order, bfac, c, d0, t,z,y,x
#ifdef QLAT_USE_ACC
template <class T, int bfac, int d0, int dirL>
__global__ void gauss_smear_global4(T* pres, const T* psrc, const T* gf, const T bw, const T norm,
  const unsigned long Nvol, const long* map_bufD)
{

  __shared__ T ls[(dirL*2)*3*3 + 1];
  __shared__ T ps[(dirL*2)*bfac*3*d0 + 1];
  ///tid should larger than 9
  unsigned int   tid   =  threadIdx.x;
  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int   ns = blockDim.x;
  const int dir_max = 4;

  if(index < Nvol){

  ///////(2*dirL) --> c0 , c1
  ///////TODO need to check dirL is not 3
  unsigned int off = tid;
  {
    const T* gf_t = &gf[index*(2*dir_max)*9];
    while(off < (2*dirL)*9){
      //int dir = off/9;int c = off%9;
      //ls[(c/3)*(2*dirL)*3 + dir*3 + c%3] = gf[index*(2*dirL)*9 + off];off += ns;
      //ls[off] = gf_t[off];
      int dir = off/9;
      int c0 = (off/3)%3;
      int c1 =  off%3;
      ls[c1*(2*dirL)*3 + dir*3 + c0 ] = gf_t[(dir + 1)*9 + c1*3 + c0];
      ////ls[c1*(2*dirL)*3 + dir*3 + c0 ] = gf_t[(dir + 1)*9 + c0*3 + c1];
      off += ns;
    }
  }

  const long res_off = index;
  T* wm = &pres[res_off*bfac*3*d0];

  ///////(2*dir) --> bi, c1 , d0
  for (int dir = -dirL; dir < dirL; ++dir){
    ////long src_off = map_bufD[(dir+4)*Nvol + index];
    long src_off = map_bufD[index* dirL*2 + (dir + dirL)];
    const T* src_t = &psrc[src_off*bfac*3*d0];
    //T* res_t     = &ps[(dir+dirL)*bfac*3*d0];
    ////T* res_t     = &ps[dir+dirL];
    off = tid;
    while(off < bfac*3*d0){
      //res_t[off] = src_t[off]; off+=ns;
      unsigned int bi = off/(3*d0);
      unsigned int c  = (off%(3*d0))/d0;
      unsigned int di = off%d0;
      ps[(bi*d0 + di)*(2*dirL*3) + (dir+dirL)*3 + c] = src_t[off];
      off+=ns;
    }
  }
  __syncthreads();

  T tmp = 0.0;
  off = tid;
  while(off < bfac*3*d0){
    unsigned int c0 = off/(bfac*d0);
    unsigned int bi = (off%(bfac*d0))/d0;
    unsigned int di = off%d0;

    T* a = &ls[c0*(2*dirL*3)];
    T* b = &ps[(bi*d0 + di)*(2*dirL*3)];
    //T* b = &ps[(bi*d0 + di)];

    tmp = 0.0;

    for(int dir=0;dir<(2*dirL)*3; dir++)
    {
      //tmp += b[dir*bfac*d0] * a[dir];
      tmp += b[dir] * a[dir];
    }

    wm[(bi*3 + c0)*d0 + di] = norm*(wm[(bi*3 + c0)*d0 + di ] + bw*tmp);
    off += ns;
  }

  }

}
#endif



template <typename Ty >
struct smear_fun{

  Geometry geo;
  Geometry geo_ext;
  int dirL;
  bool smear_in_time_dir;

  ////shift_vec *svec;
  Vec_redistribute *vec_rot;
  move_index mv_idx;

  /////for shifters and geo
  fft_desc_basic fd;
  fft_desc_basic fd_new;

  std::vector<qlat::vector_gpu<long > > mapvq_send;
  std::vector<qlat::vector_gpu<long > > mapvq_recv;
  qlat::vector<MPI_Request> send_reqs;
  qlat::vector<MPI_Request> recv_reqs;
  ///unsigned int bfac;
  ////unsigned int Tsize;
  unsigned long Nvol;
  unsigned long Nvol_ext;

  std::vector<qlat::vector_acc<long > > map_bufV;
  qlat::vector_acc<long > map_bufD;
  qlat::vector_acc<int> Nv,nv,mv;

  unsigned int NVmpi;
  int groupP;

  int fft_copy;

  qlat::vector_gpu<Ty > send_buffer;
  qlat::vector_gpu<Ty > recv_buffer;
  qlat::vector_gpu<Ty > gauge;
  qlat::vector_gpu<Ty > prop;
  qlat::vector_gpu<Ty > prop_buf;
  std::vector<qlat::vector_gpu<Ty > > vL;

  bool  gauge_setup_flag;
  void* gauge_check;
  /////buffers

  bool  mem_setup_flag;
  qlat::vector_acc<qlat::Complex > mom_factors;

  //template <class Tg>
  //void setup(const GaugeFieldT<Tg >& gf, const CoordinateD& mom, const bool smear_in_time_dir){
  //  qlat::vector_acc<qlat::Complex > momF(8);
  //  for (int i = 0; i < 8; ++i) {
  //    const int dir = i - 4;
  //    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
  //    momF[i] = std::polar(1.0, -phase);
  //  }
  //  gauge_setup(gf, momF);
  //}

  smear_fun(const Geometry& geo_, const bool smear_in_time_dir_){
    Nvol  = 0;
    fft_copy = 0;
    NVmpi = 0;
    groupP = 0;

    gauge_check = NULL;
    gauge_setup_flag = false;

    mem_setup_flag = false;
    dirL = 3;

    fft_desc_basic fd_(geo_);
    copy_fft_desc(fd, fd_);
    copy_fft_desc(fd_new, fd_);
    NVmpi = fd.mz*fd.my*fd.mx;

    smear_in_time_dir = smear_in_time_dir_;
    init_mem(geo_, smear_in_time_dir);

    ////Box smearing buffer size
    vL.resize(8);
    vec_rot = NULL;
  }

  void check_setup(){
    if(Nv.size() == 0){qassert(false);}
    if(Nvol == 0){qassert(false);}
    if(Nvol_ext == 0){qassert(false);}
    if(gauge.size() == 0){qassert(false);}
    if(gauge_check == NULL){qassert(false);}
    if(map_bufV.size() == 0){qassert(false);}
    if(map_bufD.size() == 0){qassert(false);}
    if(vL.size() != 8){qassert(false);}
    if(mem_setup_flag == false){qassert(false);}
    ////if(gauge.size() != Nvol*4*2*9){qassert(false);}
  };

  void get_mapvq(const std::vector<CommPackInfo> &pack_infos, std::vector<qlat::vector_gpu<long > >& mapvq, int dir=0)
  {
    std::vector<std::vector<long > > mapv(2);//mapv.resize(4*pack_infos.size());
    for (long i = 0; i < (long)pack_infos.size(); ++i) {
      const CommPackInfo& cpi = pack_infos[i];
      for(long j=0; j< cpi.size ; j++)
      {
        long bufi = cpi.buffer_idx + j;
        long fi   = cpi.offset     + j;
        if(dir == 0){mapv[0].push_back(bufi);mapv[1].push_back(  fi);}
        if(dir == 1){mapv[0].push_back(  fi);mapv[1].push_back(bufi);}
      }
    }
    mapvq.resize(2);for(unsigned int i=0;i<2;i++){mapvq[i].copy_from(mapv[i], 1);}
  }

  /////index for MPI buffers
  void init_mem(const Geometry& geo_, bool smear_in_time_dir_ = false)
  {
    if(geo == geo_ and smear_in_time_dir == smear_in_time_dir_ and mem_setup_flag == true){return ;}
    geo = geo_;
    ////move to default geo
    geo.resize(Coordinate(0, 0, 0, 0), Coordinate(0, 0, 0, 0));geo.multiplicity = 1;geo.eo=0;
    Geometry geo1 = geo;
    if(!smear_in_time_dir){geo1.resize(Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));}
    if( smear_in_time_dir){geo1.resize(Coordinate(1, 1, 1, 1), Coordinate(1, 1, 1, 1));}
    geo_ext = geo1;
    const CommPlan& plan = get_comm_plan(set_marks_field_1, "", geo_ext);
    //bfac  = bfac_a;
    ////Tsize = Tsize_a;

    geo_to_nv(geo, nv, Nv, mv);
    Nvol     = geo.local_volume();
    Nvol_ext = geo_ext.local_volume_expanded();
    //Nv.resize(4);nv.resize(4);mv.resize(4);
    //for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    //for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
    if(!smear_in_time_dir){dirL = 3;}
    if( smear_in_time_dir){dirL = 4;}

    get_mapvq(plan.send_pack_infos, mapvq_send, 0);
    get_mapvq(plan.recv_pack_infos, mapvq_recv, 1);
    send_reqs.resize(plan.send_msg_infos.size());
    recv_reqs.resize(plan.recv_msg_infos.size());


    /////const int dir_max = 4;
    map_bufV.resize(2);
    for(unsigned int i=0;i<map_bufV.size();i++){map_bufV[i].resize(Nvol       );}
    #pragma omp parallel for
    for(unsigned long index=0;index<Nvol;index++){
      const Coordinate xl = geo.coordinate_from_index(index);
      map_bufV[0][index]  = geo_ext.offset_from_coordinate(xl);
      map_bufV[1][index]  = index;
    }

    map_bufD.resize(Nvol*dirL*2);
    for(int dir=-dirL;dir<dirL; dir++)
    #pragma omp parallel for
    for(unsigned long index=0;index<Nvol;index++)
    {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      map_bufD[index*dirL*2 + (dir+dirL)] = geo_ext.offset_from_coordinate(xl1);
    }
    mem_setup_flag = true;
  }

  /////define the redistributed smearing kernels
  void init_distribute()
  {
    TIMERA("Gauge init_distribute");
    qassert(dirL==3);
    qassert(gauge_setup_flag);
    if(fft_copy == 1){return ;}

    ///////Redistribute the data
    ////if(int(NVmpi) < bfac){return ;}

    check_setup();
    //groupP = (bfac+NVmpi-1)/NVmpi;
    //print0("====Vec setup, NVmpi %d, groupP %d \n", NVmpi, groupP);

    vec_rot = new Vec_redistribute(fd);
    /////update gf to each MPI core
    qlat::vector_gpu<Ty > gfT;gfT.resize(NVmpi*gauge.size());
    qlat::vector_gpu<Ty > gfT_buf;gfT_buf.resize(NVmpi*gauge.size());
    const int dir_max = 4;
    const size_t Ndata = gauge.size();

    for(long vi=0;vi<NVmpi;vi++){cpy_data_thread(&(gfT.data()[vi*Ndata]), gauge.data(), Ndata);}

    vec_rot->reorder(gfT.data(), gfT_buf.data(), 1, (dir_max*2)*9 ,   0);
    gauge.copy_from(gfT);

    //Ty tmp = gauge.norm();
    //print0("v %.3e %.3e \n" , tmp.real(), tmp.imag());
    //for(long vi=0;vi<NVmpi;vi++){
    //  Ty tmp = gfT.data()[vi*Ndata];
    //  printf("v %.3e %.3e \n" , tmp.real(), tmp.imag());
    //}
    /////update gf to each MPI core


    ////qassert(false);///Need match with CPU
    ///////T distributed, spatial z,y,x
    //Nv.resize(4);nv.resize(4);mv.resize(4);
    //for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    //for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

    Nvol     = Nv[3] * nv[2]*nv[1]*nv[0];
    Nvol_ext = Nv[3] * nv[2]*nv[1]*nv[0];

    std::vector<int > Nn;Nn.resize(4);
    for(int i=0;i<3;i++){Nn[i] = nv[i];}
    Nn[3] = Nv[3];

    //std::vector<qlat::vector_acc<long > > map_bufV;
    map_bufV.resize(2);
    for(unsigned int i=0;i<map_bufV.size();i++){map_bufV[i].resize(Nvol       );}
    #pragma omp parallel for
    for(unsigned long index=0;index<Nvol;index++){
      map_bufV[0][index] = index;
      map_bufV[1][index] = index;
    }

    map_bufD.resize(Nvol*dirL*2);
    for(int dir=-dirL;dir<dirL; dir++)
    #pragma omp parallel for
    for(unsigned long index=0;index<Nvol;index++)
    {
      std::vector<int > xl;xl.resize(4);
      xl[3] =  index/(Nn[2]*Nn[1]*Nn[0]);
      xl[2] = (index%(Nn[2]*Nn[1]*Nn[0]))/(Nn[1]*Nn[0]);
      xl[1] = (index/(Nn[0]))%Nn[1];
      xl[0] =  index%(Nn[0]);

      int di = 0;
      if(dir >= 0){di=dir   ;xl[di] = (xl[di]+Nn[di]+1)%Nn[di];}
      if(dir <  0){di=-dir-1;xl[di] = (xl[di]+Nn[di]-1)%Nn[di];}

      //int ti =  index/(nv[2]*nv[1]*nv[0]);
      //int zi = (index%(nv[2]*nv[1]*nv[0]))/(nv[1]*nv[0]);
      //int yi = (index/(nv[0]))%nv[1];
      //int xi =  index%(nv[0]);
      //const Coordinate xl(xi,yi,zi,ti);
      //const Coordinate xl1 = coordinate_shifts(xl, dir);

      ///map_bufD[(dir+4)*Nvol + index] = ((xl[3]*Nn[2]+xl[2])*Nn[1] + xl[1])*Nn[0] + xl[0];
      //map_bufD[index*dirL*2 + (dir+dirL)] = geo_ext.offset_from_coordinate(xl1);
      map_bufD[index*dirL*2 + (dir+dirL)] =  ((xl[3]*Nn[2]+xl[2])*Nn[1] + xl[1])*Nn[0] + xl[0];
    }

    //geo_ext.resize(Coordinate(0, 0, 0, 0), Coordinate(0, 0, 0, 0));

    ///////fd update for box smearing
    desc_xyz_in_one(fd_new, geo);

    fft_copy = 1;
  }

  void refresh_expanded_GPU(Ty* f, const int bfac, int GPU = 1)
  {
    check_setup();
    if(fft_copy == 1){return ;}
    const CommPlan& plan = get_comm_plan(set_marks_field_1, "", geo_ext);
    const long total_bytes =
        (plan.total_recv_size + plan.total_send_size) * bfac * sizeof(Ty);
    if (0 == total_bytes) {return;}
    TIMER_FLOPS("refresh_expanded");
    timer.flops += total_bytes / 2;

    send_buffer.resize(plan.total_send_size*bfac);
    recv_buffer.resize(plan.total_recv_size*bfac);
    ////init_buf<T >(bfac);
    //T* send_buffer = (T*) send_bufferV;
    //T* recv_buffer = (T*) recv_bufferV;
    ////shift_copy(send_buffer, f , smf.mapvq_send, smf.bfac);
    cpy_data_from_index(send_buffer.data(), f , mapvq_send[0].data(), mapvq_send[1].data(), mapvq_send[0].size(), bfac, GPU);

    /////vector<MPI_Request> send_reqs(plan.send_msg_infos.size());
    /////vector<MPI_Request> recv_reqs(plan.recv_msg_infos.size());
    {
      //sync_node();
      TIMER_FLOPS("refresh_expanded-comm");
      timer.flops +=
          (plan.total_recv_size + plan.total_send_size) * bfac * sizeof(Ty) / 2;
      {
        TIMER("refresh_expanded-comm-init");
        const int mpi_tag = 10;
        for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
          const CommMsgInfo& cmi = plan.recv_msg_infos[i];
          MPI_Irecv(&recv_buffer[cmi.buffer_idx*bfac], cmi.size * bfac * sizeof(Ty), MPI_BYTE,
                    cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);
        }
        for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
          const CommMsgInfo& cmi = plan.send_msg_infos[i];
          MPI_Isend(&send_buffer[cmi.buffer_idx*bfac], cmi.size * bfac * sizeof(Ty), MPI_BYTE,
                    cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);
        }
      }
      MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
      MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
      //sync_node();
    }
    //shift_copy(f, recv_buffer , mapvq_recv, bfac);
    cpy_data_from_index(f, recv_buffer.data(), mapvq_recv[0].data(), mapvq_recv[1].data(), mapvq_recv[0].size(), bfac, GPU);
  }

  //void gauge_setup(qlat::vector_gpu<T >& gfE, const GaugeFieldT<Tg >& gf,
  //           const qlat::vector_acc<qlat::Complex >& momF = qlat::vector_acc<qlat::Complex >(), bool force_update = false);
  template <class Tg>
  void gauge_setup(const GaugeFieldT<Tg >& gf, const CoordinateD& mom_, bool force_update = false)
  {
    qlat::vector_acc<qlat::Complex > momF(8);
    for (int i = 0; i < 8; ++i) {
      const int dir = i - 4;
      const double phase = dir >= 0 ? mom_[dir] : -mom_[-dir - 1];
      momF[i] = std::polar(1.0, -phase);
    }

    bool update = false;
    if(gauge_setup_flag == false){ update = true;}
    if(force_update){ update = true;}
    if(gauge.size() == 0){  update = true;}
    if(gauge_check != (void*) qlat::get_data(gf).data()){update = true;}
    if(mom_factors.size() != 8){update = true;}
    else{for(int i=0;i<momF.size();i++){if(momF[i]!=mom_factors[i]){update = true;}}}

    if(update){
      TIMERA("gauge setup");
      qassert(geo.total_site() == gf.geo().total_site());
      gauge.resize(2*4* gf.geo().local_volume() * 9);
      mom_factors = momF;
      extend_links_to_vecs(gauge.data(), gf, momF);
      gauge_setup_flag = true;
      gauge_check = (void*) qlat::get_data(gf).data();

      ///Need to redistribute when copied
      fft_copy = 0 ;
    }
  }

  void clear_mem(){
    #ifdef __CLEAR_SMEAR_MEM__
    prop.resize(0);
    prop_buf.resize(0);
    for(int i=0;i<vL.size();i++){vL[i].resize(0);}
    mv_idx.free_mem();
    #endif
  }

  ~smear_fun(){
    clear_mem();
    if(vec_rot != NULL){delete vec_rot; vec_rot = NULL;}
  }

};

/////smear plan buffer related
struct SmearPlanKey {
  //Geometry geo;
  Coordinate total_site;
  bool smear_in_time_dir;
  int bfac;
  int civ;
  DATA_TYPE prec;
};

inline bool operator<(const SmearPlanKey& x, const SmearPlanKey& y)
{
  if(x.total_site < y.total_site ){  return true;}
  if(y.total_site < x.total_site ){  return false;}

  if(x.smear_in_time_dir < y.smear_in_time_dir ){  return true;}
  if(y.smear_in_time_dir < x.smear_in_time_dir ){  return false;}

  if(x.bfac < y.bfac ){return true;}
  if(y.bfac < x.bfac ){return false;}

  if(x.civ  < y.civ  ){return true;}
  if(y.civ  < x.civ  ){return false;}

  if(x.prec < y.prec ){return true;}
  if(y.prec < x.prec ){return false;}
  return false;
}

struct smear_fun_copy{
  DATA_TYPE prec;
  void* smfP;
  bool is_copy;  // do not free memory if is_copy=true

  smear_fun_copy(){smfP = NULL;is_copy = false;prec = Complex_TYPE;}
  smear_fun_copy(const smear_fun_copy& smf)
  {
    #ifndef QLAT_USE_ACC
    qassert(false);
    #endif
    is_copy = true;
    smfP = smf.smfP;
    prec = smf.prec;
  }
  smear_fun_copy(smear_fun_copy&& smf) noexcept
  {
    is_copy = smf.is_copy;
    smfP = smf.smfP;
    prec = smf.prec;
    smf.is_copy = true;
  }

  qacc void swap(smear_fun_copy& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    void* tmp = smfP;
    smfP   = x.smfP;
    x.smfP = tmp;
    DATA_TYPE tmp_prec = prec;
    prec   = x.prec;
    x.prec = tmp_prec;
  }
  //
  const smear_fun_copy& operator=(const smear_fun_copy& smf)
  {
    qassert(not is_copy);
    delete_pointer();
    prec = smf.prec;

    //smear_fun<T > smf;
    //smf.setup(gf, mom, smear_in_time_dir);
    if(smf.smfP == NULL){abort_r("smear fun point NULL!\n");}
    ////print0("Smear fun copy \n");
    //smf.is_copy = true;
    //smfP = smf.smfP;

    if(prec == Complex_TYPE ){
      const smear_fun<Complex  >* a =   ((const smear_fun<Complex  >*) smf.smfP);
      smfP = (void*) (new smear_fun<Complex  >(a->geo, a->smear_in_time_dir));
    }
    else if(prec == ComplexF_TYPE ){
      const smear_fun<ComplexF >* a =   ((const smear_fun<ComplexF >*) smf.smfP);
      smfP = (void*) (new smear_fun<ComplexF >(a->geo, a->smear_in_time_dir));
    }
    else{print0("Only Complex and ComplexF supported for smearing! \n");qassert(false);}

    /////swap(*this, smf);
    return *this;
  }

  void delete_pointer()
  {
    if(smfP != NULL){
    if(prec == Complex_TYPE  ){delete ((smear_fun<Complex  >*)smfP);}
    if(prec == ComplexF_TYPE ){delete ((smear_fun<ComplexF >*)smfP);}
    smfP = NULL;
    }
  }
  ~smear_fun_copy(){if(is_copy == false){delete_pointer();}}

};


inline smear_fun_copy make_smear_plan(const Geometry& geo, const bool smear_in_time_dir, const DATA_TYPE prec)
{
  TIMER_VERBOSE("make_smear_plan");
  smear_fun_copy st;st.prec = prec;
  if(prec == Complex_TYPE ){     st.smfP = (void*) (new smear_fun<Complex >(geo, smear_in_time_dir));}
  else if(prec == ComplexF_TYPE){st.smfP = (void*) (new smear_fun<ComplexF>(geo, smear_in_time_dir));}
  else{print0("Only Complex and ComplexF supported for smearing! \n");qassert(false);}
  return st;
}

inline smear_fun_copy make_smear_plan(const SmearPlanKey& skey)
{
  Geometry geo;geo.init(skey.total_site, 1);
  return make_smear_plan(geo, skey.smear_in_time_dir, skey.prec);
}

inline Cache<SmearPlanKey, smear_fun_copy >& get_smear_plan_cache()
{
  static Cache<SmearPlanKey, smear_fun_copy > cache("SmearPlanKeyCache", 16);
  return cache;
}

inline const smear_fun_copy& get_smear_plan(const SmearPlanKey& skey)
{
  if (!get_smear_plan_cache().has(skey)) {
    get_smear_plan_cache()[skey] = make_smear_plan(skey);
  }
  ////smear_fun_copy& = get_smear_plan_cache()[skey];
  ////print0("setup %5d %5d \n", skey.civ, skey.bfac);
  return get_smear_plan_cache()[skey];
}

template<typename Ty, int bfac, int civ>
inline SmearPlanKey get_smear_plan_key(const Geometry& geo, const bool smear_in_time_dir)
{
  SmearPlanKey skey;
  //skey.geo  = geo.total_site();
  skey.total_site  = geo.total_site();
  skey.bfac = bfac;
  skey.civ  = civ;
  skey.smear_in_time_dir = smear_in_time_dir;
  skey.prec = get_data_type<Ty >();
  return skey;
}
/////smear plan buffer related

////TODO need to change other parts for c0, c1
template <class T, class Tg>
void extend_links_to_vecs(T* resE, const GaugeFieldT<Tg >& gf, const qlat::vector_acc<qlat::Complex >& mom_factors=qlat::vector_acc<qlat::Complex >()){
  TIMERB("extend_links_to_vecs");
  const Geometry& geo = gf.geo();
  GaugeFieldT<Tg > gf1;
  set_left_expanded_gauge_field(gf1, gf);
  const int dir_limit = 4;
  ////set up mom factors
  qlat::Complex* momF = NULL;
  qlat::vector_acc<qlat::Complex > buf_mom;buf_mom.resize(8);
  for(int i=0;i<buf_mom.size();i++){buf_mom[i] = 1;}
  if(mom_factors.size() == 0){momF = (qlat::Complex*) qlat::get_data(buf_mom).data();}
  if(mom_factors.size() != 0){
    qassert(mom_factors.size()==8);
    momF = (qlat::Complex*) qlat::get_data(mom_factors).data();
  }
  ////set up mom factors
  //mom_factors()[dir + 4];
  qacc_for(index,  geo.local_volume(),{
    for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      const Coordinate xl = geo.coordinate_from_index(index);
      const ColorMatrixT<Tg > link =
          dir >= 0 ? gf1.get_elem(xl, dir)
                   : (ColorMatrixT<Tg >)matrix_adjoint(
                         gf1.get_elem(coordinate_shifts(xl, dir), -dir - 1));
      ////convention used in gauss_smear_kernel CPU and gauss_smear_global4
      ////also in multiply_gauge of utils_shift_vecs.h
      for(int ci=0; ci<9; ci++){
        //resE[(index*dir_limit*2+ (dir+dir_limit))*9 + (ci%3)*3 + ci/3] = momF[dir + 4] * link.p[ci];
        resE[(index*dir_limit*2+ (dir+dir_limit))*9 + ci] = momF[dir + 4] * link.p[ci];
      }
    }
  });
}

template <class T, class Tg>
void extend_links_to_vecs(qlat::vector_gpu<T >& gfE, const GaugeFieldT<Tg >& gf, const qlat::vector_acc<qlat::Complex >& mom_factors=qlat::vector_acc<qlat::Complex >()){
  gfE.resize(2*4* gf.geo().local_volume() * 9);
  extend_links_to_vecs(gfE.data(), gf, mom_factors);
}

template <class T>
void rotate_Vec_prop(Propagator4dT<T>& prop, qlat::vector_gpu<T > &propT, unsigned int NVmpi, unsigned int groupP, int dir = 0)
{
  TIMER("Rotate Vec prop");
  ///unsigned int NVmpi = fd.mz*fd.my*fd.mx;
  ///int groupP = (12+NVmpi-1)/NVmpi;
  qassert(groupP <= 12);qassert(groupP * NVmpi >= 12);
  qassert(groupP >   0);

  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();
  if(dir == 0)propT.resize(NVmpi*Nvol*groupP*12);

  qacc_for(index, long(Nvol), {
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem(index);

    for(int c1 = 0;c1 < 3; c1++)
    for(int d1 = 0;d1 < 4; d1++)
    for(int c0 = 0;c0 < 3; c0++)
    for(int d0 = 0;d0 < 4; d0++)
    {
      int off1 = c1*4+d1;
      int n0 = off1/groupP;
      int n1 = off1%groupP;

      int off0 = c0*4+d0;
      //LInt off = (c0*3+c1)*16+d0*4 + d1;
      //long offP = n0*Nvol*groupP*12 + index*groupP*12 + n1*12 + off0;
      long offP = ((n0*Nvol + index)*groupP + n1)*12 + off0;
      if(dir == 0){propT[offP] = v0(d0*3+c0, d1*3+c1);}
      if(dir == 1){v0(d0*3+c0, d1*3+c1) = propT[offP];}
    }
  });


}

void get_smear_para(std::string smear_para, double& width, int& step)
{
  if(smear_para == std::string("NONE")){width = 0.0;step = 0;return; }
  std::vector<std::string > temL = stringtolist(smear_para);
  if(temL.size() != 2){abort_r("read smear para wrong! \n");}
  width = stringtodouble(temL[0]);
  step  = stringtonum(temL[1]);
}

template <class T, int bfac, int civ>
void smear_propagator_box_kernel(qlat::vector_gpu<T >& prop, qlat::vector_gpu<T >& vp, qlat::vector_gpu<T >& vm,
  const int Bsize, const int dir, shift_vec& svec)
{
  //return ;
  TIMERC("Shift Kernel");
  std::vector<int > iDir(4);
  for(int i=0;i<4;i++){iDir[i] = 0;}

  vp.copy_from(prop);vm.copy_from(prop);
  for(int bi=0;bi<Bsize;bi++){
    iDir[dir] =  1;svec.shift_vecP(vp.p, vp.p, iDir , bfac*3*civ);
    iDir[dir] = -1;svec.shift_vecP(vm.p, vm.p, iDir , bfac*3*civ);
    prop += vp;
    prop += vm;
  }
}

template <class T, int bfac, int civ>
void smear_propagator_box(T* src, const int Bsize, smear_fun<T >& smf){
  if (Bsize <= 0) {return;}
  /////qassert(!smear_in_time_dir);
  smf.check_setup();
  long Nvol = smf.Nvol;
  //const Geometry& geo = smf.geo;
  //const int dir_limit = smear_in_time_dir ? 4 : 3;
  const int nsites = bfac*3*civ;
  const int dir_limit = 3;
  /////std::vector<int > nv,Nv,mv;geo_to_nv(geo, nv, Nv, mv);

  //qlat::vector_gpu<T > gfE; 
  //extend_links_to_vecs(gfE, gf);
  ////fft_desc_basic fd(geo);
  shift_vec svec(smf.fd_new, true);
  svec.set_gauge(smf.gauge.data(), bfac, civ);
  //rotate_prop(prop,0);

  size_t Ndata = size_t(Nvol) * nsites;
  ////T* pdata = (T*) qlat::get_data(prop).data();

  for(unsigned int i=0;i<smf.vL.size();i++){smf.vL[i].resize(Ndata);smf.vL[i].set_zero();}
  qlat::vector_gpu<T >& vprop = smf.vL[0];
  qlat::vector_gpu<T >& v1 = smf.vL[1];
  qlat::vector_gpu<T >& v2 = smf.vL[2];
  qlat::vector_gpu<T >& vp = smf.vL[3];
  qlat::vector_gpu<T >& vm = smf.vL[4];
  std::vector<qlat::vector_gpu<T >*  > vL(3);
  for(unsigned int i=0;i<vL.size();i++){vL[i] = &smf.vL[i + 5];}

  //qlat::vector_gpu<T > vprop; vprop.resize(Ndata);
  //qlat::vector_gpu<T > v1; v1.resize(Ndata);
  //qlat::vector_gpu<T > v2; v2.resize(Ndata);
  //qlat::vector_gpu<T > vp; vp.resize(Ndata);
  //qlat::vector_gpu<T > vm; vm.resize(Ndata);
  //std::vector<qlat::vector_gpu<T >  > vL(3);
  //for(unsigned int i=0;i<vL.size();i++){vL[i].resize(Ndata);vL[i].set_zero();}

  cpy_data_thread(vprop.p, src, Ndata);

  for(int dir = 0; dir < dir_limit; ++dir)
  {
    v1.copy_from(vprop);
    smear_propagator_box_kernel<T,bfac,civ>(v1, vp, vm, Bsize, dir, svec);
    int d1=(dir+1)%3;int d2=(dir+2)%3;v2.copy_from(v1);
    smear_propagator_box_kernel<T,bfac,civ>(v1, vp, vm, Bsize, d1 , svec);
    smear_propagator_box_kernel<T,bfac,civ>(v2, vp, vm, Bsize, d2 , svec);
    (*vL[d1]) += v2;
    (*vL[d2]) += v1;
  }
  vprop.set_zero();
  for(int dir = 0; dir < dir_limit; ++dir)
  {
    smear_propagator_box_kernel<T,bfac,civ>(*vL[dir], vp, vm, Bsize, dir, svec);
    vprop += (*vL[dir]);
  }

  T* data = vprop.p;
  qacc_for(index, Nvol, {
    for(int c=0;c<nsites;c++){src[index*nsites + c] = data[index*nsites + c] * (1.0/6.0);} 
  });

  smf.clear_mem();
}

template <class T, int bfac, int civ>
void gauss_smear_kernel(T* src, const double width, const int step, const T norm, smear_fun<T >& smf)
{
  double bw_fac = width*width/(4.0*step - 6.0*width*width);
  ////match convention to qlat
  if(smf.dirL == 4 ){bw_fac = bw_fac * 3.0/4.0;}
  const T bw = bw_fac;

  unsigned long Nvol = smf.Nvol;
  const int nsites = bfac*3*civ;

  smf.check_setup();
  T* gf = (T*) smf.gauge.data();

  T* prop     = NULL;
  T* prop_buf = NULL;
  ////int dirL = 3;
  #ifdef QLAT_USE_ACC
  const int GPU = 1;
  int nt  = 3*3*9;
  if(bfac*civ <= 12){ nt =         32;}
  if(bfac*civ <=  6){ nt = 3*bfac*civ;}
  //int nt = 32;
  dim3 dimBlock(nt, 1, 1);
  long sn = Nvol;
  dim3 dimGrid( sn, 1, 1);
  if(smf.fft_copy==0){
    smf.prop.resize(    smf.Nvol * bfac * 3* civ );
    smf.prop_buf.resize(smf.Nvol_ext * bfac * 3* civ );
  }
  prop = smf.prop.data();
  prop_buf = smf.prop_buf.data();
  if(smf.fft_copy==0){cpy_data_thread(prop, src, smf.prop.size());}
  #else
  const int GPU = 0;
  if(smf.fft_copy==0){smf.prop_buf.resize(smf.Nvol_ext * bfac * 3* civ );}
  prop = src;
  prop_buf = smf.prop_buf.data();
  #endif

  long* Pdir1 = (long*) qlat::get_data(smf.map_bufD).data();

  for (int i = 0; i < step; ++i) {
    {TIMERC("Copy prop");
    cpy_data_from_index(prop_buf, prop, &smf.map_bufV[0][0], &smf.map_bufV[1][0], Nvol, nsites, GPU);}

    {TIMERC("Communication");
    smf.refresh_expanded_GPU(prop_buf, nsites, GPU);}

    {TIMERC("Matrix multiply");
    #ifdef QLAT_USE_ACC
    if(smf.dirL==3){gauss_smear_global4<T, bfac, civ, 3><<< dimGrid, dimBlock >>>(prop, prop_buf, gf, bw, norm, Nvol, Pdir1);}
    if(smf.dirL==4){gauss_smear_global4<T, bfac, civ, 4><<< dimGrid, dimBlock >>>(prop, prop_buf, gf, bw, norm, Nvol, Pdir1);}
    #else

    const int dir_max = 4;
    const int dir_limit = smf.dirL;
    qacc_for(index,  long(Nvol),{
      T buf[nsites];for(int i=0;i<nsites;i++){buf[i] = 0;}
      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
        const T* wm1p = &prop_buf[size_t(Pdir1[index*dir_limit*2 + (dir + dir_limit)])*nsites];
        const T* lp = &gf[(index*dir_max*2 + dir + dir_max)*9];

        for(int bi=0;bi<bfac;bi++){
          Eigen::Matrix<T, 3, 3, Eigen::ColMajor>&     lE = *((Eigen::Matrix<T, 3, 3, Eigen::ColMajor>*) lp);
          Eigen::Matrix<T, civ, 3, Eigen::ColMajor>& wmE = *((Eigen::Matrix<T, civ, 3, Eigen::ColMajor>*)  &buf[bi*3*civ]);
          Eigen::Matrix<T, civ, 3, Eigen::ColMajor>&  pE = *((Eigen::Matrix<T, civ, 3, Eigen::ColMajor>*) &wm1p[bi*3*civ]);
          wmE += ( pE * lE);////avoid mix of Col and Row when civ == 1
        }
      }
      T* wmp = &prop[index*nsites];
      //for(int i=0;i<nsites;i++){wmp[i] += bw*buf[i];}
      for(int i=0;i<nsites;i++){wmp[i]  = norm*(wmp[i] + bw*buf[i]);}
    });
    #endif
    }
  }

  #ifdef QLAT_USE_ACC
  if(smf.fft_copy==0){cpy_data_thread(src, prop, smf.prop.size());}
  #endif
  smf.clear_mem();
}

template <class T, int bfac, int civ>
void smear_kernel_sort(T* src, const double width, const int step, smear_fun<T >& smf)
{
  if(step >= 0){
    const double norm   = (1 - 3.0*width*width/(2*step));
    gauss_smear_kernel<T, bfac, civ>(src, width, step, norm, smf);
  }

  if(step ==-1){
    if(smf.smear_in_time_dir == true){abort_r("Box smearing not supported for t smearings \n");}
    smear_propagator_box<T, bfac, civ>(src,int(width + 0.001), smf);
  }
}

template <class T>
void smear_kernel(T* src, const double width, const int step, smear_fun<T >& smf, int bfac, int civ)
{
  bool cfind = false;
  {
  TIMERC("==smear_kernel");
  #define smear_macros(ba,da) if(bfac == ba and civ ==  da and cfind == false){cfind = true; \
    smear_kernel_sort<T, ba, da >(src, width, step, smf);}

  /////max size of a prop
  ///////macros for color 3 * dirac 4 in inner prop
  smear_macros(   1,  4);
  smear_macros(   2,  4);
  smear_macros(   3,  4);
  smear_macros(   4,  4);
  smear_macros(   5,  4);
  smear_macros(   6,  4);
  smear_macros(   7,  4);
  smear_macros(   8,  4);
  smear_macros(   9,  4);
  smear_macros(  10,  4);
  smear_macros(  11,  4);
  smear_macros(  12,  4);

  ///////macros for inner color 3 and all outter prop
  smear_macros(   1,  1);
  smear_macros(   2,  1);
  smear_macros(   3,  1);
  smear_macros(   4,  1);
  smear_macros(   5,  1);
  smear_macros(   6,  1);
  smear_macros(   7,  1);
  smear_macros(   8,  1);
  smear_macros(   9,  1);
  smear_macros(  10,  1);
  smear_macros(  11,  1);
  smear_macros(  12,  1);

  smear_macros(  24,  1);
  smear_macros(  48,  1);

  smear_macros(   1, 48);
  smear_macros(   1, 24);
  smear_macros(   1, 18);
  smear_macros(   1, 16);
  smear_macros(   1, 12);
  smear_macros(   1,  9);
  smear_macros(   1,  8);
  smear_macros(   1,  6);
  smear_macros(   1,  3);
  smear_macros(   1,  2);

  smear_macros(   4, 12);
  smear_macros(   4,  3);
  smear_macros(   4,  6);
  smear_macros(   4,  8);
  smear_macros(   4,  2);
  smear_macros(   3, 16);

  qacc_barrier(dummy);
  #undef smear_macros
  }
  if(cfind == false){print0("Case bfac %5d, civ %5d \n", bfac, civ); abort_r("smear kernel not found!\n");}
  ////qassert(cfind);
}

template <class T>
void rotate_prop(Propagator4dT<T>& prop, int dir = 0)
{
  TIMERB("rotate_prop");
  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();

  T* src =  (T*) qlat::get_data(prop).data();
  qacc_for(index, long(Nvol), {
    T buf[12*12];
    T* res = &src[index*12*12];
    for(int i=0;i<12*12;i++){buf[i] = res[i];}

    for(int d1 = 0;d1 < 4; d1++)
    for(int c1 = 0;c1 < 3; c1++)
    for(int d0 = 0;d0 < 4; d0++)
    for(int c0 = 0;c0 < 3; c0++)
    {
      int off0 = c1*4*4*3 + ((d1*4+d0)*3+c0);
      int off1 = (d1*3+c1)*12 + d0*3 + c0;
      if(dir == 0){res[off0] = buf[off1];}
      if(dir == 1){res[off1] = buf[off0];}
    }
  });
}


template <class Ty, int c0,int d0, class Tg>
void smear_propagator_gwu_convension_inner(Ty* prop, const GaugeFieldT<Tg >& gf,
                      const double width, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const int mode = 1)
{
  TIMER_FLOPS("==smear propagator");
  long long Tfloat = 0;
  ///double mem       = 0.0;
  Geometry geo = gf.geo();
  geo.multiplicity = 1;geo.eo=0;
  if(c0*d0 > 48){abort_r("d0 should be smaller than 48 for gpu mem. \n");}

  {
    long long Lat = geo.local_volume();
    int nsrc = c0*d0;
    long long vGb = Lat *nsrc;
    const int n_avg = smear_in_time_dir ? 8 : 6;
    int Fcount = 3*(3*6 + 2*2); 
    if(step >= 0){
    Tfloat = step*n_avg*vGb*Fcount;
    }else{
    Tfloat = c0*d0*2*int(width)*vGb*Fcount;
    }
  }
  timer.flops += Tfloat;

  if (0 == step) {return;}
  /////smear_fun<T > smf(gf.geo(), smear_in_time_dir);

  Ty* src = prop;
  //Ty* src = (Ty*) qlat::get_data(prop).data();
  /////rotate_prop(prop, 0);

  SmearPlanKey skey = get_smear_plan_key<Ty, c0, d0>(geo, smear_in_time_dir);
  const smear_fun_copy& smf_copy = get_smear_plan(skey);
  smear_fun<Ty >& smf = *((smear_fun<Ty >*) (smf_copy.smfP));
  smf.gauge_setup(gf, mom);

  fft_desc_basic fd(geo);
  ////Vec_redistribute vec_large(fd);

  long Nvol_pre = geo.local_volume();
  bool reorder = false;
  int groupP = (d0+smf.NVmpi-1)/smf.NVmpi;
  if(smf.NVmpi <= d0 and smf.NVmpi != 1 and (d0%smf.NVmpi == 0) and smear_in_time_dir == false and mode == 1){reorder = true;}
  if(reorder ){
    int Nprop = smf.NVmpi * c0*3 * groupP;
    /////print0("====Vec setup, NVmpi %d, groupP %d \n", smf.NVmpi, groupP);
    smf.init_distribute();
    smf.prop.resize(    Nvol_pre * Nprop );
    smf.prop_buf.resize(Nvol_pre * Nprop );
    qassert(!smear_in_time_dir);

    Ty* res = smf.prop.data();
    cpy_data_thread(res, src, smf.prop.size());
    smf.mv_idx.dojob(res, res, 1, smf.NVmpi, Nvol_pre*c0*3, 1,   groupP, true);
    {TIMERC("Vec prop");smf.vec_rot->reorder(res, smf.prop_buf.data(), 1, c0*3*groupP ,   0);}
    src = smf.prop.data();
  }

  ////print0("===Case %d %d \n", c0, grouP);
  if( reorder)smear_kernel(src, width, step, smf,  c0, groupP);
  if(!reorder)smear_kernel(src, width, step, smf,  c0, d0);

  if(reorder ){
    Ty* res = smf.prop.data();
    {TIMERC("Vec prop");smf.vec_rot->reorder(res, smf.prop_buf.data(), 1, c0*3*groupP , 100);}
    smf.mv_idx.dojob(res, res, 1, smf.NVmpi, Nvol_pre*c0*3, 0,  groupP , true);
    ////src = (Ty*) qlat::get_data(prop).data();
    ////src = prop;
    cpy_data_thread(prop, smf.prop.data(), smf.prop.size());
  }
  /////rotate_prop(prop, 1);
}

template <class Ty, class Tg>
void smear_propagator_gwu_convension(qpropT& prop, const GaugeFieldT<Tg >& gf,
                      const double width, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const int mode = 1)
{
  if (0 == step) {return;}
  const long Nvol = prop.geo().local_volume();
  Ty* src = (Ty*) qlat::get_data(prop).data();
  move_index mv_civ;int flag = 0;
  flag = 0;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, false);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(int d0=0;d0<12*4;d0++)
    for(int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + c0*12*4 + d0] = buf[d0*3 + c0];
    }
  });

  smear_propagator_gwu_convension_inner<Ty, 1, 12*4, Tg>(src, gf, width, step, mom, smear_in_time_dir, mode);

  qacc_for(isp, Nvol, {
    Ty buf[12*12];
    for(int i=0;i<12*12;i++){buf[i] = src[isp*12*12 + i];}
    for(int d0=0;d0<12*4;d0++)
    for(int c0=0;c0<   3;c0++)
    {
      src[isp*12*12 + d0*3 + c0] = buf[c0*12*4 + d0];
    }
  });
  flag = 1;mv_civ.dojob(src, src, 1, 12*12, Nvol, flag, 1, false);
}

template <class Ty, class Tg>
void smear_propagator_gwu_convension(qlat::FieldM<Ty , 12>& prop, const GaugeFieldT<Tg >& gf,
                      const double width, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const int mode = 1)
{
  if (0 == step) {return;}
  Ty* src = (Ty*) qlat::get_data(prop).data();
  qacc_for(isp, prop.geo().local_volume(), {
    Ty buf[12];
    for(int i=0;i<12;i++){buf[i] = src[isp*12 + i];}
    for(int d0=0;d0<4;d0++)
    for(int c0=0;c0<   3;c0++)
    {
      src[isp*12 + c0*4 + d0] = buf[d0*3 + c0];
    }
  });

  smear_propagator_gwu_convension_inner<Ty, 1, 4, Tg>(src, gf, width, step, mom, smear_in_time_dir, mode);
  qacc_for(isp, prop.geo().local_volume(), {
    Ty buf[12];
    for(int i=0;i<12;i++){buf[i] = src[isp*12 + i];}
    for(int d0=0;d0<4;d0++)
    for(int c0=0;c0<   3;c0++)
    {
      src[isp*12 + d0*3 + c0] = buf[c0*4 + d0];
    }
  });
}


template <class Ty, class Tg>
void smear_propagator_gwu_convension(Propagator4dT<Ty>& prop, const GaugeFieldT<Tg >& gf,
                      const double width, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const int mode = 1)
{
  if (0 == step) {return;}
  rotate_prop(prop, 0);
  Ty* src = (Ty*) qlat::get_data(prop).data();
  smear_propagator_gwu_convension_inner<Ty, 1, 12*4, Tg>(src, gf, width, step, mom, smear_in_time_dir, mode);
  rotate_prop(prop, 1);
}

template <class T, class Tg>
void smear_propagator_qlat_convension(Propagator4dT<T>& prop, const GaugeFieldT<Tg >& gf,
                      const double coef, const int step, const CoordinateD& mom = CoordinateD(), const bool smear_in_time_dir = false, const int mode = 1)
{
  if(coef <= 0){return ;}
  double width = 0.0;
  if(step <  0){width = coef;}////box smearings
  if(step >= 0){width = std::sqrt(coef*2*step/(3.0));}////gauss smearings
  smear_propagator_gwu_convension(prop, gf, width, step, mom, smear_in_time_dir, mode);
}




}
#endif
