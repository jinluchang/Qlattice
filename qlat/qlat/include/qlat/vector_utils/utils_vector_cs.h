// utils_vector_cs.h
// Gen Wang
// Jun. 2023

#ifndef UTILS_VECTOR_CS_H
#define UTILS_VECTOR_CS_H

#pragma once
#include <qlat/qcd.h>
#include "utils_float_type.h"
#include "utils_COPY_data.h"

////needed for norm calculation
#include "utils_reduce_vec.h"
#include "utils_vector_GPU.h"
#include "general_funs.h"
#include "utils_fft_desc.h"
#include "utils_Matrix_prod.h"

namespace qlat{

inline Long setup_cs_b_size(LInt nsum, long long bsize0 = -1)
{
  Qassert(nsum >= 1);
  LInt bcut = 300;///default pure CPU cuts
  #ifdef QLAT_USE_ACC
  bcut = 4096;
  #endif
  if(bsize0 != -1){bcut = std::abs(bsize0);}

  std::string val = get_env(std::string("q_setup_cs_b_size"));
  if(val != ""){bcut = stringtonum(val);}

  Qassert(bcut >= 1);
  /////increase bsize_or till slightly larger than bcut
  LInt bsize_or = nsum;
  if(bsize_or > bcut)
  {
    for(LInt temb=bcut;temb<nsum;temb++)
    {
      if((nsum)%temb == 0)
      {
        bsize_or = temb;
        break;
      }
    }
  }
  return bsize_or;   ////b_size = noden*cs/bfac;
}

template <typename Ty >
struct vector_cs{
  Long nvec;
  QMEM  GPU;
  //
  ////infos
  LInt nsum;// length of each vector
  Long bfac, b_size, bfac_group;
  LInt La, Lb; // derived variables
  Long btotal; //  == bfac
  //
  std::vector<qlat::vector_gpu<Ty > > buf_g;// actual data positions
  qlat::vector<Ty* > pL; //// pointers to data
  qlat::vector_gpu<Ty > alpha_buf;
  qlat::vector<Ty* > pointersL; ////grouped with ni -> btotal 
  //
  qlat::vector_gpu<Ty > alphaG;
  qlat::vector_gpu<Ty > norm2_buf;
  //
  std::vector<Long > jobA;
  //
  bool initialized;
  //
  double flops_matrix;
  double flops_copy;
  Int work_i;
  Int work_f;
  //
  inline void inialize_para(){
    nsum   =  0;     ////sites to be summed
    nvec   =  0;     ////number of vecs
    bfac   =  0;     //volume related, outer loop, derived
    b_size = -1;     //volume related, continus in memeory
    bfac_group = -1; //create large size of data or not, for GPU memories and performance
    GPU    =  QMGPU;     ////default data on GPU,  GPU, -1 -- unified memory, 0 -- on CPU, 1 -- on GPU
    La = 0;Lb  = 0;  ////parameter from bfac, La = bfac/bfac_group, Lb = bfac_group*nvec*Long(b_size);
    btotal =  bfac;
    flops_matrix = 0.0;
    flops_copy   = 0.0;
    work_i = 0;
    work_f = 0;
    clean_mem();
  }
  //
  vector_cs(){
    inialize_para();
  }
  //
  vector_cs(Int nvec_, LInt nsum_, QMEM GPU_ = QMGPU, Int b_size_ = -1, Int bfac_group_ = -1)
  {
    resize(nvec_, nsum_, GPU_, b_size_, bfac_group_);
  }
  //
  //inline void resize(Int nvec_, Int nsum_)
  //{
  //  resize(nvec_, nsum_, 1, -1, 1, -1);
  //}
  //
  inline size_t vlength(){
    return nsum;
  }
  //
  inline void resize(Int nvec_, Int nsum_, QMEM GPU_)
  {
    resize(nvec_, nsum_, GPU_, b_size, bfac_group);
  }
  //
  inline void resize(Int nvec_, Int nsum_)
  {
    resize(nvec_, nsum_, GPU, b_size, bfac_group);
  }
  //
  inline void resize(Int nvec_)
  {
    Qassert(nsum > 0 and b_size > 0 and bfac > 0);
    resize(nvec_, nsum, GPU, b_size, bfac_group);
  }
  //
  /////==added
  template <class Ta >
  inline void resize(Int nvec_, vector_cs<Ta >& ref)
  {
    resize(nvec_, ref.nsum, ref.GPU, ref.b_size, ref.bfac_group);
  }
  //
  /////==added
  template <class Ta >
  inline void resize(Int nvec_, QMEM GPU, vector_cs<Ta >& ref)
  {
    resize(nvec_, ref.nsum, GPU, ref.b_size, ref.bfac_group);
  }
  //
  inline double get_flops_matrix(){
    double  res = flops_matrix;
    flops_matrix = 0.0;
    return res;
  }
  inline double get_flops_copy(){
    double  res = flops_copy;
    flops_copy = 0.0;
    return res;
  }
  //
  inline void resize(Int nvec_, LInt nsum_, QMEM GPU_, Int b_size_, Int bfac_group_, bool silence_mem = false)
  {
    ////Qassert(geo.node_site != Coordinate(0,0,0,0));
    //nsum = nsum_;
    //Qassert(Nt != 0 and noden != 0);
    Int flag = 0;

    if(nvec_ == nvec){
      flag += 1;
    }else{
      nvec = nvec_;
    }

    if(nsum_ == nsum){
      flag += 1;
    }else{
      nsum = nsum_;
    }

    if(GPU == GPU_){
      flag += 1;
    }else{
      GPU = GPU_;
    }

    if(b_size == b_size_){
      flag += 1;
    }else{
      b_size = b_size_;
    }

    if(bfac_group == bfac_group_){
      flag += 1;
    }else{
      bfac_group = bfac_group_;
    }

    if(b_size == -1){
      b_size = setup_cs_b_size(nsum, b_size);
    }

    Qassert(b_size > 0 and nsum > 0 and nsum % b_size == 0);
    bfac = nsum / b_size;
    Qassert(nvec >= 0 and nsum >= 1 and b_size >= 1 and bfac_group >= -1);
    if(flag == 5 and initialized){
      return ;
    }

    Long max_bfac = bfac;
    if(bfac_group <= 0){
      bfac_group = 1;
      ////more data continues on GPU
      #ifdef QLAT_USE_ACC
      double memvec = double(b_size) * sizeof(Ty) / (1024.0 * 1024.0 * 1024.0);
      bfac_group = max_bfac;
      ////if larger than 8.0 Gb
      const double max_mem = 8.0;
      if(memvec <= 0 or memvec*bfac_group > max_mem)
      {
        Long max_group = max_mem / memvec;
        for(Int bini = bfac_group;bini >= 1; bini--){
          size_t tem = get_threads(bini, max_bfac, 0);
          if(tem != 0 and max_bfac%tem == 0 and tem <= max_group){bfac_group = tem;break;}
        }
        if(max_bfac % bfac_group != 0){bfac_group = 1;}
      }
      #endif
    }
    std::string val = get_env(std::string("q_setup_cs_bfac_group"));
    if(val != ""){bfac_group = stringtonum(val);}

    if(max_bfac % bfac_group != 0){
      qmessage("total %ld, current %ld \n", long(max_bfac), long(bfac_group));
    }
    Qassert(max_bfac % bfac_group == 0);
    allocate_mem(silence_mem);
    ////qacc_barrier(dummy);
  }
  //
  inline Ty** get_pointers(Long ni)
  {
    ////qmessage("===%5d %5d \n", int(ni), int(nvec));
    Qassert(ni < nvec);
    //Ty** res = &pointersL[ni*btotal + 0]; 
    //return res; 
    return &pointersL[ni*btotal + 0];
  }
  //
  inline Ty**  get_pointers(Long ni, Long bi)
  {
    Qassert(ni < nvec and bi <= bfac and ni >=0 and bi >= 0 );
    return &pointersL[ni*btotal + bi];
  }
  //
  ////ni < nvec
  inline Ty* get_pointer_b(Long ni, Long bi)
  {
    /////qmessage("rank %5d, %5d %5d, %5d %5d %5d \n", get_id_node(), int(ni), int(nvec), int(bi), int(bfac) );
    Qassert(ni < nvec and bi < bfac and ni >=0 and bi >= 0);
    size_t t   = (bi *nvec + ni ) * size_t(b_size);
    size_t ba  = t / Lb;
    size_t bb  = t % Lb;
    return &pL[ba][bb];
  }
  //
  inline Ty* get_pointer_x(Long ni, Long xi)
  {
    Qassert(ni < nvec and ni >= 0);
    Qassert(xi < Long(nsum) and xi >= 0);
    size_t x   = xi;
    size_t t   = ((x/b_size) * nvec + ni ) * size_t(b_size) + x % b_size;
    size_t ba  = t / Lb;
    size_t bb  = t % Lb;
    /////Ty* res = &pL[ba][bb];
    return &pL[ba][bb];
  }
  //
  inline void allocate_mem(bool silence_mem = false)
  {
    TIMERA("vector cs allocate_mem");
    Qassert(bfac   >  0);
    Qassert(b_size >  0);
    Qassert(bfac_group >  0);
    jobA = job_create(bfac, bfac_group);
    if(nvec == 0){clean_mem();return ;}
    //rpV.resize(bfac_group);
    //EpV.resize(bfac_group);
    //spV.resize(bfac_group);
    //
    La = bfac/bfac_group;
    Lb = bfac_group*nvec*Long(b_size);
    if(pL.size() != 0){pL.resize(0);}
    pL.resize(La);
    //
    if(!silence_mem){
      if(buf_g.size() != La){
        buf_g.resize(0);
        buf_g.resize(La);
      }
      for(LInt i=0;i<La;i++){
        buf_g[i].resize(Lb, GPU);
        //////pL[i] = (Ty*) qlat::get_data(buf_g[i]).data();
      }
    }
    //
    for(LInt i=0;i<La;i++){
      pL[i] = (Ty*) qlat::get_data(buf_g[i]).data();
    }
    //
    btotal = bfac;
    initialized = true;
    //if(silence_mem == true){qmessage("nvec %5d, bfac %5d, size %5d %5d \n", int(nvec), int(bfac), int(b_size), int(bfac_group));}
    //
    {
      pointersL.resize(nvec * btotal);
      Ty** Pbuf = (Ty**) qlat::get_data(pointersL).data();
      Ty** pA = (Ty**) qlat::get_data(this->pL).data();
      const Long& b_size = this->b_size;
      const Long& btotal  = this->btotal;
      const Long& nvec = this->nvec;
      const unsigned long& Lb = this->Lb;
      const Int& GPU = this->GPU;
      qGPU_for(isp, btotal*nvec, GPU, {
        Long bi = isp / nvec;
        Long ni = isp % nvec;
        const size_t tf  = (bi *nvec + ni ) * size_t(b_size);
        const size_t ba  = tf / Lb;
        const size_t bb  = tf % Lb;
        Pbuf[ni*btotal + bi] = &pA[ba][bb];
      });
    }
    //
    work_f = nvec;////for projections and other operations
  }
  //
  ////==added
  inline void set_zero(Int ia = -1, QBOOL dummy = QTRUE)
  {
    TIMERA("vector_cs set_zero");
    if(!initialized){return ;}
    //bool GPU_zero = true;
    //if(GPU == 0){GPU_zero = false;}
    const Long& btotal = this->btotal;
    const Long& b_size = this->b_size;
    const Int& GPU = this->GPU;
    for(Long ni = 0; ni < nvec; ni++)
    {
      if(ia != -1 and ni != ia){ continue; }
      Ty** A =     get_pointers(ni);
      qGPU_forNB(isp, btotal*b_size, GPU, {
        const Long id = isp / b_size;
        const Long jd = isp % b_size;
        A[id][jd] = 0;
      })
    }
    if (dummy == QTRUE) { qacc_barrier(dummy); }
  }
  //
  inline void clean_mem(){
    buf_g.resize(0);
    pL.resize(0);
    //buf_V.resize(0);
    alpha_buf.resize(0);
    initialized = false;
    nvec = 0;
  }
  //
  inline double mem_report()
  {
    double memvec = double(b_size) * sizeof(Ty) / (1024.0 * 1024.0 * 1024.0);
    double total  = nvec * bfac * memvec;
    //#if PRINT_TIMER>5
    qmessage("nvec %5ld, nsum %ld, bfac %ld, b_size %ld, bfac_group %ld, %.8e GB, GPU %+1d \n", 
      nvec, nsum, bfac, b_size, bfac_group, nvec * bfac * memvec, int(GPU));
    //#endif
    return total;
  }
  //
  inline void random_cs(Int seed = 0)
  {
    Ty** res = get_pointers(0);
    for(Int bg=0;bg<btotal/bfac_group;bg++)
    {
      random_Ty(res[bg*bfac_group], bfac_group*nvec*b_size, GPU, bg + 12324 + seed, 1);
    }
  }
  //
  inline void clear(){clean_mem();}
  ~vector_cs(){
    clean_mem();
  }
  //
  ////copy from a continous memory?
  template <class Ta >
  void copy_from(Ta* src, Int ncur, bool data_GPU_ = true, QBOOL dummy = QTRUE, const Int dir = 1)
  {
    TIMER_FLOPS("vector_cs copy from A");
    if(ncur >= nvec){
      qmessage("Copy to position %8d larger than n_vec %8d ! \n", int(ncur), int(nvec));
      abort_r();
    }
    //
    QMEM data_GPU = QMGPU;if(data_GPU_ == false){data_GPU = QMCPU;}
    //
    const bool same_locate = check_GPU_same(GPU, data_GPU);
    //
    Ty** res = get_pointers(ncur);
    const Long& b_size = this->b_size;
    const Long& btotal = this->btotal;
    ////Long  Ndata = btotal * b_size;
    ////int GPU_set = this->GPU;
    const Int GPU_set = check_GPU_multi(  GPU, data_GPU);

    QMEM GPU_r = QMGPU; QMEM GPU_s = QMGPU;
    if(dir == 1){GPU_r =      GPU; GPU_s = data_GPU;}
    if(dir == 0){GPU_r = data_GPU; GPU_s =      GPU;}

    if(same_locate){
      //qGPU_for2dNB(jd, b_size, id, btotal, GPU_set, {
      //  //Long id = isp / b_size;
      //  //Long jd = isp % b_size; 
      //  const Long isp = id * b_size + jd;
      //  if(dir == 1){res[id][jd] = src[isp];}
      //  if(dir == 0){src[isp] = res[id][jd];}
      //});

      //qmessage("Check \n");
      //print_norm2(ncur);
      //print_numbers(src, 10, data_GPU);
      const Long  Ndata = btotal * b_size;
      qGPU_forNB(isp, Ndata, GPU_set, {
        const Long id = isp / b_size;
        const Long jd = isp % b_size; 
        if(dir == 1){res[id][jd] = src[isp];}
        if(dir == 0){src[isp] = res[id][jd];}
      });
      //qacc_barrier(dummy);
      //print_numbers(&res[0][0], 10, GPU);
      //{
      //print_norm2(ncur);
      //}
      //qmessage("Check \n");
    }

    if(!same_locate){
      /////if(sizeof(Ty) == sizeof(Ta))
      size_t roff = 0;size_t soff = 0;
      for(Int ig=0;ig<btotal/bfac_group;ig++)
      {
        Int id = ig * bfac_group;
        if(dir == 1){roff = nvec * size_t(b_size);soff = size_t(b_size);}
        if(dir == 0){roff = size_t(b_size);soff = nvec * size_t(b_size) ;}
        ////if(dir == 1)cpy_GPU(res[id], &src[id*b_size], b_size,     GPU, data_GPU, false);
        ////if(dir == 0)cpy_GPU(&src[id*b_size], res[id], b_size,     data_GPU, GPU, false);
        if(dir == 1)cpy_GPU2D(res[id], &src[id*b_size], 
          size_t(b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);
        if(dir == 0)cpy_GPU2D(&src[id*b_size],res[id],  
          size_t(b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);
      }
    }
    if (dummy == QTRUE) { qacc_barrier(dummy); }

    double flops = double(b_size) * btotal * sizeof(Ty);
    timer.flops  += flops;
    flops_copy += flops;
  }

  template <class Ta >
  void copy_to(Ta* src, Int ncur , bool data_GPU = true, QBOOL dummy = QTRUE)
  {
    copy_from(src, ncur, data_GPU, dummy, 0);
  }

  ////nsrc the ni's copy from the src, nres the destination ni's
  template <class Ta >
  void copy_from(vector_cs<Ta >& src, std::vector<Int >& nres, std::vector<Int >& nsrc, QBOOL dummy = QTRUE,
    QBOOL continus = QFALSE, const Int dir = 1)
  {
    TIMER_FLOPS("vector_cs copy from B");
    if(dir == 1){
      Qassert(src.initialized);
      if(!initialized){resize(src.nvec, src.nsum, src.GPU, src.b_size, src.bfac_group);}
    }
    if(dir == 0){
      Qassert(initialized);
      if(!src.initialized){src.resize(nvec, nsum, GPU, b_size, bfac_group);}
    }
    /////initialize nsrc, nres if needed
    Qassert(nsrc.size() == nres.size());
    if(nsrc.size() == 0 ){
      nsrc.resize(nvec);nres.resize(nvec);
      for(Int ni=0;ni<nvec;ni++){
        nsrc[ni] = ni;nres[ni] = ni;
      }
    }
    for(unsigned int ni=0;ni<nsrc.size();ni++){Qassert(nsrc[ni] < src.nvec);};
    for(unsigned int ni=0;ni<nres.size();ni++){Qassert(nres[ni] <     nvec);};
    if(src.nsum == nsum and src.b_size == b_size and src.bfac_group ==  bfac_group){
      (void)nsum;
    }else{abort_r("data format not matched!");}

    Long& b_size = this->b_size;
    ////Long& btotal = this->btotal;
    Long  Ndata = btotal * b_size;
    /////const bool same_locate = (GPU == src.GPU);
    const bool same_locate = check_GPU_same(GPU, src.GPU);
    Int GPU_set = check_GPU_multi(  GPU, src.GPU);
    QMEM GPU_r = this->GPU;
    QMEM GPU_s = src.GPU;
    size_t roff = 0;
    size_t soff = 0; 

    if(dir == 1)
    {
      GPU_r = this->GPU;
      roff = size_t(b_size) *     nvec ;
      GPU_s = src.GPU;
      soff = size_t(b_size) * src.nvec ;
    }
    if(dir == 0)
    {
      GPU_r = src.GPU;
      roff = size_t(b_size) * src.nvec ;
      GPU_s = this->GPU;
      soff = size_t(b_size) *     nvec ;
    }
    /////qmessage("===roff %8d, n %8d, soff %8d, n %8d \n", int(roff), int(nvec), int(soff), int(src.nvec));

    for(unsigned int n0=0;n0<nsrc.size(); n0++){
      ////qlat::vector_gpu<Ty* > va = src.get_pointers(nsrc[n0]);
      ////qlat::vector_gpu<Ty* > vb =     get_pointers(nres[n0]);
      ////Ty** A = va.data();Ty** B = vb.data();
      Ty** A =     get_pointers(nres[n0]);
      Ta** B = src.get_pointers(nsrc[n0]);
      if(same_locate and continus == QFALSE)
      {
        //qGPU_for2dNB(jd, b_size, id, btotal, GPU_set, {
        //  if(dir == 1)A[id][jd] = B[id][jd];
        //  if(dir == 0)B[id][jd] = A[id][jd];
        //});
        qGPU_forNB(isp, Ndata, GPU_set, {
          Long id = isp / b_size;
          Long jd = isp % b_size; 
          ////const Long isp = id * b_size + jd;
          if(dir == 1)A[id][jd] = B[id][jd];
          if(dir == 0)B[id][jd] = A[id][jd];
        });

      }

      if(!same_locate and continus == QFALSE)
      {
        for(Int id=0;id<btotal/bfac_group;id++)
        {
          ////cpy_GPU(A[id], B[id], nsrc.size()*b_size,  GPU_r, GPU_s, false);
          if(dir == 1)cpy_GPU2D(A[id*bfac_group], B[id*bfac_group], 
            size_t(b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);
          if(dir == 0)cpy_GPU2D(B[id*bfac_group], A[id*bfac_group], 
            size_t(b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);
        }
        //for(Int id=0;id<btotal;id++)
        //{
        //  cpy_GPU(A[id], B[id], b_size,  GPU_r, GPU_s, false);
        //}
      }
    }

    /////additional continus copy to avoid too many GPU calls
    //if(!same_locate and continus == QTRUE)
    //if(continus == QTRUE)
    //if(!same_locate and continus == QTRUE)
    //if(!same_locate and continus == QTRUE)
    if(continus == QTRUE)
    {
      Ty**  A =     get_pointers(nres[0]);
      Ta**  B = src.get_pointers(nsrc[0]);

      //Ty** A = NULL;
      //Ty** B = NULL;
      //if(dir == 1)
      //{
      //Ty**  A =     get_pointers(nres[0]);
      //Ta**  B = src.get_pointers(nsrc[0]);
      //}
      //if(dir == 0)
      //{
      //  A = src.get_pointers(nsrc[0]);
      //  B =     get_pointers(nres[0]);
      //}
      for(Int id=0;id<btotal/bfac_group;id++)
      {
        ////cpy_GPU(A[id], B[id], nsrc.size()*b_size,  GPU_r, GPU_s, false);
        ////if(sizeof(Ty) != sizeof(Ta)){
        ////  qmessage("Check %d %d %d %d!\n", int(sizeof(Ty)), int(sizeof(Ta)), int(GPU_r), int(GPU_s));return ;
        ////}
        if(dir == 1){cpy_GPU2D(A[id*bfac_group], B[id*bfac_group], 
          size_t(nsrc.size()*b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);}
        if(dir == 0){cpy_GPU2D(B[id*bfac_group], A[id*bfac_group], 
          size_t(nsrc.size()*b_size), size_t(bfac_group), roff, soff,  GPU_r, GPU_s, QFALSE);}
      }
    }
    if(dummy==QTRUE){qacc_barrier(dummy);}
    double flops = double(b_size) * nsrc.size() * btotal * sizeof(Ty);
    timer.flops  += flops;
    flops_copy += flops;
  }

  template <class Ta >
  void copy_to(vector_cs<Ta >& res, std::vector<Int >& nsrc, std::vector<Int >& nres, QBOOL dummy = QTRUE, QBOOL continus = QFALSE)
  {
    copy_from(res, nsrc, nres, dummy, continus, 0);
  }

  template <class Ta >
  void copy_from_group_same(vector_cs<Ta >& src, std::vector<Int >& nA, std::vector<Int >& nB, QBOOL dummy = QTRUE, Int dir = 1)
  {
    //if(a1-a0 == 0 and b1 - b0 == 0){return ;}
    //Qassert(a1 > a0 and b1 > b0);
    //std::vector<Int > nA;
    //std::vector<Int > nB;
    //const Int na = a1 - a0;
    //const Int nb = b1 - b0;
    //Qassert(na == nb);
    //nA.resize(na);nB.resize(nb);
    //for(Int a=0;a<na;a++){nA[a] = a0 + a;}
    //for(Int b=0;b<nb;b++){nB[b] = b0 + b;}
    const QBOOL continus = QTRUE;
    copy_from(src, nA, nB, dummy, continus, dir);
  }

  template <class Ta >
  void copy_from_group(vector_cs<Ta >& src, Int a0, Int a1, Int b0, Int b1, QBOOL dummy = QTRUE, Int dir = 1)
  {
    if(a1-a0 == 0 and b1 - b0 == 0){return ;}
    Qassert(a1 > a0 and b1 > b0);
    std::vector<Int > nA;
    std::vector<Int > nB;
    const Int na = a1 - a0;
    const Int nb = b1 - b0;
    Qassert(na == nb);
    nA.resize(na);nB.resize(nb);
    for(Int a=0;a<na;a++){nA[a] = a0 + a;}
    for(Int b=0;b<nb;b++){nB[b] = b0 + b;}

    Int do_same = 0;
    Qassert(initialized or src.initialized);
    if( initialized and !src.initialized){do_same = 1;}
    if(!initialized and  src.initialized){do_same = 1;}

    if(do_same == 0)
    if(src.nsum == nsum and src.b_size == b_size and src.bfac_group ==  bfac_group){
      do_same = 1;
    }
    if(do_same == 1){
      copy_from_group_same(src, nA, nB, dummy, dir);
    }else{
      TIMERA("vector_cs copy from group diff");
      Qassert(initialized and src.initialized);
      Int GPU_set = check_GPU_multi(  GPU, src.GPU);
      if(GPU_set == -2){GPU_set = -1;}
      //Qassert(GPU_set != -2);
      VectorGPUKey gkey(0, ssprintf("vector_cs_buf"), GPU_set);
      vector_gpu<int8_t >& tmp = get_vector_gpu_plan<int8_t >(gkey);tmp.resizeL(sizeof(Ty) * nsum);
      Ty* buf = (Ty*) tmp.data();
      if(dir == 1)
      for(unsigned int vi=0;vi<nA.size();vi++)
      {
        src.copy_to(buf, nB[vi], GPU_set, QTRUE);
        copy_from( buf , nA[vi], GPU_set, QTRUE);
      }

      if(dir == 0)
      for(unsigned int vi=0;vi<nA.size();vi++)
      {
        copy_to( buf , nA[vi], GPU_set, QTRUE);
        src.copy_from(buf, nB[vi], GPU_set, QTRUE);
      }
    }

    //if(dir == 1){
    //  Qassert(src.initialized);
    //  if(!initialized){resize(src.nvec, src.nsum, src.GPU, src.b_size, src.bfac_group);}
    //}
    //if(dir == 0){
    //  Qassert(initialized);
    //  if(!src.initialized){src.resize(nvec, nsum, GPU, b_size, bfac_group);}
    //}
    ///////initialize nsrc, nres if needed
    //Qassert(nsrc.size() == nres.size());
    //if(nsrc.size() == 0 ){
    //  nsrc.resize(nvec);nres.resize(nvec);
    //  for(Int ni=0;ni<nvec;ni++){
    //    nsrc[ni] = ni;nres[ni] = ni;
    //  }
    //}

    //const QBOOL continus = QTRUE;
    //copy_from(src, nA, nB, dummy, continus, dir);
  }

  ////===added, -1 for copy all vecs, 
  template <class Ta >
  inline void copy_from(vector_cs<Ta >& src,Int ir = -1, Int is = -1, QBOOL dummy = QTRUE, Int dir = 1)
  {
    if(src.nvec == 0 and ir == -1 and is == -1){return ;}
    if(dir == 1){
      Qassert(src.initialized); 
      if(!initialized){
        if(ir == -1){resize(src.nvec, src);}
        if(ir != -1){Qassert(ir == 0);resize(1, src);}////default resize only to 1
      }
    }
    if(dir == 0){
      Qassert(    initialized); 
      if(!src.initialized)
      {
        if(is == -1){src.resize(nvec, *this);}
        if(is != -1){Qassert(is == 0);src.resize(1, *this);}
      }
    }
    Qassert(nvec > ir);Qassert(src.nvec > is);
    std::vector<Int > ar(0);std::vector<Int > as(0);
    ////copy full
    if(ir == -1 and is == -1){
      Qassert(nvec == src.nvec);
      ar.resize(0);as.resize(0);
    }
    else
    {
      ////copy one
      if(ir == -1 and is != -1){Qassert(    nvec == 1); ir = 0;}
      if(ir != -1 and is == -1){Qassert(src.nvec == 1); is = 0;}
      /////if(ir != -1 and is != -1){ir = 0;is = 0;}
      as.resize(1);ar.resize(1);
      as[0] = is;ar[0] = ir;
    }

    const QBOOL continus = QTRUE;
    copy_from(src, ar, as, dummy, continus, dir);
  }

  template <class Ta >
  void copy_to(vector_cs<Ta >& res, Int ia = -1, Int is = -1, QBOOL dummy = QTRUE)
  {
    copy_from(res, ia, is, dummy, 0);
  }

  inline Long v_size(){
    return nvec;
  }

  ////==added
  void swap(vector_cs<Ty >& vp, Int ia = 0, Int ib = 0)
  {
    TIMERA("vector_cs swap");
    Qassert(ia < v_size() and ib < vp.v_size());
    bool GPU_set = true;if(GPU == 0){GPU_set = false;}
    VectorGPUKey gkey(0, ssprintf("vector_cs_buf"), GPU_set);
    vector_gpu<int8_t >& tmp = get_vector_gpu_plan<int8_t >(gkey);tmp.resizeL(sizeof(Ty) * nsum);
    Ty* buf = (Ty*) tmp.data();

    vp.copy_to(buf, ib, GPU_set);
    vp.copy_from(*this, ib, ia);
    copy_from(buf, ia, GPU_set);

    ///vector_cs<Ty > tmp;tmp.resize(1, vp);
    //tmp.copy_from(vp, ib);
  }

  template <typename Tf >
  Ty reduceT(Int ia, Ty** s1 = NULL)
  {
    Qassert(initialized);
    bool GPU_set = true;if(GPU == 0){GPU_set = false;}
    size_t size = btotal * b_size;
    //const Long loop = b_size;
    if(norm2_buf.size() != size){norm2_buf.resize(size, GPU_set);}
    Tf* r  = (Tf*) norm2_buf.data();
    {
      //TIMERA("vector_cs norm2_vec copy");
      //qlat::vector_gpu<Ty* > data = get_pointers(ia);
      //Ty** s = data.data();
      Ty** s = get_pointers(ia);

      const Long& b_size = this->b_size;
      const Long& btotal = this->btotal;
      const Long Ndata = btotal * b_size;
      if(s1 == NULL){
        qGPU_for(isp, Ndata, GPU_set, {
          const Long id = isp / b_size;
          const Long jd = isp % b_size;
          r[id*b_size + jd] = qlat::qnorm(s[id][jd]);
        });
      }
      if(s1 != NULL){
        qGPU_for(isp, Ndata, GPU_set, {
          const Long id = isp / b_size;
          const Long jd = isp % b_size;
          Ty tmp = qlat::qconj(s[id][jd]) * s1[id][jd];
          Tf* r0 = &r[(id*b_size + jd)*2 + 0];
          r0[0] = tmp.real();
          r0[1] = tmp.imag();
        });
      }
    }

    //reduce double / float
    Ty rsum = 0.0;
    if(s1 == NULL){
      rsum = Reduce((Tf*) r,   size, GPU_set);
    }
    if(s1 != NULL){
      rsum = Reduce((Ty*) r,   size, GPU_set);
    }

    ////Ty rsum = 0;Ty rre = 0;
    //qlat::vector<Ty > rsum;rsum.resize(2);rsum[0] = 0.0;

    //{
    //print_numbers((Ty*) &r[0], 10, GPU_set);
    //}

    //{
    ////TIMERA("norm2_vec reduce");
    //if(s1 == NULL)reduce_vecs(r, (Tf*) rsum.data(),   size, 1, GPU_set);
    //if(s1 != NULL)reduce_vecs(r, (Tf*) rsum.data(), 2*size, 1, GPU_set);
    //}

    ////qmessage("check sum %.8e %.8e \n",  rsum[0].real(), rsum[0].imag() );
    //////fflush_MPI();
    //{
    //TIMERA("norm2_vec global sum");
    //if(s1 == NULL)sum_all_size( (Tf*) rsum.data(), 1, true );
    //if(s1 != NULL)sum_all_size( (Ty*) rsum.data(), 1, true );
    //}
    return rsum;
  }

  ////norm2, ===added
  Ty norm2_vec(Int ia=0)
  {
    TIMERA("vector_cs norm2_vec");

    using D = typename IsBasicDataType<Ty>::ElementaryType;
    return reduceT<D >(ia);
  }

  inline void print_norm2(Int ic = -1)
  {
    Int nzero = 0;
    for(Int ia=0;ia<nvec;ia++)
    {
      if(ic != -1){if(ia != ic){continue ;}}
      Ty t = norm2_vec(ia);
      if(t.real() >  0){
        if(qlat::qnorm(t.real() - 1) < 1e-3){
          qmessage("==norm i %5d, v 1.0 + %+.8e %+.8e \n", ia, t.real() - 1, t.imag());
        }else{
          qmessage("==norm i %5d, v %+.8e %+.8e \n", ia, t.real(), t.imag());
        }
      }else{nzero += 1; }
    }
    if( nzero >  0)qmessage("zero norms %8d, nvec %8d \n", nzero, int(nvec));
  }

  inline void print_checksum(const bool with_sum = false){
    qmessage("eigen ");
    if(with_sum){
      Ty nsum = 0.0;
      for(unsigned int iv=0;iv<buf_g.size();iv++){
        const Ty n = buf_g[iv].norm2();
        nsum += n;
      }
      qmessage("sum %.8e ", nsum.real());
    }

    std::vector<Ty > buf;buf.resize(nsum);
    crc32_t sum = 0;
    for(Int iv=0;iv<nvec;iv++){
      copy_to(buf.data(), iv, false);
      qacc_barrier(dummy);
      sum ^= quick_checksum(buf.data(), buf.size());
    }
    qmessage("checksum %12X . \n", sum);
  }

  inline void print_prod(double cut = 1e-19, Int maxN = -1)
  {
    if(nvec <= 0){return ;}
    if(maxN == -1 or maxN > nvec){maxN = nvec;}
    qlat::vector_cs<Ty >& a = *this;
    qlat::vector<Ty > alpha;alpha.resize(maxN*maxN);
    a.vec_multi(a, alpha.data(), true, 0, maxN, 0, maxN);
    for(Int i=0;i< maxN;i++)
    {
      bool printK = false;
      for(Int j=0;j<maxN;j++)
      {
        if(std::sqrt(qnorm(alpha[i*maxN + j])) > cut and i != j){printK = true;break;}
      }

      if(printK){
        qmessage("vec i %5d norm %+.1e, ", i, alpha[i*maxN + i].real() - 1.0);
        for(Int j=0;j<maxN;j++)
        {
          if(std::sqrt(qnorm(alpha[i*maxN + j])) > cut and i != j){
            qmessage("%d %+.1e %+.1e, ", j, alpha[i*maxN + j].real(), alpha[i*maxN + j].imag());
          }
        }
        qmessage(" \n");
      }
    }
  }

  template <typename Ta, Int civ >
  void copy_from_FieldM(qlat::FieldM<Ta , civ>& src, Int ncur, Int data_GPU = 0, QBOOL dummy = QTRUE, Int dir = 1 )
  {
    Qassert(src.initialized and initialized);
    const fft_desc_basic& fd = get_fft_desc_basic_plan(src.geo());
    Qassert(fd.noden*civ == nsum);
    Ta* s0 = (Ta*) qlat::get_data(src).data();
    copy_from(s0, ncur, data_GPU, dummy, dir);
  }

  template <typename Ta, Int civ >
  void copy_to_FieldM(qlat::FieldM<Ta , civ>& src, Int ncur, Int data_GPU = 0, QBOOL dummy = QTRUE)
  {
    copy_from_FieldM(src, ncur, data_GPU, dummy, 0);
  }

  template <class T >
  vector_cs<Ty>& operator=(const vector_cs<T >& vp)
  {
    (void)vp;
    qmessage("NO SUPPORT yet!\n");
    Qassert(false);
    return *this;
  }

  ////===added
  template <class T >
  void operator_vec(const T alpha_, Int ia = -1)
  {
    TIMERA("vector_cs operator_vec");
    Qassert(initialized);
    if(ia != -1){Qassert(ia < nvec);}
    const Long& b_size = this->b_size;
    ////const Long Ndata = btotal * b_size;
    const Long& btotal = this->btotal;
    const T alpha = alpha_;
    const QMEM GPU_ = GPU;
    //std::vector<qlat::vector_gpu<Ty* > > vL;vL.resize(nvec);
    for(Long ni = 0; ni < nvec; ni++){
      if(ia != -1 and ni != ia){ continue; }
      //qlat::vector_gpu<Ty* > va =     get_pointers(ni);Ty** A = va.data();
      Ty** A = get_pointers(ni);
      const Long Ndata = btotal * b_size;

      qGPU_forNB(isp, Ndata, GPU_, {
        const Long id = isp / b_size;
        const Long jd = isp % b_size;
        A[id][jd] *= alpha;
      });
    }
    qacc_barrier(dummy);
  }


  ////===added
  template <typename Ta, typename T >
  void operator_vec_to(Ta* b, const T alpha_, Int ia=0, bool GPU_src = true, Int dir = 1)
  {
    TIMERA("vector_cs operator_vec");
    Qassert(ia < nvec);
    //const Long& loop = b_size;
    const Long& b_size = this->b_size;
    const Long& btotal = this->btotal;
    const T alpha = alpha_;
    const Int GPU_ = GPU;
    Qassert(GPU_ == int(GPU_src));
    Qassert(ia >=0);
    Ty** A =   get_pointers(ia);
    const Long Ndata = btotal * b_size;
    qGPU_forNB(isp, Ndata, GPU_, {
      const Long id = isp / b_size;
      const Long jd = isp % b_size;
      if(dir == 1){b[id*b_size + jd] += alpha * A[id][jd];}
      if(dir == 0){A[id][jd] += alpha * b[id*b_size + jd];}
    });
    qacc_barrier(dummy);
  }


  ////===added
  template <typename T >
  void operator_vec(vector_cs<Ty>& b, const T alpha_, Int ia=0, Int ib=0)
  {
    TIMERA("vector_cs operator_vec");
    Qassert(initialized and b.initialized);
    Qassert(b.GPU ==   GPU and b.nsum ==   nsum and b.b_size ==   b_size and b.bfac_group ==  bfac_group);
    Qassert(ia < nvec and ib < b.nvec);
    //const Long& loop = b_size;
    const Long& b_size = this->b_size;
    const Long& btotal = this->btotal;
    const T alpha = alpha_;
    const Int GPU_ = GPU;
    Qassert(ia >=0 and ib >=0 );
    Ty** A =   get_pointers(ia);
    Ty** B = b.get_pointers(ib);
    const Long Ndata = btotal * b_size;
    qGPU_forNB(isp, Ndata, GPU_, {
      ////A[i][j] *= alpha;
      const Long id = isp / b_size;
      const Long jd = isp % b_size;
      A[id][jd] += alpha * B[id][jd];
    });

    qacc_barrier(dummy);

  }

  /////*this^\dagger * b, ===added
  inline Ty dot_vec(vector_cs<Ty>& b, Int ia=0, Int ib=0)
  {
    TIMERA("vector_cs dot_vec");
    Qassert(initialized and b.initialized);
    Qassert(b.GPU ==   GPU and b.nsum ==   nsum and b.b_size ==   b_size);

    Ty** s1 = b.get_pointers(ib);
    using D = typename IsBasicDataType<Ty>::ElementaryType;
    return reduceT<D >(ia, s1);
    //return 0.0;
  }

  /////alpha_ij = a_i^* x b_j
  inline qlat::vector_gpu<Ty*> get_alpha_pointers(Ty* aP, Long offA)
  {
    QMEM& GPU = this->GPU;
    ///Long& bfac_group = this->bfac_group;
    //const Long bcut = bfac_group;
    qlat::vector_gpu<Ty*> aV;aV.resize(btotal, GPU);
    qGPU_for(isp, btotal, GPU, { aV[isp] = &aP[(isp) * offA]; });
    return aV;
  }

  /*
    sum the flops of vec_multi, alpha ai --> continus in bi
    nA = aend - amin
    nB = bend - Bmin
    nA * nB global size 
  */
  template<typename Ta>
  inline void vec_multi(vector_cs<Ty >& b, Ta* alpha, bool Conj = true, 
    Int aini = 0, Int aend = -1, Int bini = 0, Int bend = -1, bool do_sum = true)
  {
    TIMER_FLOPS("==vec_multi");
    //Long nvec;
    //int  cs;
    //int  GPU;
    vector_cs<Ty >& a = *this;
    Qassert(check_GPU_same(a.GPU, b.GPU));
    //if(! (b.nsum == a.nsum and b.b_size == a.b_size and b.bfac_group == a.bfac_group))
    //{
    //qmessage("%d %d %d %d %d %d  \n", int(b.nsum), int(a.nsum) , int(b.b_size), int(a.b_size), int(b.bfac_group), int(a.bfac_group));
    //abort_r();}
    Qassert(b.nsum == a.nsum and b.b_size == a.b_size and b.bfac_group == a.bfac_group);
  
    //int GPU_multi = 0;///CPU multiplication
    //if(a.GPU == QMSYNC or a.GPU == QMGPU){GPU_multi = 1;}////GPU multiplication
    Int GPU_multi = check_GPU_multi(a.GPU, b.GPU);

    if(aend == -1){aend = a.nvec;}
    if(bend == -1){bend = b.nvec;}
    Qassert(aend > aini and bend > bini and aend <= a.nvec and bend <= b.nvec);
    const Long nA   = aend - aini;
    const Long nB   = bend - bini;
  
    const bool& GPU  = a.GPU;
    const Long& bfac = a.bfac;
    const Int& b_size= a.b_size;
    ///const Int bfac_group = a.bfac_group;
    ////const Int Nt    = a.Nt;
    //const Int cs = a.cs;

    Long Nres = bfac * nA * (nB);
    //Ty* pa =   get_pointer(na, ba);
    //Ty* pb = b.get_pointer(nb, ba);

    alpha_buf.resizeL(Nres , GPU);
    alpha_buf.set_zero(QTRUE);
    //VectorGPUKey gkey(0, ssprintf("vector_cs_alpha_buf"), GPU);
    //vector_gpu<int8_t >& alpha_V = get_vector_gpu_plan<int8_t >(gkey);
    //alpha_V.resizeL(size_t(Nres)* sizeof(Ty));
    //Ty* alpha_bufP = (Ty*) alpha_V.data();
    ////clean the results
    //zero_Ty(alpha_bufP, Nres, GPU);
    Ty* alpha_bufP = (Ty*) alpha_buf.data();
  
    bool trans = false;///A(m, w) and B(n, w)
    {
    TIMERA("vec_multi Matrix");
    Long m =   nA;
    Long n =   nB;
    Long w = b_size;
  
    qlat::vector_gpu<Ty*> apV = a.get_alpha_pointers(alpha_bufP, nA*nB);
    for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
    {
      Long bji = jobA[jobi*2 + 0]; Long bcut = jobA[jobi*2+1];
      Ty** rpV = &apV[bji];
      Ty** EpV = a.get_pointers(aini, bji);
      Ty** spV = b.get_pointers(bini, bji);
      matrix_prodP(EpV, spV, rpV, m, n, w , bcut, Conj, trans, GPU_multi, QFALSE);
    }
    qacc_barrier(dummy);
    }
  
    if(do_sum)
    {
    const Long Lalpha = nA * nB ;
    zero_Ty(alpha, Lalpha, GPU);

    {
    //TIMERA("vec_multi Reduce alpha");


    Ta* alphaP     = alpha; //(Ty*) alpha.data();
    qGPU_for(xi, nA*nB, GPU,{
      Ty buf = 0.0;
      for(Long bi=0;bi<bfac;bi++){
        buf += alpha_bufP[(bi)*nA*nB + xi];
      }
      alphaP[xi] += buf;
    });
    }
  
    {
    TIMERA("vec_multi alpha Global sum");
    sum_all_size(alpha, Lalpha, GPU);
    }

    }

    long long vGb = btotal * nA * nB * b_size;
    Int Fcount0   = 6 + 2;
    double flops  = vGb*Fcount0;
    timer.flops  += flops;
    flops_matrix += flops;
  }

  template<typename Ta>
  inline void vec_multi(vector_cs<Ty >& b, qlat::vector_gpu<Ta >& alpha, bool Conj = true, 
    Int aini = 0, Int aend = -1, Int bini = 0, Int bend = -1)
  {
    vector_cs<Ty >& a = *this;
    if(aend == -1){aend = a.nvec;}
    if(bend == -1){bend = b.nvec;}
    Qassert(aend > aini and bend > bini and aend <= a.nvec and bend <= b.nvec);
    const Long nA   = aend - aini;
    const Long nB   = bend - bini;
    alpha.resize(nA * nB, a.GPU);
    vec_multi(b, alpha.data(), Conj, aini, aend, bini, bend);
  }

  /*
    b_i += \sum_j alpha_ij * a_j, alpha i--> j continus
    alpha size nA * nB
  */
  template<typename Ta>
  inline void vec_sums(vector_cs<Ty >& b, Ta* alpha, bool Conj = false,
    Int aini = 0, Int aend = -1, Int bini = 0, Int bend = -1)
  {
    TIMER_FLOPS("==vec_sums");
    ////dup only for Nt, chi=2 of a within nvec
    vector_cs<Ty >& a = *this;

    Qassert(check_GPU_same(a.GPU, b.GPU));
    Qassert(b.nsum == a.nsum and b.bfac == a.bfac and b.bfac_group == a.bfac_group);

    if(aend == -1){aend = a.nvec;}
    if(bend == -1){bend = b.nvec;}
    Qassert(aend > aini and bend > bini and aend <= a.nvec and bend <= b.nvec);
    const Long nA   = aend - aini;
    const Long nB   = bend - bini;
    const Long bfac_group = a.bfac_group;
  
    //int GPU_multi = 0;///CPU multiplication
    //if(a.GPU == QMSYNC or a.GPU == QMGPU){GPU_multi = 1;}////GPU multiplication
    Int GPU_multi = check_GPU_multi(a.GPU, b.GPU);
  
    bool trans = true;///A(m, w) and B(w, n)
    const Long m = nB;
    const Long n = b.b_size;
    const Long w = nA;

    alphaG.resizeL(m*w, GPU_multi);
    Ty* alpha_copy = alphaG.data();
    qGPU_for(isp, m*w, GPU_multi, {alpha_copy[isp] = alpha[isp];});

    qlat::vector_gpu<Ty*> apV;apV.resize(bfac_group, GPU_multi);Ty** aP = apV.data();
    qGPU_for(isp, bfac_group, GPU_multi, {aP[isp] = alpha_copy;});
    for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
    {
      Long bji = jobA[jobi*2 + 0]; Long bcut = jobA[jobi*2+1];
      Ty** ApV = a.get_pointers(aini, bji);
      Ty** BpV = b.get_pointers(bini, bji);
      matrix_prodP((Ty**) apV.data(), ApV, BpV, m, n, w , bcut, Conj, trans, GPU_multi, QFALSE);
    }
    qacc_barrier(dummy);

    long long vGb = btotal * nA * nB * b_size;
    Int Fcount0   = 6 + 2;
    double flops = vGb*Fcount0;
    timer.flops  += flops;
    flops_matrix += flops;
  }

  /////////alpha_ij = a_i^* x b_j
  //////alpha_i = v^dagger * b_i,  res = v - \sum alpha_i b_i, --> res, alpha
  //inline void vec_projections(vector_cs<Ty >& b, qlat::vector<Ty >& alpha, Int Nv = -1, bool Conj = true)
  //{
  //  vector_cs<Ty >& a = *this;
  //  ////dup only for Nt, chi=2 of a within nvec
  //  Qassert(b.GPU == a.GPU and b.nsum == a.nsum and b.bfac == a.bfac and b.bfac_group == a.bfac_group);
  //  Qassert(Long(alpha.size()) == (a.nvec * b.nvec));
  //  ////if(T_keep){Qassert(a.bfac % a.Nt == 0);}
  //
  //  Int GPU_multi = 0;///CPU multiplication
  //  if(a.GPU == QMSYNC or a.GPU == QMGPU){GPU_multi = 1;}////GPU multiplication
  //
  //  bool trans = true;///A(m, w) and B(w, n)
  //  const Long m = a.nvec;
  //  const Long n = b.b_size;
  //  const Long w = b.nvec;
  //
  //  ////TODO need improvements
  //  Ty* alp = &alpha[0];
  //  for(Long bi=0;bi<a.bfac;bi++){
  //     Ty* Ap = a.get_pointer(0, bi);
  //     Ty* Bp = b.get_pointer(0, bi);
  //     matrix_prod(alp, Bp, Ap, m, n, w ,1, Conj, trans, GPU_multi, true);
  //  }
  //  qacc_barrier(dummy);
  //}

  template<typename Ta>
  inline bool is_same(vector_cs<Ta >& b)
  {
    vector_cs<Ty >& a = *this;
    bool same = true;
    if(sizeof(Ta) != sizeof(Ty)){same = false;}
    if(!(b.nsum == a.nsum and b.b_size == a.b_size and b.bfac_group == a.bfac_group)){same = false;}
    if(a.nvec == 0 or b.nvec == 0){return false;}
    void* A = (void*) a.get_pointers(0);
    void* B = (void*) b.get_pointers(0);
    if(A != B){same = false;}
    return same;
  }

  ////Nv, end vector of basis, m start vector of basis
  ////return the projection coefficient
  ////remove_last donot remove a1-1 to vec if add_self == 1
  template<typename Ta>
  inline void Projections(qlat::vector<Ta >& alpha, vector_cs<Ty >& vec, Int a0, Int a1, Int b0, Int b1, Int add_self = 0){
    TIMERA("vector_cs projections");
    if(nvec == 0 or a1 - a0 <= 0){alpha.resize(0);return ;}
    Qassert(a1 > a0 and b1 > b0 and a1 <= nvec and b1 <= vec.nvec);
    Int Na = a1 - a0;
    Int Nb = b1 - b0;

    ////if(Nv < 0){Nv = nvec; }
    //Qassert(Nv <= int(nvec) );
    //Qassert(vini <= int(nvec) and vini >= 0 );
    /////need to be self prod for the extra ones
    ///int Nv_multi = a1;
    if(add_self == 1){
      Qassert(Nb == 1);////only when vec is within self
      Na += 1;
      ///Nv_multi += 1;
      Qassert(Na + a0 <= this->nvec and a1 == b0);
      Qassert(this->is_same(vec));
      //Ty** A = this->get_pointers(a1);
      //Ty** B =   vec.get_pointers(b0);
      //Qassert(A == B);
    }

    /////+1 for other buffers
    //alpha_i^* = v_vb^\dagger * v_i
    //v_vb - alpha_i * v_i
    Int GPU_work = GPU;
    const Long Ndata = Na * Nb;
    if(Long(alphaG.size())  < Ndata * 2){alphaG.resizeL(Ndata * 2, GPU);}
    if(Long(alpha.size() )  < Ndata){ alpha.resize( Ndata);}
    zero_Ty(alpha.data(), alpha.size(), true, QTRUE);
    /////for(Int vi=m;vi<Nv;vi++)
    /////{alpha[]}

    Ty* alphaP = alphaG.data();
    Ta* alphaA = alpha.data();

    vec_multi(vec, &alphaP[0], true , a0, a0+Na , b0, b1);
    /////change order of alpha
    if(Na != a1 - a0){Qassert(Nb == 1);}
    qGPU_for(vi, Na*Nb, GPU_work, {
      const Long ai = vi / Nb;
      const Long bi = vi % Nb;

      alphaA[bi*Na + ai]  =      alphaP[vi];////continus within each eigen vectors
      alphaP[Na*Nb + bi*Na + ai] = -1.0*alphaP[vi] ;
    });
    vec_sums( vec, &alphaP[Na*Nb], false, a0, a1, b0, b1);

    //for(Int vi=vini;vi<Nv;vi++){alpha[vi] = -1 * alpha[vi] ;}

    //for(Int vi=vini;vi<Nv;vi++)
    //{
    //  //alpha[vi] = qlat::qconj( vec.dot_vec(*this, vb, vi) );
    //  vec.operator_vec(*this, -1 * alpha[vi], vb, vi );
    //}
  }

  template<typename Ta>
  inline void Projections_all(qlat::vector<Ta >& alpha, vector_cs<Ty >& vec, Int Nm = -1, Int bi = 0){
    if(Nm == -1){Nm = nvec;}
    Projections(alpha, vec, 0, Nm, bi, bi+1);
  }

  inline void projections(vector_cs<Ty >& vec, Int a0, Int a1, Int b0, Int b1, Int add_self = 0){
    qlat::vector<Ty > alpha;
    Projections(alpha, vec, a0, a1, b0, b1, add_self );
  }

  inline void projections_all(vector_cs<Ty >& vec, Int Nm = -1, Int bi = 0){
    qlat::vector<Ty > alpha;
    Projections_all(alpha, vec, Nm, bi);
  }

  ////Nv, end vector of basis, m start vector of basis
  ////return the projection coefficient
  //qlat::vector<Ty > projections(vector_cs<Ty >& vec, Int Nv = -1, Int vb = 0, Int vini = 0, Int add_self = 0){
  //  qlat::vector<Ty > alpha;
  //  TIMERA("vector_cs projections");
  //  if(nvec == 0 or Nv == 0){alpha.resize(0);return alpha;}
  //  if(Nv < 0){Nv = nvec; }
  //  Qassert(Nv <= int(nvec) );
  //  Qassert(vini <= int(nvec) and vini >= 0 );
  //  /////need to be self prod for the extra ones
  //  Int Nv_multi = Nv;
  //  if(add_self == 1){
  //    Nv_multi += 1;
  //    Qassert(Nv_multi <= this->nvec and Nv_multi-1 == vb);
  //    Ty** a = this->get_pointers(Nv_multi-1);
  //    Ty** b =   vec.get_pointers(vb);
  //    Qassert(a == b);
  //  }

  //  /////+1 for other buffers
  //  //alpha_i^* = v_vb^\dagger * v_i
  //  //v_vb - alpha_i * v_i
  //  if(Long(alphaG.size())  < Nv_multi  ){alphaG.resizeL(Nv_multi, GPU);}
  //  if(Long(alpha.size() )  < Nv_multi-vini){ alpha.resize( Nv_multi - vini);}
  //  zero_Ty(alpha.data(), alpha.size(), true, QTRUE);
  //  /////for(Int vi=m;vi<Nv;vi++)
  //  /////{alpha[]}

  //  Ty* alphaP = alphaG.data();
  //  Ty* alphaA = alpha.data();

  //  vec_multi(vec, &alphaP[0], true , vini, Nv_multi, vb, vb+1);
  //  qacc_for(vi, Nv_multi-vini, {
  //    alphaA[vi]  =      alphaP[vi];
  //    alphaP[vi] *= -1 ;
  //  });
  //  vec_sums( vec, &alphaP[0], false, vini, Nv, vb, vb+1);

  //  //for(Int vi=vini;vi<Nv;vi++){alpha[vi] = -1 * alpha[vi] ;}

  //  //for(Int vi=vini;vi<Nv;vi++)
  //  //{
  //  //  //alpha[vi] = qlat::qconj( vec.dot_vec(*this, vb, vi) );
  //  //  vec.operator_vec(*this, -1 * alpha[vi], vb, vi );
  //  //}
  //  return alpha;
  //}
  inline void normalize(Int ai = -1){
    if(nvec == 0){return ;}
    for(Int vi=0;vi<nvec;vi++)
    {
      if(ai != -1 and vi != ai){continue ; }
      Ty norm = Ty(1.0/std::sqrt( norm2_vec(vi).real() ), 0.0);
      operator_vec(norm, vi );
    }
  }

  inline void orthogonalize(Int repeat = 2, Int i0=0, Int i1 = -1){
    TIMERA("orthogonalize");
    if(nvec == 0){return ;}
    Long Nv = nvec;
    if(i1 != -1){Nv = i1;}
    Ty norm = Ty(1.0/std::sqrt( norm2_vec(i0).real() ), 0.0);
    operator_vec(norm,  i0);
    for(Int vi=i0 + 1;vi<Nv;vi++)
    {
      for(Int i=0;i<repeat;i++){
        projections(*this, i0, vi, vi, vi+1);
      }
      Ty norm = Ty(1.0/std::sqrt( norm2_vec(vi).real() ), 0.0);
      operator_vec(norm, vi );
    }
  }

  /* 
    \sum_ni Qts[ni] * this->ni
  */
  inline void linear_combination(vector_cs<Ty >& vec, Ty* Qts, const Int Nsum, Int vb = 0)
  {
    TIMERA("vector_cs linear_combination");
    if(!vec.initialized){vec.resize(vb+1, *this);}
    vec.set_zero(vb);
    vec_sums(vec, Qts, false, 0, Nsum, vb, vb + 1);
    ////const Long Nsize = vec.a[vb].size();
    ////for(Long i=0;i<Nsize;i++){vec.a[vb][i] = 0;}
    //for(Long ni = 0; ni< Nsum;ni++)
    //{
    //  Ty tmp = Qts[ni];
    //  vec.operator_vec(*this, tmp, vb, ni);
    //}
  }

  template<typename Td>
  void rotateQ(const Td* Q,
    const Int c0, const Int c1, const Int N0, const Int Nm, bool transpose = false, Int max_Q_dim = -1)
  {
    TIMER_FLOPS("==rotateQ ");
    Qassert(initialized );
    Int GPU_ = 1;if(GPU == QMCPU){GPU_ = 0;}

    Qassert(c1 <= nvec and Nm <= nvec);Qassert(c0 < c1 and N0 < Nm);
    if(max_Q_dim == -1){max_Q_dim = Nm;}

    const Long Nsize = Nm-N0;
    const Long Csize = c1-c0;
    const bool GPU_Q = false;
    Long Qsize = Csize * Nsize;
    const Long Npass = Qsize;
    if(transpose){Qsize = Qsize * 2;}
    VectorGPUKey gkey(size_t(Qsize)*sizeof(Ty), ssprintf("rotateQ"), GPU_);
    vector_gpu<int8_t >& tmp = get_vector_gpu_plan<int8_t >(gkey);
    ////qlat::vector_gpu<Ty > tmp;tmp.resize(Qsize);
    Ty* Qb = (Ty*) tmp.data();
    {
    TIMERA("Copy Q coefficient");
    if(!transpose)
    for(Int ci=c0;ci<c1;ci++){
      const Td* s  = &Q[ ci * max_Q_dim + 0];
      Ty* r  = &Qb[(ci-c0) * Nsize + 0 ];
      cpy_GPU(r, s, Nsize,  GPU_, GPU_Q, QFALSE);
    }
    if( transpose)
    for(Int ni=N0;ni<Nm;ni++){
      const Td* s  = &Q[ni * max_Q_dim + 0];
      Ty* r  = &Qb[Npass + (ni-N0) * Csize + 0];
      //Ty tmp = 
      cpy_GPU(r, s, Csize,  GPU_, GPU_Q, QFALSE);
    }
    qacc_barrier(dummy);
    /////need to rotate Qb if necessary
    if( transpose)
    {
      qGPU_for(isp, Npass, GPU_, {
        Long ni = isp / Csize;
        Long ci = isp % Csize;
        Qb[ci*Nsize + ni] = Qb[Npass + isp];
      });
    }
    }

    {
    TIMERA("Rotate Matrix prod");
    Qassert(nvec >= Nsize);

    VectorGPUKey gkey(0, ssprintf("vector_cs_buf"), GPU_);
    vector_gpu<int8_t >& buf_V = get_vector_gpu_plan<int8_t >(gkey);
    ///qlat::vector_gpu<Ty > buf_V;
    ////if(buf_V.size() < nvec*b_size){buf_V.resize(nvec*b_size, GPU_);}////resize to maximium to reduce resize
    buf_V.resizeL(size_t(nvec)*b_size * sizeof(Ty));
    Ty* bufP = (Ty*) buf_V.data();

    /////could make use of groups if needed, buf_V may be very large
    for(Int bi=0;bi<btotal;bi++)
    {
      Ty* rs  =   get_pointer_b(N0, bi);
      Ty* rp  =   get_pointer_b(c0, bi);
      Ty* b   =  &bufP[0];
      cpy_GPU(b, rs, Nsize * b_size,  GPU_, GPU_, QTRUE);
      zero_Ty(rp, Csize*b_size, GPU_, QTRUE);

      const bool trans = true;///A(m, w) and B(w, n)
      const bool Conj  = false;
      const Long m = Csize;
      const Long n = b_size;
      const Long w = Nsize;

      Ty* Ap = Qb;
      Ty* Bp = bufP;
      matrix_prod(Ap, Bp, rp, m, n, w ,1, Conj, trans, GPU_, QTRUE);
      qacc_barrier(dummy);
    }

    }
    long long vGb = btotal * Csize * b_size * Nsize;
    Int Fcount0   = 6 + 2;
    timer.flops  += vGb*Fcount0;
  }

};

template <typename Ty >
qacc void set_zero(vector_cs<Ty>& vec)
{
  vec.set_zero();
}

/*
  A = A + B 
  B will be resized
*/
template <typename Ty, typename Tb >
void vector_cs_append(vector_cs<Ty >& A, vector_cs<Tb >& B, Int b0, Int b1, bool clearB = false, Int GPU = 0)
{
  TIMERA("vector_cs_append");
  const Long NB  = b1 - b0;
  const Long NA  = A.nvec;
  const Long Nstop = NA + NB;
  if(B.nvec == 0 or NB <= 0){return ;}
  Qassert(b1 <= B.nvec);
  if(!A.initialized){
    Qassert(GPU == 0 or GPU == 1 or GPU == -1);
    if(GPU == 0){A.resize(1, QMCPU, B);}
    if(GPU == 1){A.resize(1, QMGPU, B);}
    if(GPU ==-1){A.resize(1, QMSYNC, B);}
  }
  Qassert(A.nsum == B.nsum and A.bfac == B.bfac and A.bfac_group == B.bfac_group);

  const Long bfac_group = A.bfac_group;
  const Long btotal     = A.btotal;
  const Long b_size     = A.b_size;
  if(NB == 1){qmessage("=====WARNING! vector_cs_append should be used with group of vectors!" ); }

  //if(btotal / bfac_group <  8){
  //  qmessage("=====WARNING! ROTATION USED A LOT OF MEMORY, REDUCE bfac_group %8d with bfac %8d \n",
  //    bfac_group, btotal);}

  ////int CPU_buf = 0;
  QMEM CPU_buf = B.GPU;
  VectorGPUKey gkey(0, ssprintf("vector_cs_buf"), CPU_buf);
  vector_gpu<int8_t >& buf_V = get_vector_gpu_plan<int8_t >(gkey);
  buf_V.resize(size_t(Nstop)*bfac_group*b_size * sizeof(Ty));
  Ty* bufa = (Ty*) buf_V.data();

  if(A.buf_g.size() == 0){A.buf_g.resize(btotal/bfac_group);}

  for(Int Bi=0;Bi<btotal/bfac_group;Bi++)
  {
    const Int bi = Bi * bfac_group;
    /////for(Int GPUi = 0;GPUi < 2;GPUi++)
    //for(Int bj=0;bj<bfac_group;bj++)
    {
      Ty* pa   = A.get_pointer_b( 0, bi + 0);
      Tb* pb   = B.get_pointer_b(b0, bi + 0);
      if(NA > 0){
        //Ty* pr  =  &bufa[((bj)*(NB + NA) + 0)*b_size ];
        //cpy_GPU(pr, pa, NA * b_size, CPU_buf, A.GPU, QTRUE);
        Ty* pr  =  &bufa[0];
        cpy_GPU2D(pr, pa, NA * b_size, bfac_group, (NB + NA)*b_size, NA*b_size, CPU_buf, A.GPU, QTRUE);
      }
      if(NB > 0){
        //Ty* pr  =  &bufa[((bj)*(NB + NA) + NA)*b_size ];
        //cpy_GPU(pr, pb, NB * b_size,  CPU_buf, B.GPU, QTRUE);
        //qGPU_for(isp, NB * b_size, B.GPU, {pr[isp] = pb[isp];});
        //qGPU_for(isp, NB * b_size, 0, {pr[isp] = 0*pr[isp];});
        Ty* pr  =  &bufa[NA*b_size ];
        cpy_GPU2D(pr, pb, NB * b_size, bfac_group, (NB + NA)*b_size, B.nvec*b_size,  CPU_buf, B.GPU, QTRUE);
      }
      ///if(GPUi == 0 and NA > 0){cpy_GPU(pr, pa, NA * b_size,  CPU_buf, A.GPU, QTRUE);} ///GPU use for memory cpy, 0 CPU, 1 GPU
      ///if(GPUi == 1 and NB > 0){cpy_GPU(pr, pb, NB * b_size,  CPU_buf, B.GPU, QTRUE);} ///GPU use for memory cpy, 0 CPU, 1 GPU
    }

    if(clearB){B.buf_g[Bi].resize(0);}///clear GPU mem
    A.buf_g[Bi].resize(bfac_group * Nstop * b_size, A.GPU);
    ////evecsC.GPU = CPU_buf;evecsC.nvec = Nstop;
    cpy_GPU(A.buf_g[Bi].data(), bufa, bfac_group * Nstop * b_size, A.GPU, CPU_buf);
  }

  {
    const Long nsum   = A.nsum;;
    const Long b_size = A.b_size;;
    const Long bfac_group = A.bfac_group;
    if(clearB){B.clean_mem();}
    A.resize(Nstop, nsum, A.GPU, b_size, bfac_group, true);
  }
  buf_V.resize(0);
}


template <typename Ty>
struct vector_cs_mat{
  qlat::vector<Ty > mat;
  std::vector< qlat::vector<Ty > > bufs;
  qlat::vector<Ty* > mat_src;
  Long vsize;
  Int ncount;
  double mass2_neg;
  double dslash_flops;

  vector_cs_mat(){
    mat.resize(0);
    vsize = 0;
    ncount = 1;
    mass2_neg = 0.0;
    dslash_flops = 0.0;
  }
  Int flops;

  inline Long get_vsize(){
    return vsize;
  }

  template<typename Tc>
  inline void multi(qlat::vector_cs<Tc >& vr, qlat::vector_cs<Tc >& vs, Int ir = 0, Int is = 0){
    vs.copy_to(mat_src[0], is, true);
    multi(mat_src[1], mat_src[0]);
    vr.copy_from(mat_src[1], ir, true);
  }

  //inline void multi(qlat::vector_cs<Tc >& vr, qlat::vector_cs<Tc >& vs, Int ir = 0, Int is = 0){
  //  if(vr.v_size() == 0){vr.copy_from(vs);}
  //  const Long btotal = vr.btotal;
  //  const Long b_size = vr.b_size;
  //  const Long Nd = btotal * b_size;
  //  Qassert(Nd*Nd == Long(mat.size()) );
  //  vr.set_zero(ir);
  //  qlat::vector<Ty > buf;buf.resize(b_size);
  //  qlat::vector<Ty > res;res.resize(b_size);

  //  for(Long bi=0;bi< btotal;bi++){
  //    qlat::set_zero(res);
  //    Ty* sr = vr.get_pointer_b(ir, bi);
  //    for(Long mi=0;mi<b_size;mi++)
  //    {
  //      Long i = bi*b_size + mi;
  //      for(Long bj=0;bj< btotal;bj++){
  //        Ty* ss = vs.get_pointer_b(is, bj);

  //        Ty* matP = &mat[i*Nd + bj*b_size];
  //        ///qlat::set_zero(buf);
  //        qacc_for(mj, b_size, {
  //          buf[mj] = matP[mj] * ss[mj];
  //        } );
  //        for(Long mj=0;mj<b_size;mj++)
  //        {
  //          res[mi] += buf[mj];
  //        }
  //      }
  //    }
  //    qacc_for(mi, b_size, {
  //      sr[mi] += res[mi];
  //    } );
  //  }
  //}

  inline void init_Ndata(const Long Ndata){
    bufs.resize(6);
    mat_src.resize(6);
    vsize = Ndata;
    for(unsigned int i=0;i<bufs.size();i++){
      bufs[i].resize(Ndata);
      mat_src[i] = bufs[i].data();
    }
  }

  //inline void get_src(qlat::vector<Ty* >& mat_src_)
  //{
  //  mat_src_.resize(5);
  //  if(bufs.size() == 0){init_Ndata(std::sqrt( mat.size() ));}
  //  ////if(s0.size() == 0){s0.resize(mat.size());}
  //  ////if(r0.size() == 0){r0.resize(mat.size());}
  //  mat_src_[0] = bufs[0].data();
  //  mat_src_[1] = bufs[1].data();
  //  mat_src_[2] = bufs[2].data();
  //  mat_src_[3] = bufs[3].data();
  //  mat_src_[4] = bufs[4].data();
  //}

  //inline void copy_from_src(qlat::vector_cs<Ty >& src, Int is)
  //{
  //  Qassert(src.initialized);
  //  if(bufs.size() == 0){init_Ndata(std::sqrt( mat.size() ));}
  //  const Long btotal = src.btotal;
  //  const Long b_size = src.b_size;
  //  Ty* r = bufs[0].data();
  //  for(Long bi=0;bi< btotal;bi++){
  //    Ty* ss = src.get_pointer(is, bi);
  //    qacc_for(mi, b_size, {
  //      r[bi*b_size + mi] = ss[mi];
  //    });
  //  }
  //}

  //inline void copy_to_res(Ty* dst, qlat::vector_cs<Ty >& res, Int ir)
  //{
  //  Qassert(res.initialized);
  //  const Long btotal = res.btotal;
  //  const Long b_size = res.b_size;
  //  Ty* s = dst;
  //  for(Long bi=0;bi< btotal;bi++){
  //    Ty* rr = res.get_pointer(ir, bi);
  //    qacc_for(mi, b_size, {
  //      rr[mi] = s[bi*b_size + mi];
  //    });
  //  }
  //}

  inline void multi(Ty* vr,  Ty* vs){
    const Long loop = vsize;////std::sqrt( mat.size() );
    zero_Ty(vr, loop, true );
    if(ncount == 1)
    for(Long i=0;i<loop;i++)
    for(Long j=0;j<loop;j++)
    {
      if(i!=j){vr[i] += mat[i*loop + j] * vs[j];}
      if(i==j){vr[i] += (mat[i*loop + j] + mass2_neg) * vs[j];}
    }

    if(ncount == 2)
    {
      Ty* buf = mat_src[5];
      zero_Ty(buf, loop, true );
      for(Long i=0;i<loop;i++)
      for(Long j=0;j<loop;j++)
      {
        if(i!=j){buf[i] += mat[i*loop + j] * vs[j];}
        if(i==j){buf[i] += (mat[i*loop + j] + mass2_neg) * vs[j];}
      }
      for(Long i=0;i<loop;i++)
      for(Long j=0;j<loop;j++)
      {
        if(i!=j){vr[i] += mat[i*loop + j] * buf[j];}
        if(i==j){vr[i] += (mat[i*loop + j] + mass2_neg) * buf[j];}
      }
    }

    ////estimation of flops, wrong
    dslash_flops += loop * loop * 3.0;

  }



};


}


#endif
