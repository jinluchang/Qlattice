// general_funs.h
// Gen Wang
// Jan. 2021

#ifndef GENERAL_FUNS_H
#define GENERAL_FUNS_H
#pragma once


#include <string.h>
#include <sys/resource.h>
#include <mpi.h>
#include <time.h>
#include <typeinfo>

#include "utils_float_type.h"
#include <type_traits>

#include <iterator>
#include "utils_read_txt.h"
#include "utils_vector_GPU.h"
#include "utils_mpi.h"

namespace qlat
{

////Only cpu verstion
////flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
inline void reorder_civ(int8_t* src, int8_t* res,Int biva,Int civ,size_t sizeF,Int flag,Int size_inner)
{
  //TIMER("reorder_civ vectors int8_t");
  if(biva == 0 or civ == 0 or sizeF == 0 or size_inner == 0){return ;}

  std::vector<int8_t* > psrc;psrc.resize(civ);
  std::vector<int8_t > tmp;tmp.resize(biva*sizeF*civ*size_inner);
  if(size_inner <= 1){abort_r("size_innter too small ! \n");return;}

  if(flag == 1){memcpy((int8_t*)&tmp[0],(int8_t*)&src[0],sizeof(int8_t)*biva*sizeF*civ*size_inner);}
 
  for(size_t bi=0;bi<size_t(biva);bi++)
  {
    for(Int ci=0;ci<civ;ci++)
    {
      if(flag==0)psrc[ci] = &src[(bi*sizeF*civ + ci*sizeF + 0)*size_inner];
      if(flag==1)psrc[ci] = &res[(bi*sizeF*civ + ci*sizeF + 0)*size_inner];
    }

    #pragma omp parallel for
    for(LInt si=0;si<sizeF;si++)
    for(Int ci=0;ci<civ;ci++)
    {
      if(flag==0){
        memcpy(&tmp[(bi*sizeF*civ + si*civ + ci)*size_inner],&psrc[ci][si*size_inner],sizeof(int8_t)*size_inner);
      }
      if(flag==1){
        memcpy(&psrc[ci][si*size_inner],&tmp[(bi*sizeF*civ + si*civ + ci)*size_inner],sizeof(int8_t)*size_inner);
      }
    }
  }
 
  if(flag == 0){memcpy((int8_t*)&res[0],(int8_t*)&tmp[0],biva*sizeF*civ*size_inner);}
}

template<typename Ty>
inline void print_numbers(Ty* src, Int size, Int GPU)
{
  qlat::vector<Ty > buf;buf.resize(size);
  cpy_GPU(buf.data(), src, size, 0, GPU);
  for(Int i=0;i<size;i++)
  {
    qmessage("i %5d, value %+.8e %+.8e \n", i, buf[i].real(), buf[i].imag());
  }
}

///flag = 1 --> biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
#ifdef QLAT_USE_ACC
template <typename Ty, bool flag, Int Threads, Int Biva>
__global__ void move_index_global(Ty* src, Ty* res, Long sizeF, Int civ, Int inner)
{
  __shared__ Ty buf[Threads*Biva];

  const Int    tid = threadIdx.x;
  const Long s0    = blockIdx.x*blockDim.x;

  Int Total = Threads*civ*inner;
  if(s0 + Threads > sizeF){Total = (sizeF - s0) * civ*inner;}

  const Int nT    = Total / (civ * inner);
  const Int nB    = (Total + Threads-1)/Threads;
  const Int nC    = (Total + Biva*Threads-1)/(Biva*Threads);

  Int ci, si, i0;
  Long z0 = 0;Long off = 0;Long off1 = 0;
  for(Int ni=0;ni < nC; ni++)
  {
    if(z0 >= Total){break;}
    if(flag){
    off = z0 + tid;
    for(Int xi=0;xi<Biva;xi++)
    {
      if(off < Total){buf[xi*Threads + tid] = src[s0*civ*inner + off];off += Threads;}
    }
    __syncthreads();
    }

    off = tid;
    for(Int xi=0;xi<nB;xi++)
    {
      ci = off/(nT*inner);
      si = (off/inner)%nT;
      i0 = off%inner;

      off1 = (si*civ + ci)*inner + i0 - z0;
      if(ci < civ and off1 >= 0)
      if((off1 < Threads*Biva) and (off1 < (Total - z0)) )
      {
        if( flag){res[(ci*sizeF+s0+si)*inner + i0] = buf[off1];}
        if(!flag){buf[off1] = src[(ci*sizeF+s0+si)*inner + i0];}
      }
      off += Threads;
    }
    __syncthreads();

    if(!flag){
    off = z0 + tid;
    for(Int xi=0;xi<Biva;xi++)
    {
      if(off < Total){res[s0*civ*inner + off] = buf[xi*Threads + tid];off += Threads;}
    }
    __syncthreads();
    }

    z0 += Threads*Biva;
  }

}
#endif

////TODO change into Ty*
struct move_index
{
  //bool GPU;

  //qlat::vector_gpu<int8_t > buf;
  ////size_t buf_size;
  //qlat::vector<int8_t* > pciv;

  //move_index(bool GPU_set=false){
  //  #ifndef QLAT_USE_ACC
  //  GPU = false;
  //  #else
  //  GPU = GPU_set;
  //  #endif
  //  buf = NULL;
  //  buf_size = 0;
  //}

  //void set_mem(Int civ, size_t Bsize)
  //{
  //  TIMERA("move_index set_mem");
  //  if(buf_size != Bsize){
  //    free_mem();
  //    if(GPU){gpuMalloc(buf, Bsize, int8_t);}
  //    else{buf = (void *)aligned_alloc_no_acc(Bsize);}
  //    buf_size = Bsize;
  //  }
  //  //////psrc.resize(civ);
  //}

  /////order follow src memory order
  template <typename Ty >
  void move_civ_out(Ty* src,Ty* res,Int biva, Long sizeF,Int civ, Int size_inner, bool GPU = false)
  {
    dojob(src, res, biva, civ, sizeF, 1, size_inner, GPU);
  }

  /////order follow src memory order
  template <typename Ty >
  void move_civ_in(Ty* src,Ty* res,Int biva, Int civ, Long sizeF, Int size_inner, bool GPU = false)
  {
    dojob(src, res, biva, civ, sizeF, 0, size_inner, GPU);
  }

  ////flag == 1 : biva * sizeF * civ * size_inner --> biva * civ * sizeF * size_inner
  template <typename Ty >
  void dojob(Ty* src,Ty* res,Int biva,Int civ,Long sizeF,Int flag, Int size_inner, bool GPU = false)
  {
  if(biva == 0 or civ == 0 or sizeF == 0 or size_inner == 0){return ;}
  /////size_t sizeF = sizeF0;

  ////size_t bufN = biva*civ*size_inner*sizeof(Ty)*sizeF;
  const size_t Off = civ*sizeF*size_inner;
  #if PRINT_TIMER>5
  TIMER_FLOPS("reorder index");
  timer.flops += biva*Off*sizeof(Ty);
  #endif

  ////TIMERB("reorder index");
  if(size_inner < 1){qlat::displayln_info(qlat::ssprintf("size_inner too small %d !\n", size_inner));
    MPI_Barrier(get_comm());fflush(stdout);Qassert(false);
  }

  int8_t* tmp_buf = NULL;

  if(src == res){
    VectorGPUKey gkey(Off*sizeof(Ty) / sizeof(int8_t), std::string("move_index_buf"), GPU);
    const vector_gpu<int8_t >& buf = get_vector_gpu_plan<int8_t >(gkey);
    tmp_buf = buf.p;
    ////buf.resize(Off*sizeof(Ty), GPU);
  }
  //pciv.resize(civ);
  Ty* s0;Ty *s1;
  //#ifdef QLAT_USE_ACC
  //if(GPU)
  if(src == res)if((Off*sizeof(Ty)) % sizeof(qlat::ComplexF) != 0){
    qlat::displayln_info(qlat::ssprintf("size not divided by 16, too small. \n"));Qassert(false);}
  ///#endif

  for(Int bi=0;bi<biva;bi++){
    s0 = &src[bi*Off];
    if(src == res){s1 = (Ty*)tmp_buf;}else{s1 = (Ty*) &res[bi*Off];}
    #ifdef QLAT_USE_ACC
    if(GPU){

      {
      const Int Threads = 32;const Int Biva =  (16*16+sizeof(Ty)-1)/sizeof(Ty);
      Long Nb = (sizeF + Threads -1)/Threads;
      dim3 dimBlock(    Threads,    1, 1);
      dim3 dimGrid(     Nb,    1, 1);
      //qmessage("sizeF %d, civ %d, size_inner %d, Biva %d, Nb %d , char %d \n", int(sizeF), int(civ), int(size_inner), Biva, int(Nb), int(sizeof(int8_t)));
      if(flag==0)move_index_global<Ty, false , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      if(flag==1)move_index_global<Ty, true  , Threads, Biva><<< dimGrid, dimBlock >>>(s0, s1, sizeF, civ, size_inner);
      qacc_barrier(dummy);
      }

      //print_numbers((qlat::ComplexF*) s1, 100, GPU);
      if(src == res){
        const Long Nvol = Long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
        //cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 1);
        cpy_GPU((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 1, 1);
      }

    continue ;}
    #endif

    #pragma omp parallel for
    for(Long   si=0;si<sizeF;si++)
    for(Int    ci=0;ci<civ;ci++)
    {
      Ty* p0=NULL;Ty* p1=NULL;
      if(flag == 1){
        p0 = (Ty*) &s0[(si*civ   + ci)*size_inner];
        p1 = (Ty*) &s1[(ci*sizeF + si)*size_inner];
      }
      if(flag == 0){
        p0 = (Ty*) &s0[(ci*sizeF + si)*size_inner];
        p1 = (Ty*) &s1[(si*civ   + ci)*size_inner];
      }
      memcpy(p1, p0, sizeof(Ty)*size_inner);
    }

    if(src == res){
      Long Nvol = Long(Off*sizeof(Ty)/sizeof(qlat::ComplexF));
      cpy_data_thread((qlat::ComplexF*) &res[bi*Off], (qlat::ComplexF*) s1, Nvol, 0);
    }

  }

  }

  void free_mem(){
    VectorGPUKey gkey0(0, std::string("move_index_buf"), false);
    VectorGPUKey gkey1(0, std::string("move_index_buf"),  true);
    safe_free_vector_gpu_plan<int8_t >(gkey0, false);
    safe_free_vector_gpu_plan<int8_t >(gkey1, false);
    //buf.resize(0);
  }

  ~move_index(){
    free_mem();
  }

};

template <typename Ty>
inline void quick_move_civ_in(Ty* res,
  const Int biva, const Int civ, const Long sizeF, const Long size_inner, const Int dir = 1, bool GPU = true)
{
  VectorGPUKey gkey(0, ssprintf("quick_move_civ_buf"), GPU);
  const long nsum = civ * sizeF * size_inner;
  const long V    =       sizeF * size_inner;
  vector_gpu<int8_t >& tmp = get_vector_gpu_plan<int8_t >(gkey);tmp.resizeL(sizeof(Ty) * nsum);
  Ty* buf = (Ty*) tmp.data();
  for(Int bi = 0; bi < biva; bi++){
    Ty* src = &res[bi * nsum + 0];
    if(dir == 1){cpy_GPU(buf, src, nsum, GPU, GPU);}
    qGPU_for(isp, V, GPU, {
      const long bi = isp / size_inner;
      const long bj = isp % size_inner;
      if(dir == 1)
      for(Int ci = 0; ci < civ; ci++){
        src[(bi*civ + ci)*size_inner + bj] = buf[ci*V + bi*size_inner + bj   ];
      }

      if(dir == 0)
      for(Int ci = 0; ci < civ; ci++){
        buf[ci*V + bi*size_inner + bj    ] = src[(bi*civ + ci)*size_inner + bj];
      }
    });
    if(dir == 0){cpy_GPU(src, buf, nsum, GPU, GPU);}
  }
}

inline void clear_move_index_mem(){
  VectorGPUKey gkey0(0, std::string("move_index_buf"), false);
  VectorGPUKey gkey1(0, std::string("move_index_buf"),  true);
  VectorGPUKey qkey0(0, std::string("quick_move_civ_buf"),  false);
  VectorGPUKey qkey1(0, std::string("quick_move_civ_buf"),  true);
  safe_free_vector_gpu_plan<int8_t >(gkey0, false);
  safe_free_vector_gpu_plan<int8_t >(gkey1, false);
  safe_free_vector_gpu_plan<int8_t >(qkey0, true );
  safe_free_vector_gpu_plan<int8_t >(qkey1, true );
}

//#ifdef QLAT_USE_ACC
//#define MAX_CUDA_STEAM 100
//static qacc_Stream_t Qstream[MAX_CUDA_STEAM];
//#endif

inline void set_GPU(Int set_gpu_id = 1){
  (void)set_gpu_id;
  #ifdef QLAT_USE_ACC
  TIMER("setup GPU");
  Int num_node;MPI_Comm_size(get_comm(), &num_node);
  Int id_node;MPI_Comm_rank(get_comm(), &id_node);

  //cuInit(0);
  Int num_gpus = 0;
  qacc_ErrCheck(qacc_GetDeviceCount(&num_gpus));
  //qmessage("numG %5d \n", num_gpus);fflush(stdout);fflush_MPI();
  Qassert(num_gpus != 0);

  // get local rank
  Int localRank  = -1;
  Int localSize  =  0;
  Int globalRank =  0;
  {
  MPI_Comm_rank(get_comm(), &globalRank);
  //Qassert(globalRank == get_id_node());
  // node local comm
  MPI_Comm nodeComm;
  MPI_Comm_split_type(get_comm(), MPI_COMM_TYPE_SHARED, globalRank,
                      MPI_INFO_NULL, &nodeComm);

  // id within the node
  MPI_Comm_rank(nodeComm, &localRank);
  MPI_Comm_size(nodeComm, &localSize);


  }

  Int qlat_set_device = 1;
  //get environment variable to set gpu device
  {
    std::string val = qlat::get_env(std::string("qlat_set_GPU_device"));
    if(val != ""){qlat_set_device = stringtonum(val);}
  }

  //Qassert(localSize == num_gpus and localRank >= 0);//same number of GPUs
  if(set_gpu_id == 1 and qlat_set_device == 1)
  {
    qacc_Errchk(qacc_SetDevice(localRank % num_gpus));
  }
  Int gpu_id = -1; 
  qacc_ErrCheck(qacc_GetDevice(&gpu_id));

  Int gpu_verbos = 0;
  std::string val = qlat::get_env(std::string("qlat_GPU_verbos"));
  if(val != ""){gpu_verbos = stringtonum(val);}

  if(gpu_verbos){
    Int masterRank = -1;
    Int masterSize =  0;

    // comm across node (each node select one process with the same local rank)
    MPI_Comm masterComm;
    MPI_Comm_split(get_comm(), localRank, globalRank, &masterComm);

    MPI_Comm_rank(masterComm, &masterRank);
    // size of each master comm
    MPI_Comm_size(masterComm, &masterSize);

    char host_name[500];
    Qassert(gethostname(host_name, 500) == 0);
    printf("node info lR %d, lS %d, mR %d, mS %d, gi %d, Ng %d, Ni %3d / %3d , host %s \n",
      localRank, localSize, masterRank, masterSize, gpu_id, num_gpus, id_node, num_node, host_name);
  }

  {
  if(gpu_verbos){
    printf("CPU node %d (of %d) uses CUDA device %d (of %d) \n", id_node, num_node, gpu_id, num_gpus);
  }
  }

  fflush(stdout);
  MPI_Barrier(get_comm());
  #endif

}

inline void set_GPU_threads(Int mode=0){
  //////Set up gpu map to cpu
  (void)mode;
  #ifdef QLAT_USE_ACC
  Int num_gpus = 0;
  qacc_ErrCheck(qacc_GetDeviceCount(&num_gpus));
  qacc_ErrCheck(qacc_DeviceReset());
  if(mode == 0){
  #pragma omp parallel
  {
    unsigned int cpu_thread_id = omp_get_thread_num();
    unsigned int num_cpu_threads = omp_get_num_threads();
    qacc_ErrCheck(qacc_SetDevice(cpu_thread_id % num_gpus));
    Int gpu_id = -1; 
    qacc_ErrCheck(qacc_GetDevice(&gpu_id));
    printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);
  }}

  if(mode == 1){
  #pragma omp parallel
  {
    unsigned int cpu_thread_id = omp_get_thread_num();
    unsigned int num_cpu_threads = omp_get_num_threads();
    if(cpu_thread_id%num_gpus != 0){qmessage("Wrong mapping of omp !\n");Qassert(false);}
    Int Nthreads = cpu_thread_id/num_gpus;

    qacc_ErrCheck(qacc_SetDevice(cpu_thread_id / Nthreads));
    Int gpu_id = -1; 
    printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);
    qacc_ErrCheck(qacc_GetDevice(&gpu_id));
  }}

  #endif

}

inline Int init_mpi_thread(int* argc, char **argv[], Int mode = 3)
{
  Int provided;
  if(mode == 0)MPI_Init_thread(argc, argv, MPI_THREAD_SINGLE, &provided);
  if(mode == 1)MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
  if(mode == 2)MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);
  if(mode == 3)MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

  Int num_node;
  MPI_Comm_size(MPI_COMM_WORLD, &num_node);
  Int id_node;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
  if (0 == id_node) {
    displayln("qlat::begin(): " +
              ssprintf("MPI Initialized. num_node = %d", num_node));
  }

  return num_node;
}

inline void begin_thread(
    Int* argc, char** argv[],
    const std::vector<Coordinate>& size_node_list = std::vector<Coordinate>())
// begin Qlat and initialize a new comm
{
  const Int num_node = init_mpi_thread(argc, argv);
  Coordinate size_node;
  for (Int i = 0; i < (int)size_node_list.size(); ++i) {
    size_node = size_node_list[i];
    if (num_node == product(size_node)) {
      break;
    }   
  }
  if (num_node != product(size_node)) {
    size_node = plan_size_node(num_node);
  }
  begin_comm(get_comm(), size_node);
}

inline double get_mem_GPU_info()
{
  fflush_MPI();
  double freeD = 0;double totalD=0;
  #ifdef QLAT_USE_ACC
  size_t freeM = 0;size_t totalM = 0;
  qacc_ErrCheck(qacc_MemGetInfo(&freeM, &totalM));
  freeD = freeM*pow(0.5,30);
  totalD = totalM*pow(0.5,30);
  #endif
  if(totalD == 0){totalD = 1.0;}
  return freeD / totalD ;
}

inline double get_mem_CPU_info()
{
  fflush_MPI();
  double freeD = 0;double totalD=0;
  #ifndef QLAT_NO_SYSINFO
  struct sysinfo s_info;
  sysinfo(&s_info);
  freeD  = s_info.freeram*pow(0.5,30);
  totalD = s_info.totalram*pow(0.5,30);
  #endif
  if(totalD == 0){totalD = 1.0;}
  return freeD / totalD ;
}

inline void print_mem_info(std::string stmp = "")
{
  fflush_MPI();
  qmessage("%s, ",stmp.c_str());
  #ifdef QLAT_USE_ACC
  double freeD = 0;double totalD=0;
  size_t freeM = 0;size_t totalM = 0;
  qacc_ErrCheck(qacc_MemGetInfo(&freeM,&totalM));
  freeD = freeM*pow(0.5,30);totalD = totalM*pow(0.5,30);
  #endif
  #ifndef QLAT_NO_SYSINFO
  struct sysinfo s_info;
  sysinfo(&s_info);
  #ifdef QLAT_USE_ACC
  qmessage("===CPU free %.8e GB, total %.3e GB; GPU free %.8e GB, total %.3e GB. \n"
          , s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30),freeD,totalD);
  //Qassert(freeD / totalD > 0.08);//unified memeory expand need more memory
  #else
  qmessage("===CPU free %.8e GB, total %.3e GB. \n"
          , s_info.freeram*pow(0.5,30),s_info.totalram*pow(0.5,30));
  #endif
  #else
  qmessage("===MAC free infinity! \n");
  #endif
}


inline Int read_vector(const char *filename, std::vector<double > &dat)
{
  Int prods = 0; 
  unsigned long Vsize = 0; 
  ////{synchronize();fflush(stdout);}

  if(qlat::get_id_node() == 0)
  {
    Vsize = get_file_size_o(filename);
    if(Vsize == 0){prods = 0;}else{prods=Vsize;}
    Vsize = Vsize/8;
  }
  sum_all_size(&prods,1);
  if(prods==0)return prods;

  sum_all_size((int*)&Vsize,1);
  dat.resize(Vsize);

  if(qlat::get_id_node() == 0)
  {
    FILE* filer = fopen(filename, "rb");
    unsigned long count = 1024*1024;
    unsigned long sizec = 0; 
    unsigned long offr  = 0; 
    for(unsigned int iv=0;iv<Vsize;iv++)
    {    
      if(offr >= Vsize*8)break;
      char* buf = (char *)&dat[offr/8];
      if((offr + count) <= (Vsize*8)){sizec = count;}
      else{sizec = Vsize*8 - offr;}

      //fread(buf, 1, sizec, filer);
      long sizec_read = fread(buf, 1, sizec, filer);
      Qassert(sizec_read == (long)sizec);
      offr = offr + sizec;
    }    

    fclose(filer);
  }

  sum_all_size(&dat[0],Vsize);
  return prods;

}

template<typename Yl>
void p_vector(const Yl teml)
{
  std::cout << teml << " ";
}

template<typename Ty>
void p_vector(const std::vector<Ty> teml)
{
  for(unsigned long i=0;i< teml.size();i++)
  {
    p_vector(teml[i]);
  }
  std::cout << std::endl;
}

template<typename Ty>
void p_vector(const qlat::vector<Ty> teml)
{
  for(unsigned long i=0;i< teml.size();i++)
  {
    p_vector(teml[i]);
  }
  std::cout << std::endl;
}

template<typename Ty>
inline void random_Ty(Ty* a, Long N0,Int GPU=0, Int seed = 0, const Int mode = 0)
{
  (void)GPU;
  TIMERA("random_Ty");
  if(N0 == 0)return;
  qlat::RngState rs(qlat::get_id_node() + 1 + seed );

  double ini = qlat::u_rand_gen(rs);
  Long bfac = Long(std::sqrt(N0));
  Long Nuse = Long(N0/bfac + 1);
  qGPU_for(isp, Nuse, GPU, {
    for(Long i=0;i<bfac;i++){
      size_t off = isp*bfac + i;
      if(off < size_t(N0)){
        if(mode==0){a[off] = Ty(std::cos((ini+isp+i)*0.5) , (5.0/(isp+1+i))*ini*0.1);}
        if(mode==1){a[off] = Ty(std::cos((ini+isp+i)*0.5) , std::sin(5.0*(isp+1+i))*ini*0.1);}
      }
    }
  });
}

template<typename Ty>
inline void random_vec(vector<Ty >& buf, Int seed = 0, const Int mode = 0)
{
  Int GPU = 1;
  if(buf.mem_type == MemType::Cpu){GPU = 0;}
  random_Ty((Ty*) buf.data(), buf.size(), GPU, seed, mode);
}

template<typename Ty>
inline double norm_vec(Ty* buf, const size_t Nd, const MemType mem_type = MemType::Acc)
{
  const Int off  = 16;
  const Long  Nf = (Nd + off - 1) / off;
  vector<double > tmp;
  tmp.set_mem_type(mem_type);
  tmp.resize(Nf);
  set_zero(tmp);
  qmem_for(isp, Nf, mem_type, {
    for(Int i=0;i<off;i++)
    {
      const Long idx = isp * off + i;
      if(idx > (Long)Nd){continue;}
      tmp[isp] += qnorm(buf[idx]);
    }
  });
  const double r = reduce_vecs(tmp);
  return r;
}

template<typename Ty>
inline double norm_vec(vector<Ty >& buf)
{
  if(buf.size() == 0){return 0;}
  const Long  Nd = buf.size();
  const double r = norm_vec((Ty*) buf.data(), Nd, buf.mem_type);
  return r;
}

template<typename Ty>
inline void random_EigenM(qlat::vector<Ty >& a,Int GPU=0, Int seed = 0)
{
  Ty* buf = a.data();
  random_Ty(buf, a.size(), GPU, seed);
}

template<typename Ty>
inline void random_EigenM(std::vector<qlat::vector<Ty > >& a, Int GPU=0, Int seed = 0)
{
  Int N0 = a.size();if(N0 == 0)return ;
  for(size_t i=0;i < size_t(N0);i++){random_EigenM(a[i], GPU,  seed + i);}
}

template<typename Ty>
inline void zeroE(qlat::vector<Ty >& a,Int GPU=0, QBOOL dummy=QTRUE)
{
  zero_Ty(a.data(), a.size(), GPU, dummy);
}

template<typename Ty>
inline void zeroE(std::vector<qlat::vector<Ty > >& a,Int GPU=0, QBOOL dummy=QTRUE)
{
  for(LInt iv=0;iv<a.size();iv++){zeroE(a[iv], GPU, false);}
  if(dummy==QTRUE){qacc_barrier(dummy);}
}

template<typename Ty, Int civ>
inline void random_FieldM(qlat::FieldM<Ty , civ>& a,Int GPU=0, Int seed = 0)
{
  Qassert(a.initialized);
  const Geometry& geo = a.geo();
  Ty* buf = (Ty*) qlat::get_data(a).data();
  random_Ty(buf, geo.local_volume() * civ, GPU, seed);
}

template<class Fieldy>
inline void random_FieldG(Fieldy& a,Int GPU=0, Int seed = 0)
{
  Qassert(a.initialized);
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());

  //const Geometry& geo = a.geo();
  const Long Nd = GetFieldSize(a);
  Qassert(Nd % (2 * sizeof(D)) == 0);
  ComplexT<D>* buf = (ComplexT<D>*) qlat::get_data(a).data();
  random_Ty(buf, Nd / (2 * sizeof(D)), GPU, seed);
}

template<typename Ty, Int civ>
inline void copy_FieldM(qlat::FieldM<Ty , civ>& res, qlat::FieldM<Ty , civ>& src )
{
  Qassert(src.initialized);
  const Geometry& geo = src.geo();
  if(!res.initialized){res.init(geo);}
  Ty* a = (Ty*) qlat::get_data(src).data();
  Ty* b = (Ty*) qlat::get_data(res).data();
  const Long  V = geo.local_volume() ;
  cpy_GPU(b, a, V*civ, 1, 1);
}

template<typename Ty>
inline double norm_FieldG(qlat::FieldG<Ty>& a)
{
  Qassert(a.initialized);
  const Geometry& geo = a.geo();
  const Int civ = a.multiplicity;
  //Ty* buf = (Ty*) qlat::get_data(a).data();
  const Long  V = geo.local_volume() ;
  qlat::vector_gpu<Ty > tmp;tmp.resize(V*civ);
  Ty* srcP = (Ty* ) qlat::get_data(a).data();
  Ty* resP = tmp.data();
  qacc_for(isp,V*civ,{ resP[isp] = srcP[isp];});
  Ty  norm = tmp.norm2();
  return norm.real();
  //qmessage("norm %.3e %.3e \n", norm.real(), norm.imag());
}

template<typename Ty, Int civ>
inline Ty norm_FieldM(qlat::FieldM<Ty , civ>& a)
{
  Qassert(a.initialized);
  const Geometry& geo = a.geo();
  //Ty* buf = (Ty*) qlat::get_data(a).data();
  const Long  V = geo.local_volume() ;
  qlat::vector_gpu<Ty > tmp;tmp.resize(V*civ);
  Ty* srcP = (Ty* ) qlat::get_data(a).data();
  Ty* resP = tmp.data();
  qacc_for(isp,V*civ,{ resP[isp] = srcP[isp];});
  Ty  norm = tmp.norm2();
  return norm;
  //qmessage("norm %.3e %.3e \n", norm.real(), norm.imag());
}

template <class Td>
void random_prop(Propagator4dT<Td >& prop, Int seed = -1)
{
  Qassert(prop.initialized);
  ////qmessage("print time %.3f\n", tm.tv_sec);
  Int rand_seed = qlat::get_id_node() + 1;
  if(seed == -1){timeval tm;gettimeofday(&tm, NULL);rand_seed += int(tm.tv_sec);}else{rand_seed += seed;}

  qlat::RngState rs(rand_seed);
  double ini = qlat::u_rand_gen(rs);

  /////int dir_limit = 4;
  const Geometry& geo = prop.geo();

  qacc_for(isp,  geo.local_volume(),{
    qlat::WilsonMatrixT<Td>& v0 =  prop.get_elem_offset(isp);
    for(Int ci=0;ci<12*12;ci++){
      v0.p[ci] = (ci/(12*12.0))* qlat::ComplexT<Td>(std::cos((ini+isp + ci*2)*0.5 + ci) , (ci+(5.0+ci)/(isp+1))*ini*0.1 + 0.2); 
    }
  }); 
}

template <class Td>
void random_link(GaugeFieldT<Td> &gf, const Int seed = -1)
{
  if(seed == -1)
  {
    qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem_offset(isp), 1.0);});
  }else{
    const Geometry& geo = gf.geo();
    qlat::ComplexT<Td>* res = (qlat::ComplexT<Td>*) qlat::get_data(gf).data();
    random_Ty(res, geo.local_volume()*gf.multiplicity*sizeof(ColorMatrixT<Td>)/(sizeof(Td)*2), 1, seed);

    //qacc_for(isp, gf.field.size(), { set_unit(gf.get_elem_offset(isp), 1.0);});
    //ColorMatrixT<Td>& unit = gf.get_elem_offset(0);set_unit(unit, 1.0);
    ColorMatrixT<Td> unit;set_unit(unit, 1.0);
    /////TODO This function cannot be done on GPU
    /////Eigen normalize/normalized problem
    qacc_for(isp, gf.field.size(), {
      gf.get_elem_offset(isp) = gf.get_elem_offset(isp) * (1/2.0) + unit;
      unitarize(gf.get_elem_offset(isp));
    });
  }
}

inline size_t get_threads(size_t thread, size_t Total, Int def=1)
{
  for(size_t temb=thread;temb<Total;temb++)
  {
    if(Total%temb == 0)
    {
      return temb;
    }
  }
  return def;
}

inline std::vector<unsigned int > get_num_power(const size_t x,const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > re; re.resize(a.size());
  for(unsigned int ai=0;ai<a.size();ai++)
  {
    re[ai] = 0;
    unsigned int    fac   = a[ai];
    unsigned int    count = 0;
    size_t use   = x;
    for(unsigned int j = 0;j < 300; j++)
    {
      if(use % fac == 0){
        count += 1;
        use = use/fac;}
      else{
        re[ai] = count;
        break;
      }
    }
  }
  size_t test = 1;
  for(unsigned int ai=0;ai<a.size();ai++){
    test = test * std::pow(a[ai], re[ai]);
  }
  assert(test == x);
  return re;
}

//////Even node in xyzT directions
inline Coordinate spread_even(const Int n, const Coordinate& Lat, const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > Mpow = get_num_power(n, a);
  std::vector<std::vector<unsigned int > > Lpow;
  Lpow.resize(4);for(Int i=0;i<4;i++){Lpow[i] =  get_num_power(Lat[i], a);}
  Coordinate re;
  for(Int i=0;i<4;i++)re[i] = 1;

  std::vector<Int > nL(4);

  for(LInt i=0;i<Mpow.size();i++){
    Int fac = a[i];
    Int num = Mpow[i];
    if(Mpow[i] != 0){
      Int suma = Lpow[0][i] + Lpow[1][i] + Lpow[2][i] + Lpow[3][i] ;
      assert(int(Mpow[i]) <= suma);

      for(Int j=0;j<4;j++){nL[j] = Lpow[j][i];}
      //std::vector<Int>::iterator result = std::max_element(nL.begin(), nL.end());
      //int pos_max = std::distance(nL.begin(), result);
      for(unsigned int ni = 0; ni< Mpow[i] + 1 ;ni++)
      {
        for(Int j=3;j>=0;j--)
        {
          if(nL[j] > 0){num -= 1; nL[j] -= 1; re[j] *= fac;}
          if(num == 0){break;}
        }
        if(num == 0){break;}
      }
    }
  }

  assert(re[0]*re[1]*re[2]*re[3] == n);
  return re;
}


//////Most power in T direction
inline Coordinate spread_powT(const Int n, const Coordinate& Lat, const std::vector<unsigned int >& a)
{
  std::vector<unsigned int > Mpow = get_num_power(n, a);
  std::vector<std::vector<unsigned int > > Lpow;
  Lpow.resize(4);for(Int i=0;i<4;i++){Lpow[i] =  get_num_power(Lat[i], a);}
  Coordinate re;
  for(Int i=0;i<4;i++)re[i] = 1;

  for(LInt i=0;i<Mpow.size();i++){
    Int num = Mpow[i];
    if(num != 0){
      Int suma = Lpow[0][i] + Lpow[1][i] + Lpow[2][i] + Lpow[3][i] ;
      assert(num <= suma);
      unsigned int tem = num;
      for(unsigned int ni=0;ni<4;ni++){
        if(tem >= Lpow[4-ni-1][i]){
          re[4-ni-1] *= int(std::pow(a[i], Lpow[4-ni-1][i]));
          tem = tem - Lpow[4-ni-1][i];
        }
        else{
          re[4-ni-1] *= (unsigned int)(std::pow(a[i], tem));
          tem = 0;
        }
        ////if(tem >= Lpow[4-ni-1][i]){}
        if(tem == 0)break;
      }
    }
  }

  assert(re[0]*re[1]*re[2]*re[3] == n);
  return re;

}

inline Coordinate guess_nodeL(Int n, const Coordinate& Lat, const Int mode = 0)
{

  std::vector<unsigned int > a;a.resize(8);
  a[0] = 2;a[1] = 3;a[2] = 5;a[3] = 7;a[4] =11;a[5] =13;a[6] =17;a[7] =19;
  Coordinate re;
  if(mode == 0){re = spread_powT(n, Lat,a);}
  if(mode == 1){re = spread_even(n, Lat,a);}
  return re;
}

inline void add_nodeL(std::vector<Coordinate>& size_node_list)
{
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 2, 1,  1));
  //size_node_list.push_back(Coordinate(1, 1, 1,  2));
  //size_node_list.push_back(Coordinate(1, 2, 2,  1));
  size_node_list.push_back(Coordinate(1, 2, 1,  2));
  size_node_list.push_back(Coordinate(2, 2, 2,  1));
  size_node_list.push_back(Coordinate(2, 2, 2,  2));
  size_node_list.push_back(Coordinate(1, 1, 3,  1));
  ////size_node_list.push_back(Coordinate(1, 1, 2,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 3,  2));
  size_node_list.push_back(Coordinate(1, 2, 3,  1));
  size_node_list.push_back(Coordinate(1, 2, 4,  1));
  size_node_list.push_back(Coordinate(1, 1, 1, 12));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 6,  4));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(2, 4, 4, 2 ));

  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
  size_node_list.push_back(Coordinate(4, 4, 4,  6));
  size_node_list.push_back(Coordinate(4, 4, 4, 12));
  size_node_list.push_back(Coordinate(4, 4, 8, 12));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));
}

//#include <cuda.h>
//#include <cuda_runtime.h>
//void initialisation_cuda()
//{
//    char* local_rank_env;
//    Int local_rank;
//    qacc_Error_t cudaRet;
// 
//     /* Recovery of the local rank of the process via the environment variable
//        set by Slurm, as  MPI_Comm_rank cannot be used here because this routine
//        is used BEFORE the initialisation of MPI*/
//    local_rank_env = getenv("SLURM_LOCALID");
// 
//    if (local_rank_env) {
//        local_rank = atoi(local_rank_env);
//        /* Define the GPU to use for each MPI process */
//        cudaRet = qacc_ErrCheck(qacc_SetDevice(local_rank));
//        if(cudaRet != CUDA_SUCCESS) {
//            printf("Erreur: qacc_SetDevice has failed\n");
//            exit(1);
//        }
//    } else {
//        printf("Error : impossible to determine the local rank of the process\n");
//        exit(1);
//    }
//}

////mode_dis % 1 == 0, t in diff node, mode_dis % 2 == 1, t in single node
////mode_dis < 2, T, mode_dis >= 2 even
inline void begin_Lat(int* argc, char** argv[], inputpara& in, const Int with_GPU = 1, const Int read_input_Lat = 0){
  //initialisation_cuda();
  if(read_input_Lat >= 0)
  {
    Int n_node = init_mpi(argc, argv);
    in.load_para(*argc, *argv);
    Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
    Coordinate spreadT;
    if(in.layout != std::string("NONE")){
      spreadT = string_to_Coordinate(in.layout);
      if(spreadT[0] * spreadT[1] * spreadT[2] * spreadT[3] != n_node)
      {printf("Wrong input layout\n");abort_r();}
    }
    else{
      if(in.mode_dis%10 >= 0 and in.mode_dis%10 < 2){spreadT = guess_nodeL(n_node, Lat, 0);}
      if(in.mode_dis%10 >= 2 and in.mode_dis%10 < 4){spreadT = guess_nodeL(n_node, Lat, 1);}
      ////if(in.mode_dis >= 10 and in.mode_dis < 12){spreadT = guess_nodeL(n_node, Lat, 1);}
    }
    ///3D begin
    ////begin_comm(MPI_COMM_WORLD , spreadT);

    ///4D begin
    Int id_node, n;
    Int old_id;MPI_Comm_rank(MPI_COMM_WORLD, &old_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    if(in.mode_dis <  10){MPI_Comm_rank(MPI_COMM_WORLD, &id_node);}
    if(in.mode_dis >= 10){
      id_node = get_mpi_id_node_close();
    }
      ///abort_r();
    Int t =  id_node/(spreadT[0]*spreadT[1]*spreadT[2]);
    Int z = (id_node/(spreadT[0]*spreadT[1]))%(spreadT[2]);
    Int y = (id_node/(spreadT[0]))%(spreadT[1]);
    Int x = (id_node%(spreadT[0]));
    ///int new_id = ((z*spreadT[1] + y)*spreadT[0] + x)*spreadT[3] + t;
    Int new_id = ((x*spreadT[1] + y)*spreadT[2] + z)*spreadT[3] + t;
    if(in.mode_dis % 2 == 0)begin(id_node, spreadT);
    if(in.mode_dis % 2 == 1)begin(new_id, spreadT);
    ////printf("new id %5d, old id %5d\n", get_id_node(), old_id);
  }

  if(read_input_Lat == -1)
  {
    std::vector<Coordinate> size_node_list;
    add_nodeL(size_node_list);
    begin(argc, argv, size_node_list);
    in.load_para(*argc, *argv);
  }

  Qassert(with_GPU == 0 or with_GPU == 1);
  if(with_GPU == 1){set_GPU();}

  omp_set_num_threads(omp_get_max_threads());
  qmessage("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  fflush_MPI();
  print_mem_info();

}

inline Int end_Lat(const Int with_mpi = 0, const Int with_timer = 1)
{
  if(with_timer == 1){
    qlat::Timer::display();
  }
  fflush_MPI();
  qlat::end();
  if(with_mpi == 1 and qlat::is_MPI_initialized()){MPI_Finalize();}
  return 0;
}

inline std::vector<Long > job_create(Long total, Long each)
{
  std::vector<Long > a;a.resize(0);
  if(total == 0){return a;}
  if(total < 1 or each < 1){
    qmessage("===Give me valid job types total %ld, each %ld \n", (long)total, (long)each);
    abort_r();}
  /////std::vector<Long > a = job_create(total, each);
  Long jobN  = (total + each - 1)/each;
  Int i0 = 0; Int dj = each;
  for(Int ji = 0; ji < jobN ; ji++)
  {
    if(i0 >= total){break;}
    if(i0 + dj > total){dj = total - i0;}
    a.push_back(i0);
    a.push_back(dj);
    i0 += dj;
  }

  return a;
}

template<typename Ty>
inline void allocate_buf(std::vector<qlat::vector_gpu<Ty > > & buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    buf[i].resize(n1);
  }
}

template<typename Ty>
inline void allocate_buf(std::vector<qlat::vector<Ty > > & buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(0);
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    buf[i].resize(0);
    buf[i].resize(n1);
  }
}

template<typename Ty>
inline void allocate_buf(std::vector<Ty* >& buf, size_t n0, size_t n1)
{
  TIMERA("CUDA Buf mem allocation");
  buf.resize(n0);
  for(LInt i=0;i<buf.size();i++){
    gpuMalloc(buf[i], n1, Ty, 1);
  }
}

template<typename Ty>
inline Ty inv_self(const Ty& lam, double m, double rho,Int one_minus_halfD=1)
{
  ComplexT<double > tem(lam.real(),lam.imag());
  ComplexT<double > v0 = (one_minus_halfD>0)?(1.0-tem/2.0)/(rho*tem+m*(1.0-tem/2.0)):1.0/(rho*tem+m*(1.0-tem/2.0));
  Ty res(v0.real(),v0.imag());
  return res;
}

template<typename Ty>
vector<Ty* > EigenM_to_pointers(std::vector<qlat::vector_gpu<Ty > >& src, Long Nvol = -1)
{
  vector< Ty* >  res;
  const size_t Nvec = src.size();
  if(Nvec == 0){return res;}

  if(Nvol != -1){
    Qassert(src[0].size() % Nvol == 0);
  }else{
    Nvol = src[0].size();
  }

  const size_t Nt = src[0].size() / Nvol;
  res.resize(Nvec * Nt);
  for(size_t iv=0;iv<Nvec;iv++)
  {
    Qassert(src[iv].size() == Nt * Nvol);
    for(size_t it=0;it<Nt;it++)
    {
      res[iv*Nt + it] = &src[iv].data()[it* Nvol];
    }
  }
  return res;
}

/*
  Default each vector into GPU pointers
  Nvol small : devide each src into chucks
*/
template<typename Ty>
vector<Ty* > FieldG_to_pointers(std::vector<FieldG<Ty > >& src, Long Nvol = -1)
{
  qlat::vector<Ty* >  res;
  const size_t Nvec = src.size();
  if(Nvec == 0){return res;}

  Qassert(src.size() > 0 and src[0].initialized);
  const Geometry& geo = src[0].geo();
  const Long Nd = geo.local_volume() * src[0].multiplicity;

  if(Nvol != -1){
    Qassert(Nd % Nvol == 0);
  }else{
    Nvol = Nd;
  }
    
  const size_t Nt = Nd / Nvol;
  res.resize(Nvec * Nt);
  for(size_t iv=0;iv<Nvec;iv++)
  {
    Qassert(src[iv].initialized and src[iv].multiplicity == src[0].multiplicity);
    Ty* sP = (Ty*) get_data(src[iv]).data();
    for(size_t it=0;it<Nt;it++)
    {
      res[iv*Nt + it] = &sP[it* Nvol];
    }
  }
  return res;
}

// Each FieldM to gpu pointers
template<typename Ty>
qlat::vector<Ty* > FieldM_to_pointers(std::vector<qlat::FieldM<Ty, 12*12> >& src)
{
  qlat::vector< Ty* >  res;
  res.resize(src.size());
  for(LInt iv=0;iv<src.size();iv++)
  {
    res[iv] = (Ty*) qlat::get_data(src[iv]).data();
  }
  return res;
}

// Each src ant NT to gpu pointers
template<typename Ty>
qlat::vector<Ty* > FieldM_to_Tpointers(std::vector<qlat::FieldM<Ty, 12*12> >& src, Long Nvol = -1)
{
  qlat::vector< Ty* >  res;
  const size_t Nvec = src.size();
  if(Nvec == 0){return res;}
  for(size_t iv=0;iv<Nvec;iv++){
    Qassert(src[iv].initialized);
  }
  const Geometry& geo = src[0].geo();
  const long V = geo.local_volume();

  if(Nvol != -1){
    Qassert((12*12*V) % Nvol == 0);
  }else{
    Nvol = 12 * 12 * V;
  }

  const size_t Nt = (12 * 12 * V) / Nvol;
  res.resize(Nvec * Nt);
  for(size_t iv=0;iv<Nvec;iv++)
  {
    Ty* p = (Ty*) qlat::get_data(src[iv]).data();
    for(size_t it=0;it<Nt;it++)
    {
      res[iv*Nt + it] = &p[it* Nvol];
    }
  }
  return res;
}

template<typename Ty>
qlat::vector<Ty* > EigenM_to_pointers(std::vector<qlat::vector<Ty > >& src)
{
  qlat::vector<Ty* >  res;
  res.resize(src.size());
  for(LInt iv=0;iv<src.size();iv++)
  {
    res[iv] = src[iv].data();
  }
  return res;
}

/////Ty should have size(), resize(), and data()
template<class Ty>
Ty sum_local_to_global_vector(Ty src, MPI_Comm* commp=NULL)
{
  ////int Nt = geo.node_site[3];
  Int Nmpi  = qlat::get_num_node();
  Int rank  = qlat::get_id_node();
  if(commp != NULL){MPI_Comm_size(*commp, &Nmpi);MPI_Comm_rank(*commp, &rank);}

  qlat::vector<Long > size_global;size_global.resize(Nmpi);
  zero_Ty(size_global.data(), size_global.size(), 0);

  size_global[rank] = src.size();
  sum_all_size(size_global.data(), size_global.size(), 0, commp);

  Long total = 0;Long current = 0;
  for(Int i=0;i<Nmpi;i++){total += size_global[i];if(i < rank){current += size_global[i];}}

  /////for(unsigned int i=0;i<size_global.size();i++){printf("rank %d, size %ld \n", rank, size_global[i]);}
  /////printf("rank %d, Total %ld, current %ld \n", rank, total, current);

  Ty res;res.resize(total);
  for(unsigned long pos=current;pos<src.size();pos++){res[pos] = src[pos - current];}

  sum_all_size(res.data(), res.size(), 0, commp);

  return res;

}

inline std::vector<Int > num_to_site(const Long num, const std::vector<Int > key_T)
{
  Qassert(key_T.size() > 0);
  Int dim = key_T.size();
  std::vector<Int > site;site.resize(dim);
  for(Int iv=0;iv<dim;iv++){site[iv] = 0;}

  Long tem_i = num;
  for(Int Ni=0; Ni < dim; Ni++)
  {
    Long N_T = 1;
    for(Int numi = Ni+1;numi < dim;numi++)
    {
      N_T = N_T*key_T[numi];
    }
    site[Ni] = tem_i/N_T;
    tem_i = tem_i%N_T;
  }
  return site;
}

inline std::vector<Long > random_list(const Long n, const Long m, const Int seed)
{
  std::vector<Long > r;
  std::vector<Long > b;
  qlat::RngState rs(13021 + seed);
  for(Long ri = 0; ri < n; ri++){r.push_back(ri);}
  if(m >= n or m <  0){return r;}
  if(m == 0){r.resize(0);return r;}

  for(Long ri=0;ri<m;ri++){
    Long u = int(qlat::u_rand_gen(rs) * r.size());
    b.push_back(r[u]);
    r.erase(r.begin()+ u);
  }
  return b;

}

template<typename Ty, typename Int>
inline Ty Reduce(Ty* buf, Int Ndata, Int GPU = 1)
{
  TIMERB("Reduce");
  qlat::vector<Ty > rsum;rsum.resize(1);rsum[0] = 0.0;
  reduce_vecs(buf, rsum.data(), Ndata, 1, GPU);
  Ty tmp = rsum[0];
  sum_all_size( &tmp, 1, 0 );
  return tmp;
}

template<typename Ty>
inline Ty vec_norm2(Ty* s0, Ty* s1, Long Ndata, QMEM GPU = QMGPU, const Long Ngroup = 4)
{
  TIMERB("vec_norm2 single");

  VectorGPUKey gkey(0, ssprintf("vec_norm2_buf"), GPU);
  Qassert(Ndata % Ngroup == 0);
  const Long Nvol = Ndata / Ngroup;

  vector_gpu<int8_t >& Buf = get_vector_gpu_plan<int8_t >(gkey);
  Buf.resizeL(size_t(Nvol)* sizeof(Ty));
  Ty* buf = (Ty*) Buf.data();

  qGPU_for(isp, Nvol, GPU, {
    buf[isp] = 0.0;
    for(Int gi=0;gi<Ngroup;gi++){
      const Long index = gi * Nvol + isp;
      buf[isp] += qlat::qconj(s0[index]) * s1[index];
    }
  });

  return Reduce(buf, Nvol, GPU);
}

template<typename Ty>
std::vector<Long > get_sort_index(Ty* src, Long size)
{
  //std::vector<Ty > copy;copy.resize(size);
  //for(Long it=0;it<size;it++){copy[it] = src[it];}
  std::vector<std::pair<Ty, Long > > vec; 
  for(Long it=0;it<size;it++) 
  {    
    vec.push_back( std::make_pair(src[it], it) );
  }    

  std::sort(vec.begin(), vec.end(), [](std::pair<Ty, Long >& a, std::pair<Ty, Long >& b) { 
    return a.first < b.first;
  });  //check order

  std::vector<Long > index;index.resize(size);
  for(Long it=0;it<size;it++) 
  {
    index[it] = vec[it].second;
  }
  return index;
}

template<typename Ty>
void sort_vectors_by_axis(std::vector<std::vector<Ty > >& src, std::vector<std::vector<Ty > >& res, Int axis = 0)
{
  Int Naxis = src.size();
  Qassert(Naxis > axis);
  res.resize(src.size());

  const Long Ndata = src[0].size();

  std::vector<Long > index = get_sort_index(&src[axis][0], Ndata);
  for(Int ai=0;ai<Naxis;ai++)
  {
    Qassert(Long( src[ai].size() ) == Ndata);
    res[ai].resize(src[ai].size());
    for(Long ni=0;ni<Ndata;ni++)
    {
      res[ai][ni] = src[ai][index[ni]];
    }
  }
}

// find position of first elements
template<typename Ty>
Long std_find(const std::vector<Ty >& src, const Ty& elem)
{
  Long pos = -1;
  if(src.size() == 0){return pos;}
  typename std::vector<Ty >::const_iterator it = std::find(src.begin(), src.end(), elem);;
  //it = std::find(src.begin(), src.end(), elem);
  if(it !=  src.end()){
    pos = std::distance(src.begin(), it);;
  }else{
    pos = -1;
  }
  return pos;
}

template<typename Ty>
void std_erase(std::vector<Ty >& src, const Long offset)
{
  Qassert(offset < Long(src.size()));
  src.erase(src.begin() + offset);
}


template<typename Ty>
qacc void decodeT(Ty& src)
{
  (void)src;
  return ;
}

template<>
qacc void decodeT(RealDD& src)
{
  RealDD a;
  a.Y() = src.Y() + src.X();
  a.X() = 0.0;
  a.X() = src - a;
  src = a;
}

template<>
qacc void decodeT(qlat::ComplexT<RealDD>& src)
{
  RealDD a = src.real();
  RealDD b = src.imag();
  decodeT(a);
  decodeT(b);
  src = qlat::ComplexT<RealDD>(a, b);
}

qacc Long Coordinate_to_index(const Coordinate& sp, const Coordinate& Lat){
  return ( ( Long(sp[3]) * Lat[2] + sp[2] ) * Lat[1] + sp[1] ) * Lat[0] + sp[0];
}

qacc Coordinate index_to_Coordinate(const Long& idx, const Coordinate& Lat){
  Coordinate sp;Long tem = idx;
  for(Int i=0;i<4;i++){
    sp[i] = tem % Lat[i];
    tem   = tem / Lat[i];
  }
  return sp;
}

//inline qlat::vector<Int > Coordinates_to_list(std::vector<Coordinate >& moms)
//{
//  qlat::vector<Int > momN;
//  momN.resize(moms.size() * 4);
//  for(unsigned int i=0;i<moms.size;i++){
//    for(Int j=0;j<4;j++){
//      momN[i*4 + j] = moms[i][j];
//    }
//  }
//  return momN;
//}

}

#endif

