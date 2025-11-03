// utils_io_vec.h
// Gen Wang
// Jan. 2021

#ifndef UTILS_IO_VEC_H
#define UTILS_IO_VEC_H
#pragma once

#include "general_funs.h"
#include "utils_fft_desc.h"
#include "utils_eo_copies.h"
#include "utils_props_type.h"

#define IO_DEFAULT  0
#define IO_ENDIAN false
#define IO_GN 1
#define IO_THREAD -1
#define __IO_SMALL_MEM__

////q_io_vec_ionum=io_num
////q_io_vec_thread=0

////Target, read double prop, double eigensystem
////read single prop, single eigensystem
////all in cps base,
////prop --> 12 *vol * d * c * complex, eigensytem, n * volv * d * c * complex
////Define gammas in PS, chiral, cps, base , ga
//

namespace qlat
{

struct io_vec
{
////public:

  std::vector<Int > node_ioL;
  Int ionum;
  Int nx,ny,nz,nt;
  Int Nx,Ny,Nz,Nt;
  size_t vol,noden;
  Int rank,Nmpi;
  std::vector<Int > nv,Nv;
  box<Geometry> geoB;
  Int threadio;
  std::vector<FILE* > file_omp;
  std::vector<size_t > currsend,currrecv,currspls,currrpls;
  Int MPI_size_c;
  std::vector<char > res;
  std::vector<char > tmp;
  std::vector<char > buf;
  ////char* tmp;
  bool do_checksum;
  std::vector<crc32_t > io_crc;
  crc32_t full_crc;

  size_t end_of_file;

  ///char* buf;size_t size_buf;

  //////only used with eigen system

  //////each noden has the memory order to t,z,y,x
  std::vector<std::vector<LInt > > map_Nitoi;
  ////fft_desc_basic fd;
  ///void get_pos(Int i);
  //FILE* io_read(const char *filename,const char* mode);
  //void io_close(FILE *file);

  inline FILE* io_read(const char *filename,const char* mode){
    Int do_thread_io = 0;
    std::string val = get_env(std::string("q_io_vec_thread"));
    if(val != ""){do_thread_io = stringtonum(val);}

    Int curr_threadio = threadio;
    Int Nv = omp_get_max_threads();
    if(curr_threadio < 1 or curr_threadio > Nv){curr_threadio = Nv;}
    if(do_thread_io == 0){curr_threadio = 1;}
    if(do_thread_io >  0){curr_threadio = do_thread_io;}
    if(curr_threadio * ionum > 64 ){curr_threadio = 1;}
    threadio = curr_threadio;

    if(do_checksum)ini_crc(true);
    if(node_ioL[rank]>=0){
      /////Currently can only open one file with openmp
      Qassert(file_omp.size() == 0);
      if(threadio!=1)for(Int i=0;i<threadio;i++){file_omp.push_back(fopen(filename, mode));}
      return fopen(filename, mode);
    }else{return NULL;}

  }

  // io geo always local
  inline const Geometry& geo(){
    //const Coordinate total_site = Coordinate(nx, ny, nz, nt);
    //return get_geo_cache(total_site);
    return geoB();
  }
  
  inline void io_close(FILE *file){
    if(do_checksum){sum_crc();end_of_file = 0;do_checksum=false;}
    if(node_ioL[rank]>=0){
      if(threadio!=1)for(Int i=0;i<threadio;i++){fclose(file_omp[i]);}
      file_omp.resize(0);
      fclose(file);
    }
    file = NULL;
    //if(buf != NULL){free(buf);size_buf = 0;buf = NULL;}
  }

  /////default initial file pos
  inline void io_off(FILE *file,size_t off, bool default_pos = true){
    if( default_pos){
      if(node_ioL[rank]>=0){
      fseek(file , off , SEEK_CUR );}
    }
    if(!default_pos){
      if(node_ioL[rank]>=0){
      fseek(file , off , SEEK_SET );}
    }
  }

  inline void ini_MPI(Int size_c0)
  {
    if(MPI_size_c == size_c0)return;
    MPI_size_c = size_c0;
    size_t size_c = size_c0 * noden;
    if(node_ioL[rank]>=0){
      ////if(tmp != NULL){delete []tmp;tmp=NULL;}
      ////tmp = new char[size_c0 * vol];
      //if(tmp != NULL){free(tmp);tmp=NULL;}
      //tmp = (char*) aligned_alloc_no_acc(size_c0 * vol);
      if(tmp.size() != size_c0 * vol){tmp.resize(size_c0 * vol);}
    }

    currsend.resize(Nmpi);
    currrecv.resize(Nmpi);
    currspls.resize(Nmpi);
    currrpls.resize(Nmpi);

    std::fill(currsend.begin(), currsend.end(), 0);
    std::fill(currrecv.begin(), currrecv.end(), 0);
    std::fill(currspls.begin(), currspls.end(), 0);
    std::fill(currrpls.begin(), currrpls.end(), 0);

    if(node_ioL[rank]>=0){
      for(Int n=0;n<Nmpi;n++){
        currsend[n] = size_c;
        currspls[n] = size_c*n;
      }
    }
    for(Int n=0;n<Nmpi;n++){
      if(node_ioL[n]>=0){
        currrecv[n] = size_c;
        currrpls[n] = size_c*node_ioL[n];
      }
    }
  }

  inline void ini_crc(bool new_file){
    io_crc.resize(Nmpi);
    for(LInt i=0;i<io_crc.size();i++){io_crc[i] = 0;}
    if(new_file)full_crc = 0;
  }
  inline crc32_t sum_crc(){
    sum_all_size(&io_crc[0], io_crc.size());
    for(LInt i=0;i<io_crc.size();i++){full_crc ^= io_crc[i];}
    return full_crc;
  }

  io_vec(){
    /////x,y,z,t
    threadio = IO_THREAD;
    MPI_size_c = 0;////tmp = NULL;
    do_checksum = 0;
    end_of_file = 0;

    nx = 0;ny = 0;nz = 0;nt = 0;
    Nx = 0;Ny = 0;Nz = 0;Nt = 0;

    vol   =  0;
    noden =  0;

    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();

    ////buf = NULL;size_buf = 0;
  }


  ////fd = &fds;
  io_vec(const Geometry& geo,Int ionum_or, Int threadio_set = IO_THREAD, bool do_checksum_set=false){
    TIMERA("Create io_vec");
    /////x,y,z,t
    threadio = threadio_set;
    MPI_size_c = 0;////tmp = NULL;
    do_checksum = do_checksum_set;
    end_of_file = 0;

    //int Nv = omp_get_max_threads();
    //if(threadio < 1 or threadio > Nv){threadio = Nv;}

    nv.resize(4);Nv.resize(4);
    for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
    Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];

    // set mapping to global geometry
    geoB.set_view(get_geo_cache(geo));

    vol   =  nx*ny*nz*nt;
    noden =  Nx*Ny*Nz*Nt;

    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();

    //for(Int ri=0;ri<Nmpi;ri++)MPI_Bcast(&map_Nitoi[ri][0], (noden/Nx)*sizeof(LInt), MPI_CHAR, ri, get_comm());
    //if(node_ioL[rank]==0){map_Nitoi[rank].resize(0);}

    node_ioL.resize(Nmpi);
    if(ionum_or > 0 and ionum_or < Nmpi){ionum=ionum_or;}else{
      std::string val = get_env(std::string("q_io_vec_ionum"));
      if(val == ""){ionum = 32;}else{ionum = stringtonum(val);}
      if(ionum < 0 or ionum > Nmpi){ionum = Nmpi;}
      ////qmessage("==io number %d \n", ionum);
    }

    for(Int ni=0;ni<Nmpi;ni++){node_ioL[ni]=-1;}
    /////distribute ranks out nodes
    {
      /////sorted id according to nodes
      std::vector<Int> id_to_node = qlat::get_id_node_in_shuffle_list();
      //int off = Nmpi/ionum;
      Int countN = 0;
      for(Int ni=0;ni<Nmpi;ni++){
        //if(ni%off==0 and countN < ionum)
        if(countN < ionum){
          ///node_ioL[ni]=countN;
          node_ioL[id_to_node[ni]]=countN;
          countN++;
        }
      }
    }

    ////fd = fft_desc_basic(geo); ///new ways
    #ifndef __IO_SMALL_MEM__
    const fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
    if(node_ioL[rank]>=0){
      map_Nitoi.resize(Nmpi);
      for(Int ri=0;ri<Nmpi;ri++){
        map_Nitoi[ri].resize(noden/Nx);

        #pragma omp parallel for
        for(size_t isp=0;isp<size_t(noden/Nx);isp++){
          //qlat::Coordinate ts = geo.coordinate_from_index(isp * Nx);
          //qlat::Coordinate gs = geo.coordinate_g_from_l(ts);
          //LInt offv = ((gs[3]*nz+gs[2])*ny+gs[1])*nx+gs[0];
          size_t offv = fd.index_g_from_local(isp*Nx +0 , ri);
          map_Nitoi[ri][isp] = offv;
        }
      }
    }
    #endif

    ///io_crc.resize(ionum);
    ini_crc(true);

    ////buf = NULL;size_buf = 0;
  }

  inline void clear_buf(){
    res.resize(0);////ionum*noden*gN*dsize, each nodes
    tmp.resize(0);////gN*dsize * vol, each io's
    buf.resize(0);////gN*dsize * vol, each io's
  }


  ~io_vec(){
    //if(node_ioL[rank]>=0){if(tmp != NULL){delete []tmp;tmp=NULL;}}
    //if(node_ioL[rank]>=0){if(tmp != NULL){free(tmp);tmp=NULL;}}
    map_Nitoi.resize(0);
    node_ioL.resize(0);
    currsend.resize(0);
    currrecv.resize(0);
    currspls.resize(0);
    currrpls.resize(0);
    io_crc.resize(0);
  }


};

/////io_vec buffers related
struct IOvecKey {
  Coordinate total_site;
  Int ionum;
  bool do_checksum_set;
  void init(const Coordinate& total_site_, Int ionum_ = 0, bool do_checksum_set_=false){
    total_site = total_site_;
    ionum = ionum_;
    do_checksum_set = do_checksum_set_;
  }

  IOvecKey(const Coordinate& total_site_, Int ionum_ = 0, bool do_checksum_set_=false)
  {
    init(total_site_, ionum_, do_checksum_set_);
  }

  IOvecKey(const Geometry& geo, Int ionum_ = 0, bool do_checksum_set_=false)
  {
    init(geo.total_site(), ionum_, do_checksum_set_);
  }
};

inline bool operator<(const IOvecKey& x, const IOvecKey& y)
{
  if(x.total_site < y.total_site ){  return true;}
  if(y.total_site < x.total_site ){  return false;}
  if(x.ionum < y.ionum ){  return true;}
  if(y.ionum < x.ionum ){  return false;}
  if(x.do_checksum_set < y.do_checksum_set ){  return true;}
  if(y.do_checksum_set < x.do_checksum_set ){  return false;}
  return false;
}

inline Cache<IOvecKey, io_vec >& get_io_vec_cache()
{
  static Cache<IOvecKey, io_vec > cache("IOvecKey", 16);
  return cache;
}

inline io_vec& get_io_vec_plan(const IOvecKey& fkey)
{
  if (!get_io_vec_cache().has(fkey)) {
    get_io_vec_cache()[fkey] = io_vec(fkey.total_site, fkey.ionum, IO_THREAD, fkey.do_checksum_set);
  }
  return get_io_vec_cache()[fkey];
}

inline io_vec& get_io_vec_plan(const Geometry& geo, Int ionum_ = 0)
{
  IOvecKey fkey(geo, ionum_);
  return get_io_vec_plan(fkey);
}

inline io_vec& get_io_vec_plan_with_checksum(const Geometry& geo)
{
  IOvecKey fkey(geo, 0, true);
  return get_io_vec_plan(fkey);
}

inline void clear_io_vec_cache()
{
  get_io_vec_cache().clear();
}

/////io_vec buffers related

inline void send_vec_kentucky(char* src,char* res,Int dsize,Int gN, io_vec& io, bool read=true)
{
  const std::vector<Int >& node_ioL = io.node_ioL;
  Int rank = io.rank;size_t vol = io.vol;size_t noden = io.noden;
  ///int Nmpi=io.Nmpi;

  Int Nx = io.Nx;

  Int size_c0 = gN*dsize;
  ////int size_c = gN*noden*dsize;

  #ifdef __IO_SMALL_MEM__
  const fft_desc_basic& fd = get_fft_desc_basic_plan(io.geo());
  #endif

  ////char* tmp=NULL;
  io.ini_MPI(size_c0);
  MPI_Barrier(get_comm());

  ///position p;
  if(read==true)
  {
  TIMERA("IO sort mem time");
  if(node_ioL[rank]>=0)
  for(Int ri=0;ri<io.Nmpi;ri++)for(Int gi=0;gi<gN;gi++){
      #pragma omp parallel for
      for(size_t isp=0;isp<size_t(noden/Nx);isp++){
        char* pres=&io.tmp[(ri*gN+gi)*noden*dsize];
        char* psrc=&src[gi*vol*dsize];
        #ifndef __IO_SMALL_MEM__
        ///size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
        size_t offv = io.map_Nitoi[ri][isp];
        #else
        size_t offv = fd.index_g_from_local(isp*Nx +0 , ri);
        #endif
        memcpy(&pres[isp*Nx*dsize],&psrc[offv*dsize],Nx*dsize);
      }
  }
  }

  {
  TIMERA("IO MPI time");
  /////for(Int io=0;io<node_ioL.size();io++)if(node_ioL[io]>=0)
  /////{
  /////  MPI_Scatter(tmp,size_c,MPI_CHAR,&res[node_ioL[io]*size_c],size_c,MPI_CHAR,io,get_comm());
  /////  ////if(rank==io){memcpy(&res[node_ioL[io]*size_c],tmp,size_c);}
  /////}


  //if(qlat::get_num_node() != 1)
  //{
  //if(read==true)
  //{MPI_Alltoallv(io.tmp,(int*) &io.currsend[0],(int*) &io.currspls[0], MPI_CHAR,
  //                  res,(int*) &io.currrecv[0],(int*) &io.currrpls[0], MPI_CHAR, get_comm());}
  //if(read==false)
  //{MPI_Alltoallv(   res,(int*) &io.currrecv[0],(int*) &io.currrpls[0], MPI_CHAR,
  //               io.tmp,(int*) &io.currsend[0],(int*) &io.currspls[0], MPI_CHAR, get_comm());}

  MPI_Comm tem_comm = get_comm();
  if(read==true)
  {MPI_Alltoallv_Send_Recv(
     (char*) io.tmp.data(), &io.currsend[0], &io.currspls[0], 
     (char*)    res       , &io.currrecv[0], &io.currrpls[0], tem_comm);}

  if(read==false)
  {MPI_Alltoallv_Send_Recv(
     (char*)    res       , &io.currrecv[0], &io.currrpls[0],
     (char*) io.tmp.data(), &io.currsend[0], &io.currspls[0], tem_comm);}

  //}
  //else{
  //  //memcpy(res,tmp, size_c);
  //  #pragma omp parallel for
  //  for(size_t isp=0;isp<size_c;isp++){res[isp] = tmp[isp];}
  //}

  }

  if(read == false)
  {
  TIMERA("IO sort mem time");
  if(node_ioL[rank]>=0)
  for(Int ri=0;ri<io.Nmpi;ri++)for(Int gi=0;gi<gN;gi++)
  {
    #pragma omp parallel for
    for(size_t isp=0;isp<size_t(noden/Nx);isp++){
      char* pres=&io.tmp[(ri*gN+gi)*noden*dsize];
      char* psrc=&src[gi*vol*dsize];
      #ifndef __IO_SMALL_MEM__
      //size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
      size_t offv = io.map_Nitoi[ri][isp];
      #else
      size_t offv = fd.index_g_from_local(isp*Nx +0 , ri);
      #endif
      memcpy(&psrc[offv*dsize],&pres[isp*Nx*dsize],Nx*dsize);
    }
  }
  }


}

////file for read/write,
////props read into pointer, 
////Nvec number of reading, Rendian change endian, 
////dsize inner size of vec can be 12*(8/4),
////single_file file format single/double,
////gN number of vec read for each io rank, read read/write of file
inline void read_kentucky_vector(FILE *file,char* props,Int Nvec,io_vec& io,bool Rendian=false, Int dsize=8, bool single_file=false,  Int gN=1, bool read = true)
{
  TIMER("IO VECS");

  timeval tm0,tm1,tm2,tm3;
  gettimeofday(&tm0, NULL);
  double mpi_t = 0.0;

  const std::vector<Int >& node_ioL = io.node_ioL;
  Int rank = io.rank;size_t vol = io.vol;int ionum=io.ionum;size_t noden=io.noden;

  if(io.do_checksum){
    if(io.end_of_file == 0){abort_r("io_vec need end of file for check sum! \n ");}
  }

  size_t size_buf = gN*vol*dsize;
  char* buf=NULL;
  ////char* buf=NULL;if(node_ioL[rank]>=0){buf = new char[gN*vol*dsize];}
  //if(buf.size() != size_buf){buf.resize(size_buf);}
  if(node_ioL[rank]>=0){
    if(io.buf.size() != size_buf){
      io.buf.resize(size_buf);
      //if(io.buf != NULL){free(io.buf);io.size_buf = 0;io.buf = NULL;}
      //io.buf = (char *)aligned_alloc_no_acc(size_buf);io.size_buf = size_buf;
      }
  }
  buf = io.buf.data();
  ////char res[ionum*noden*gN*dsize];
  if(io.res.size() != ionum*noden*gN*dsize){io.res.resize(ionum*noden*gN*dsize);}
  char* res = &io.res[0];
  size_t sizec = vol*dsize;

  size_t off_file = 0;////size_t sem=0;
  std::vector<size_t > semL;semL.resize(io.threadio);for(LInt i=0;i<semL.size();i++){semL[i]=0;}
  if(node_ioL[rank]>=0){off_file = ftell(file);}

  Int curr_v = 0;
  for(Int readi=0;readi<Nvec;readi++){

    /////From res to buf
    if (read == false) {
      TIMERA("IO copy mem time");
      for (Int iou = 0; iou < ionum; iou++)
        for (Int gi = 0; gi < gN; gi++)
          if (curr_v + iou * gN + gi < Nvec) {
            Int offv = curr_v + iou * gN + gi;
            // memcpy(&res[(iou*gN+gi)*noden*dsize +
            // 0],&props[(offv)*noden*dsize + 0],noden*dsize);

            char* pres = &res[(iou * gN + gi) * noden * dsize + 0];
            char* psrc = &props[(offv)*noden * dsize + 0];
#pragma omp parallel for
            for (size_t isp = 0; isp < size_t(noden * dsize); isp++) {
              pres[isp] = psrc[isp];
            }
          }
      gettimeofday(&tm2, NULL);
      send_vec_kentucky((char*)&buf[0], (char*)&res[0], dsize, gN, io, read);
      gettimeofday(&tm3, NULL);
      mpi_t += tm3.tv_sec - tm2.tv_sec;
      mpi_t += (tm3.tv_usec - tm2.tv_usec) / 1000000.0;
    }

    {
    TIMERA("IO disk");
    if(node_ioL[rank]>=0)
    if(curr_v + node_ioL[rank]*gN < Nvec)
    {
      size_t offv = curr_v + node_ioL[rank]*gN;
      fseek ( file , off_file + offv*sizec , SEEK_SET );
      size_t pos_file_cur = off_file + offv*sizec;

      Int rN = gN;if(offv+rN>=size_t(Nvec)){rN = Nvec-offv;}

      /////Switch endian of the file write
      if(read==false){
        Int dsize_single = sizeof(RealD);int fac = dsize/sizeof(RealD);
        if(single_file == true){dsize_single = sizeof(float);fac = dsize/sizeof(float);}
        //////Write Prop and link
        if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
        //////Write source noise
        if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
      }

      if(rN>0){
        Int curr_threadio = io.threadio;
        if(curr_threadio==1){
          semL[0] = file_operation(&buf[0], rN*sizec, 1, file, read);
          //if(read==true ){semL[0] =  fread(&buf[0], rN*sizec, 1, file);}
          //if(read==false){semL[0] = fwrite(&buf[0], rN*sizec, 1, file);}
          if(semL[0] != 1){printf("Reading/Writing error %zu 1 \n", semL[0] );}
        }


        ////omp reading
        if(curr_threadio>1){
          size_t Nvec = rN*sizec;
          size_t Group = (Nvec-1)/curr_threadio + 1;
          size_t off_file_local = ftell(file);
          #pragma omp parallel for
          for(Int tid=0;tid<curr_threadio;tid++)
          {
            //FILE *file_omp = file;
            /////FILE *fp2 = fdopen (dup (fileno (fp)), "r");
            size_t currN = Group; if((tid+1)*Group > Nvec){currN = Nvec - tid*Group;}
            if(currN > 0){
              size_t iniN  = tid*Group;//size_t endN = tid*Group + currN;
              //FILE *file_omp = fdopen (dup (fileno (file)), "r");
              fseek(io.file_omp[tid] , off_file_local + iniN, SEEK_SET );
              semL[tid] = file_operation(&buf[iniN], currN, 1, io.file_omp[tid], read);
              //if(read==true ){semL[tid] =  fread(&buf[iniN], currN, 1, io.file_omp[tid]);}
              //if(read==false){semL[tid] = fwrite(&buf[iniN], currN, 1, io.file_omp[tid]);}
            }
          }
          for(Int tid=0;tid<curr_threadio;tid++){if(semL[tid] != 1){printf("Reading/Writing error %zu 1 \n", semL[tid] );}}
          //fseek(file , rN*sizec , SEEK_CUR );
          fseek(file , off_file_local + rN*sizec , SEEK_SET );
        }

        ///////Get crc32 for buf (char *)aligned_alloc_no_acc(gN*vol*dsize)
        //crc32_z(initial, (const unsigned char*)data, );
        if(io.do_checksum){
          crc32_t crc32_tem = crc32_par((void*) buf, rN*sizec);
          /////crc32_shift(crc, off_file + offv*sizec);
          ////Shift current crc32 to the end of file
          crc32_tem = crc32_combine(crc32_tem, 0, io.end_of_file - pos_file_cur - rN*sizec);
          io.io_crc[rank] ^= crc32_tem;
        }

      }

      /////Switch endian of the file read
      if(read==true){
        Int dsize_single = sizeof(RealD);int fac = dsize/sizeof(RealD);
        if(single_file == true){dsize_single = sizeof(float);fac = dsize/sizeof(float);}
        //////Read Prop and link
        if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
        //////Read source noise
        if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
      }

    }
    }

    ////////Reordersrc buf
    ///if(gN!=1){if(node_ioL[rank]>=0)reorder_civ((char*)&buf[0],(char*)&buf[0],1,gN,vol,0 ,sizeof(RealD));}
    /////
      
    if(read==true){
      gettimeofday(&tm2, NULL);
      send_vec_kentucky((char*) &buf[0],(char*) &res[0], dsize,gN, io, read);
      gettimeofday(&tm3, NULL);
      mpi_t += tm3.tv_sec - tm2.tv_sec;
      mpi_t += (tm3.tv_usec - tm2.tv_usec)/1000000.0;
    }

    if(read==true)
    {
    TIMERA("IO copy mem time");
    for(Int iou=0;iou<ionum;iou++)
    for(Int gi=0;gi<gN;gi++)if(curr_v + iou*gN + gi<Nvec){
      Int offv = curr_v + iou*gN + gi;
      //memcpy(&props[(offv)*noden*dsize + 0],&res[(iou*gN+gi)*noden*dsize + 0],noden*dsize);

      char* pres=&props[(offv)*noden*dsize + 0];
      char* psrc=&res[(iou*gN+gi)*noden*dsize + 0];
      #pragma omp parallel for
      for(size_t isp=0;isp<size_t(noden*dsize);isp++){pres[isp] = psrc[isp];}
    }
    }
    curr_v += ionum*gN;
    if(curr_v >= Nvec)break;

  //fseek ( pFile , 9 , SEEK_SET );
  //fseek ( pFile , 9 , SEEK_CUR );
  //fseek ( pFile , 9 , SEEK_END );
  }

  if(node_ioL[rank]>=0){fseek ( file , off_file + Nvec*sizec , SEEK_SET );}
  ////if(node_ioL[rank]>=0){delete []buf;}
  ///if(node_ioL[rank]>=0){free(buf);}
  buf = NULL;

  gettimeofday(&tm1, NULL);
  double time0 = tm1.tv_sec - tm0.tv_sec;
  time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
  double total_R = vol*Nvec*dsize/(1024*1024*1024.0);

  if(read==true) qmessage("READ_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);
  if(read==false)qmessage("WRIT_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);

  MPI_Barrier(get_comm());
  fflush(stdout);
}

inline void save_txt_eigenvalues(std::vector<double > &values,std::vector<double > &errors,const char* filename, const char* sDescription)
{
  if(qlat::get_id_node() == 0)
  {
    Int nvec = values.size()/2;
    FILE* filew = fopen(filename, "w");

    fprintf(filew, "Eigenvalues and eigenvectors for the %s\n", sDescription);
    fprintf(filew, "Each eigenvector is preceded by a line describing the eigenvalue.\n");
    fprintf(filew, "The residue is defined as norm(mat.vec-lambda.vec).\n");
    fprintf(filew, "The format is: a tag EIGV, the real and imaginary part of the eigenvalue and the residue.\n");

    for(Int iv=0;iv<nvec;iv++)fprintf(filew, "EIGV %+.15le\t%+.15le\t%.10le\n", values[iv*2+0],values[iv*2+1],errors[iv]);

    fclose(filew);
  }
}

inline void save_txt_eigenvalues(std::vector<double > &values,std::vector<double > &errors,const std::string& filename, const std::string& sDescription)
{
  save_txt_eigenvalues(values, errors, filename.c_str(), sDescription.c_str());
}

inline void load_txt_eigenvalues(std::vector<double > &values,std::vector<double > &errors,const char* filename)
{
  Int nvec = 0;
  values.resize(0);errors.resize(0);
  if(qlat::get_id_node() == 0)
  {
    FILE* filer = fopen(filename, "r");
    char sTemp[300],tem[300];
    char* ftem = NULL;
    for(Int i=0;i<4;i++){
      ftem = fgets(sTemp, 300, filer);
      if(ftem == NULL){qmessage("Read eigenvalues error!\n");}
    }

    //while(!feof(filer))
    while(fgets(tem, 300, filer) != NULL){
      //fgets(tem, 300, filer);
      std::string re = std::string(tem);
      Int size = re.size();
      if(size > 8 and strcmp(re.substr(0,4).c_str(), "EIGV")==0){
        double v[3];
        sscanf(re.substr(4,size).c_str(),"%le %le %le",&v[0],&v[1],&v[2]);
        values.push_back(v[0]);
        values.push_back(v[1]);
        errors.push_back(v[2]);
        nvec += 1;
      }
      /////printf(tem,"");

    }

    fclose(filer);
  }

  MPI_Bcast(&nvec, 1, MPI_INT, 0, get_comm());
  if(nvec != 0 ){
  if(qlat::get_id_node() != 0){values.resize(nvec*2);errors.resize(nvec);}

  MPI_Bcast(&values[0], 2*nvec, MPI_DOUBLE, 0, get_comm());
  MPI_Bcast(&errors[0],   nvec, MPI_DOUBLE, 0, get_comm());
  }

}

inline void load_txt_eigenvalues(std::vector<double > &values,std::vector<double > &errors, const std::string& filename)
{
  load_txt_eigenvalues(values, errors, filename.c_str());
}

//
//order C, iDir, C_col, C_row, RealIm, t,z,y,x
  
template<typename Ty>
void rotate_gwu_vec_file(Ty* src,Int n_vec,size_t noden,bool single_file,bool read=true){
  TIMERB("Rotate gwu file vec");

  if(single_file == true){
    return;
  }

  //size_t noden = (*props[0]).desc->sites_on_node;
  size_t Np=noden*12*2;
  for(Int dc1=0;dc1<n_vec;dc1++){
    std::vector<Ty > p;Ty *q;
    p.resize(Np);
    //memcpy(&p[0],&src[dc1*Np],sizeof(Ty)*p.size());
    #pragma omp parallel for
    for(size_t isp=0;isp<p.size();isp++){p[isp] = src[dc1*Np + isp];}

    q = &src[dc1*Np];
    if(read==true){
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++){
      for(Int dc0=0;dc0<12;dc0++)
      {   
        q[isp*12*2 + dc0*2 + 0] = p[(0*12 + dc0)*noden + isp];
        q[isp*12*2 + dc0*2 + 1] = p[(1*12 + dc0)*noden + isp];
      }
    }}

    if(read==false){
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++){
      for(Int dc0=0;dc0<12;dc0++)
      {   
        q[(0*12 + dc0)*noden + isp] = p[isp*12*2 + dc0*2 + 0];
        q[(1*12 + dc0)*noden + isp] = p[isp*12*2 + dc0*2 + 1];
      }
    }}
  }
}

///////Ty is real as the transformation have only real numbers
template<typename Ty>
void gwu_to_cps_rotation_vec(Ty* src,Int n_vec,size_t noden,bool source=false,bool PS0=true,bool PS1=true, bool read=true){
  TIMERB("Rotate gwu file vec");
  ///prop
  if(source == true)if(n_vec%12!=0){qmessage("n_vec %8d.\n",n_vec);abort_r("source vector size wrong.\n");}
  //size_t noden = geo.node_site[0]*geo.node_site[1]*geo.node_site[2]*geo.node_site[3];

  Int signT = 1;if(read==false){signT=-1;}
  //std::vector<Typ0 > q;q.resize(12*2)
  const double sqrt2=std::sqrt(2.0);
  //////Rotate source dirac
  if(source == true){
    Int dr,d0,d1;Ty *q;
    Int Nprop = n_vec/12;
    
    size_t Np=noden*12*2;
    for(Int ip=0;ip<Nprop;ip++){
      q = &src[(ip*12+0)*Np];
      if(PS0 == true){
        std::vector<Ty > p;
        p.resize(12*Np);
        //memcpy(&p[0],&src[(ip*12+0)*Np],sizeof(Ty)*p.size());
        #pragma omp parallel for
        for(size_t isp=0;isp<p.size();isp++){p[isp] = src[(ip*12+0)*Np + isp];}


        dr=0;d0=1;d1=3;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = (-p[d0*3*Np+c]*signT + p[d1*3*Np+c])/sqrt2;
        dr=1;d0=0;d1=2;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = ( p[d0*3*Np+c]*signT - p[d1*3*Np+c])/sqrt2;
        dr=2;d0=1;d1=3;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = (-p[d0*3*Np+c] - p[d1*3*Np+c]*signT)/sqrt2;
        dr=3;d0=0;d1=2;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = (+p[d0*3*Np+c] + p[d1*3*Np+c]*signT)/sqrt2;
      }else{
        dr = 2;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = -1.0*q[dr*3*Np+c];
        dr = 3;
        #pragma omp parallel for
        for(size_t c=0;c<3*Np;c++)q[dr*3*Np+c] = -1.0*q[dr*3*Np+c];
      }
    }
  }

  /////Rotate sink dirac
  for(Int dc1=0;dc1<n_vec;dc1++){
  #pragma omp parallel for
  for(size_t isp=0;isp<noden;isp++){
    Int dr,d0,d1;Ty *q;
    std::vector<Ty > p;p.resize(12*2);
    q = &src[(dc1*noden + isp)*12*2 + 0];

    if(PS1 == true){
      memcpy(&p[0],&src[(dc1*noden + isp)*12*2 + 0],sizeof(Ty)*p.size());
      dr=0;d0=1;d1=3;for(Int c=0;c<6;c++)q[dr*6+c] = (-p[d0*6+c]*signT + p[d1*6+c])/sqrt2;
      dr=1;d0=0;d1=2;for(Int c=0;c<6;c++)q[dr*6+c] = ( p[d0*6+c]*signT - p[d1*6+c])/sqrt2;
      dr=2;d0=1;d1=3;for(Int c=0;c<6;c++)q[dr*6+c] = (-p[d0*6+c] - p[d1*6+c]*signT)/sqrt2;
      dr=3;d0=0;d1=2;for(Int c=0;c<6;c++)q[dr*6+c] = (+p[d0*6+c] + p[d1*6+c]*signT)/sqrt2;
    }else{
      dr = 2;for(Int c=0;c<6;c++)q[dr*6+c] = -1.0*q[dr*6+c];
      dr = 3;for(Int c=0;c<6;c++)q[dr*6+c] = -1.0*q[dr*6+c];
    }
  }}

}

template<typename Ty>
Ty get_norm_vec(Ty *src,size_t noden){
  Ty res = 0.0;

  ComplexT<Ty > *p = (ComplexT<Ty >*) src;
  /////need sum 12 first to reduce float sum error
  #pragma omp parallel for reduction(+: res)
  for(size_t isp=0;isp<noden;isp++)
  {
    Ty a = 0.0;
    for(unsigned int dc=0;dc<12;dc++){
       a += qnorm(p[isp*12+dc]);
    }
    res += a;
  }

  //qmessage("==omp_get_max_threads %8d \n ",omp_get_max_threads());
  sum_all_size(&res,1);
  return res;

}

inline Int test_single(const char *filename,io_vec &io_use,Int iv=0){
  
  double normd=0.0;double normf=0.0;
  Int n_vec = 1;
  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  FILE* file;
  double err = 1e-3;
  bool do_checksum_buf = io_use.do_checksum;

  {
  io_use.do_checksum = false;
  std::vector<double > prop_E;
  prop_E.resize(n_vec*12*noden*2);

  file = io_use.io_read(filename,"rb");
  io_use.io_off(file,iv*Fsize*2, false);
  read_kentucky_vector(file,(char*) &prop_E[0], n_vec*12*2,io_use, false, sizeof(RealD), false, 12*2);
  io_use.io_close(file);

  rotate_gwu_vec_file(&prop_E[0],n_vec,noden, false);
  gwu_to_cps_rotation_vec(&prop_E[0],n_vec,noden, false, true, true);

  normd = get_norm_vec(&prop_E[0],noden);
  }

  io_use.do_checksum = do_checksum_buf;
  if(fabs(normd - 1.0) < err)return 0;

  {
  io_use.do_checksum = false;
  std::vector<RealF > prop_E;
  prop_E.resize(n_vec*12*noden*2);
  file = io_use.io_read(filename,"rb");
  io_use.io_off(file,iv*Fsize, false);

  read_kentucky_vector(file,(char*) &prop_E[0], n_vec,io_use, true, 6*2*sizeof(RealD), true);
  io_use.io_close(file);

  rotate_gwu_vec_file(&prop_E[0],n_vec,noden, true);
  gwu_to_cps_rotation_vec(&prop_E[0],n_vec,noden, false, false, false);

  normf = get_norm_vec(&prop_E[0],noden);
  }


  io_use.do_checksum = do_checksum_buf;
  if(fabs(normf - 1.0) < err)return 1;

  qmessage("Norm of vector double %.3e, %.3e.\n",normd,normd-1.0);
  qmessage("Norm of vector single %.3e, %.3e.\n",normf,normf-1.0);
  return -1;
}
/////================END of NO checksum read/write

inline Int check_Eigen_file_type(const char *filename, io_vec &io_use,Int n1,bool check){

  size_t noden = io_use.noden;
  size_t sizen = get_file_size_MPI(filename);

  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  if(sizen <  Fsize){qmessage("%s \n",filename);abort_r("Eigen file size too small. \n");}
  if(sizen <2*Fsize){
    if(n1 != 1){abort_r("Eigen file size too small 1. \n");}
    /////read with single
    return 1;
  }

  ////Default use single
  if(check == false){
    if(sizen < n1*Fsize){abort_r("Eigen file size too small 2. \n");}
    return 1;
  }

  Int tmp = test_single(filename,io_use,0);
  if(tmp == -1){abort_r("Eigen system file norm not correct.\n");}

  if(tmp == 1){
    if(sizen <  n1*Fsize){abort_r("Eigen file size too small 2. \n");}
    return 1 ;
  }
  if(tmp == 0){
    if(sizen <2*n1*Fsize){abort_r("Eigen file size too small 3. \n");}
    return 0 ;
  }

  return 0;
}


/////Read the data into the point resp which have memory allocated already
/////check --> if false, abort if file is double,
/////read --> flag for writing and reading, may be can only write single currently
template<typename Ty>
void load_gwu_eigen(FILE* file,std::vector<Ty* > resp,io_vec &io_use,Int n0,Int n1,
  bool check=true, bool read=true, bool read_single=true)
{
  (void)check;
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(resp.size() < size_t(n1-n0)){abort_r("Final point size wrong!\n");}
  if(sizeof(Ty) != sizeof(RealD) and sizeof(Ty) != sizeof(float)){abort_r("Need double or float pointer! \n");}
  ////if(n0!=0 and read==false){abort_r("Write vectors with off not zero!");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  bool single = read_single;
  //bool single = true;
  //////read==false, write only single vectors
  //if(read==true){
  //  Int type = check_Eigen_file_type(filename,io_use,n1,check);
  //  if(type == 0){single = false;}
  //  if(type == 1){single = true ;}
  //  if(type !=0 and type != 1){abort_r("Eigen system not in gwu format! \n");}
  //}
  ////FILE* file;

  if(single == true){
    /////Load single precision
      
    Int count = 0;int off = io_use.ionum;
    std::vector<RealF > prop_E;
    Int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    ////if(read==true )file = io_use.io_read(filename,"rb");
    ////if(read==false)file = io_use.io_read(filename,"wb");
    io_use.io_off(file, n0*Fsize, true);

    for(Int iv=0;iv<n_vec;iv++){
      Int ri = off;if(count + off > n_vec){ri = n_vec - count;}
      //int offE = count*12*noden*2;

      if (read == false) {
        for (Int ic = 0; ic < ri; ic++)
#pragma omp parallel for
          for (size_t isp = 0; isp < noden * 12 * 2; isp++) {
            prop_E[ic * noden * 12 * 2 + isp] = resp[count + ic][isp];
          }
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0], ri, noden, false, false, false,
                                read);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0], ri, noden, true, read);
      }

      /////qmessage("WRITE!!!!!!\n");
      read_kentucky_vector(file,(char*) &prop_E[0], ri,io_use, true, 6*2*sizeof(RealD), true, 1, read);
      fflush_MPI();

      if(read==true){
        //rotate_qlat_to_gwu(prop_E, &src_new.vec[0],geo, true);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, true, read);
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, false, false, read);
        for(Int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
          resp[count + ic][isp] = prop_E[ic*noden*12*2 + isp];}
      }

      count += ri;if(count >= n_vec){break;}
    }
    /////io_use.io_close(file);

    /////rotate_qlat_to_gwu(prop_E,&src_new.vec[0],geo);
    /////////From ky to milc
    /////for(Int iv=0;iv<n_vec;iv++)ga.ga[4][0].multiply(*src_new.vec[iv],*src_new.vec[iv]);
  }else{
    /////Load double precision

    Int count = 0;int off = io_use.ionum;
    std::vector<double > prop_E;
    Int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    //////file = io_use.io_read(filename,"rb");
    io_use.io_off(file,n0*Fsize*2, true);
      
    for(Int iv=0;iv<n_vec;iv++){
      Int ri = off;if(count + off > n_vec){ri = n_vec - count;}
      //int offE = count*12*noden*2;

      if (read == false) {
        for (Int ic = 0; ic < ri; ic++)
#pragma omp parallel for
          for (size_t isp = 0; isp < noden * 12 * 2; isp++) {
            prop_E[ic * noden * 12 * 2 + isp] = resp[count + ic][isp];
          }
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0], ri, noden, false, true, true, read);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0], ri, noden, false, read);
      }

      read_kentucky_vector(file,(char*) &prop_E[0], ri*12*2,io_use, false, sizeof(RealD), false, 12*2, read);
      fflush_MPI();

      if(read==true){
        /////double precision eigen vector in ps base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, false, read);
        ////Do not rotate source, rotate sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, true, true, read);

        for(Int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
          resp[count + ic][isp] = prop_E[ic*noden*12*2 + isp];
        }
      }

      count += ri;if(count >= n_vec){break;}
    }

    //////io_use.io_close(file);

  }

}

//template<typename Ty>
//void load_gwu_eigen(FILE* file, std::vector<qlat::FieldM<Ty , 12> > res,io_vec &io_use,Int n0,Int n1,
//  bool check=true, bool read=true, bool read_single=true)
//{
//  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
//  Int n_vec = n1-n0;
//  size_t noden = io_use.noden;
//  if(read == true)
//  if(res.size() != n_vec)
//  {
//    res.resize(0);
//    res.resize(n_vec);
//    for(Int iv=0;iv<n_vec;iv++){res[iv].init(io_use.geo());}
//  }
//  std::vector<Ty* > resp;resp.resize(n_vec);
//  for(Int iv=0;iv<n_vec;iv++){resp[iv] = (Ty*) qlat::get_data(res[iv]).data();}
//  load_gwu_eigen(file, resp, io_use,n0,n1,check, read, read_single);
//}

inline FILE* open_gwu_eigen(const char *filename,io_vec &io_use, bool read=true)
{
  if(read==true ){return io_use.io_read(filename,"rb");}
  if(read==false){return io_use.io_read(filename,"wb");}
  return NULL;
}

template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty* > resp,io_vec &io_use,Int n0,Int n1,
  bool check=true, bool read=true)
{
  /////FILE* file = NULL;
  bool read_single = true;
  if(read==true){read_single = check_Eigen_file_type(filename, io_use, n1, check);}
  FILE* file = open_gwu_eigen(filename, io_use, read);
  load_gwu_eigen(file, resp, io_use, n0, n1, check, read, read_single);
  io_use.io_close(file);
}

template<typename Ty>
void load_gwu_eigen(FILE *file,std::vector<Ty > &res,io_vec &io_use,Int n0,Int n1,bool check=true,bool read=true, bool read_single=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  Int n_vec = n1-n0;
  size_t noden = io_use.noden;
  if(read == true)res.resize(n_vec*noden*12*2);

  std::vector<Ty* > resp;resp.resize(n_vec);
  for(Int iv=0;iv<n_vec;iv++){resp[iv] = &res[iv*noden*12*2];}
  load_gwu_eigen(file, resp, io_use,n0,n1,check,read, read_single);
}

//////Ty should be float or double
template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty > &res,io_vec &io_use,Int n0,Int n1,bool check=true,bool read=true){
  /////FILE* file = NULL;
  bool read_single = true;
  if(read==true){read_single = check_Eigen_file_type(filename, io_use, n1, check);}
  FILE* file = open_gwu_eigen(filename, io_use, read);
  load_gwu_eigen(file, res, io_use,n0,n1,check,read, read_single);
  io_use.io_close(file);
}


template<typename Td>
void load_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Td> > &eigen,io_vec &io_use,Int n0,Int n1,bool check=true,bool read=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(sizeof(Td) != sizeof(double ) and sizeof(Td) != sizeof(float )){abort_r("Cannot understand the input format! \n");}

  Int n_vec = n1-n0;
  if(read==true){
    eigen.resize(0);
    eigen.resize(n_vec);
    for(Int iv=0;iv<n_vec;iv++)eigen[iv].init(io_use.geo());
  }

  if(sizeof(Td) == sizeof(double )){
    std::vector<double* > resp;resp.resize(n_vec);
    for(Int iv=0;iv<n_vec;iv++){
      resp[iv]=(double*)(qlat::get_data(eigen[iv]).data());
    }
    load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
  }

  if(sizeof(Td) == sizeof(float )){
    std::vector<RealF* > resp;resp.resize(n_vec);
    for(Int iv=0;iv<n_vec;iv++){
      resp[iv]=(RealF*)(qlat::get_data(eigen[iv]).data());
    }
    load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
  }

  ////return eigen;
}

template<typename Td>
void save_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Td> > &eigen,io_vec &io_use,Int n0,Int n1,bool check=true){
  load_gwu_eigen(filename,eigen,io_use,n0,n1,check,false);
}

inline void load_gwu_eigen(const char *filename,EigenM &Mvec,io_vec &io_use,Int n0,Int n1,bool check=true,bool read=true)
{
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  Int n_vec = n1-n0;
  if(n_vec > int(Mvec.size())){abort_r("Read number of eigen larger than memory. \n");}

  Long Nsize = io_use.noden*12;

  std::vector<Ftype* > resp;resp.resize(n_vec);
  for(Int iv=0;iv<n_vec;iv++){
    if(Mvec[iv].size() != Nsize){
      if(read==false)abort_r("Mvec not initialized! \n");
      if(read==true ){Mvec[iv].resize(Nsize);set_zero(Mvec[iv]);}
    }
    resp[iv]= (Ftype*)(&Mvec[iv][0]);
  }
  load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
}

inline void save_gwu_eigen(const char *filename,EigenM &Mvec,io_vec &io_use,Int n0,Int n1,bool check=true)
{
  load_gwu_eigen(filename, Mvec,io_use,n0,n1,check,false);
}


template<typename Td>
void load_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Td> > &prop,io_vec &io_use,bool read=true){

  if(sizeof(Td) != sizeof(double ) and sizeof(Td) != sizeof(float )){abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);

  FILE* file;bool single = true;
  if(read==true){
  size_t sizen = get_file_size_MPI(filename);

  if(sizen != 2*Fsize*12 and sizen != Fsize*12){qmessage("File %s \n",filename);abort_r("prop size wrong! \n");}
  prop.resize(0);
  prop.resize(12);
  for(Int iv=0;iv<12;iv++)prop[iv].init(io_use.geo());
  if(sizen == 2*Fsize*12){single=false;}
  }
    
  ////Can only write with single!
  if(read==false){single = true;}
  
  if(single == true)
  {
    ////Single vector read
    std::vector<RealF > prop_qlat;
    prop_qlat.resize(12*noden*12*2);

    if(read==false){
    for(Int iv=0;iv<12;iv++){
      qlat::ComplexT<Td>* res   = (qlat::ComplexT<Td>*) qlat::get_data(prop[iv]).data();
      ComplexT<RealF> *src = (ComplexT<RealF>*) &prop_qlat[iv*noden*12*2];
      #pragma omp parallel for
      for(size_t isp=0;isp<noden*12;isp++)src[isp] = ComplexT<RealF>(res[isp].real(),res[isp].imag());
    }

    ////Do not rotate source, in ps/ky base
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true, true, read);
    /////double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, true, read);
    }

    if(read==true )file = io_use.io_read(filename,"rb");
    if(read==false)file = io_use.io_read(filename,"wb");
    read_kentucky_vector(file,(char*) &prop_qlat[0], 12,io_use, true, 6*2*sizeof(RealD), true, 1 , read);
    io_use.io_close(file);
      
    if(read==true){
    /////double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, true);
    ////Do not rotate source, in ps/ky base
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true, true);

    for(Int iv=0;iv<12;iv++){
      qlat::ComplexT<Td>* res   = (qlat::ComplexT<Td>*) qlat::get_data(prop[iv]).data();
      ComplexT<RealF> *src = (ComplexT<RealF>*) &prop_qlat[iv*noden*12*2];
      #pragma omp parallel for
      for(size_t isp=0;isp<noden*12;isp++)res[isp]= qlat::ComplexT<Td>(src[isp].real(),src[isp].imag());
    }
    }
  }

  if(single == false)
  {
    std::vector<double > prop_qlat;
    prop_qlat.resize(12*noden*12*2);

    file = io_use.io_read(filename,"rb");
    read_kentucky_vector(file,(char*) &prop_qlat[0], 12*2*12,io_use, false, sizeof(RealD), false, 1);
    io_use.io_close(file);
      
    ///double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, false);
    //Do not rotate source, 
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true,true);

    for(Int iv=0;iv<12;iv++){
      qlat::ComplexT<Td>* res = (qlat::ComplexT<Td>*) qlat::get_data(prop[iv]).data();
      ComplexT<double> *src = (ComplexT<double>*) &prop_qlat[iv*noden*12*2];
      for(size_t isp=0;isp<noden*12;isp++)res[isp]= qlat::ComplexT<Td>(src[isp].real(),src[isp].imag());
    }

  }
}

template<typename Td>
void save_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Td> > &prop,io_vec &io_use){
  load_gwu_prop(filename,prop,io_use,false);
}

//////final result 12*12 --> Nt*Nxyz
template<typename Td>
void load_gwu_prop(const char *filename, qlat::FieldM<qlat::ComplexT<Td>, 12*12>& res,io_vec &io_use,bool read=true){
  const Geometry& geo = io_use.geo();
  if(read == true ){res.init();res.init(geo);}
  if(read == false){abort_r("Not supported! \n");}

  Long sizeF = geo.local_volume();

  std::vector<qlat::FermionField4dT<Td> > prop;
  load_gwu_prop(filename, prop, io_use, read);

  move_index mv_civ;
  qlat::ComplexT<Td>* p0; qlat::ComplexT<Td>* p1;qlat::ComplexT<Td>* pt;
  pt = (qlat::ComplexT<Td>*) qlat::get_data(res).data();

  for(Int iv=0;iv<12;iv++){
    p0 = (qlat::ComplexT<Td>*) qlat::get_data(prop[iv]).data();
    p1 = (qlat::ComplexT<Td>*) &pt[iv*12*sizeF + 0];
    mv_civ.dojob(p0, p1, 1, 12, sizeF, 1, 1, false);
  }
}


template <typename Td>
void save_gwu_prop(const char *filename,Propagator4dT<Td>& prop){
  Qassert(prop.initialized);
  io_vec& io_use = get_io_vec_plan(prop.geo());
  std::vector<qlat::FermionField4dT<Td > > prop_qlat;
  prop4d_to_Fermion(prop_qlat, prop, 1);
  save_gwu_prop(filename,prop_qlat,io_use);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <typename Td>
void save_gwu_prop(std::string &filename,Propagator4dT<Td>& prop){
  // char tem[500];
  // ssprintf(tem,filename.c_str());
  // save_gwu_prop(tem,prop);
  save_gwu_prop(filename.c_str(),prop);
}

template <typename Td>
void load_gwu_prop(const char *filename,Propagator4dT<Td>& prop){
  Qassert(prop.initialized);
  io_vec& io_use = get_io_vec_plan(prop.geo());
  std::vector<qlat::FermionField4dT<Td > > prop_qlat;
  load_gwu_prop(filename,prop_qlat,io_use);
  prop4d_to_Fermion(prop_qlat, prop, 0);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <typename Td>
void load_gwu_prop(std::string &filename,Propagator4dT<Td>& prop){
  // char tem[500];
  // ssprintf(tem,filename.c_str());
  // load_gwu_prop(tem,prop);
  load_gwu_prop(filename.c_str(),prop);
}

template <typename Td>
void load_gwu_link(const char *filename,GaugeFieldT<Td> &gf, bool read = true){
  io_vec io_use(gf.geo(),8);
  //if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ))
  //{abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(4*9*noden*2)*sizeof(RealD);

  std::vector<double > link_qlat;
  link_qlat.resize(4*9*noden*2);

  if(read==true)
  {
    size_t sizen = get_file_size_MPI(filename);
    if(sizen != Fsize){abort_r("Link size wrong! \n");}
  }


  if(read == false)
    for (size_t index = 0; index < noden; ++index) {
      ColorMatrixT<Td>& res = gf.get_elem_offset(index * gf.multiplicity + 0);

      for (Int dir = 0; dir < 4; dir++) {
        for (Int c0 = 0; c0 < 3; c0++)
          for (Int c1 = 0; c1 < 3; c1++) {
            Int dir0 = ((dir * 3 + c0) * 3 + c1);
            Int dir1R = ((dir * 3 + c1) * 3 + c0) * 2 + 0;
            Int dir1I = ((dir * 3 + c1) * 3 + c0) * 2 + 1;
            link_qlat[dir1R * noden + index] = res.p[dir0].real();
            link_qlat[dir1I * noden + index] = res.p[dir0].imag();
          }
      }
    }

  FILE* file;
  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  read_kentucky_vector(file,(char*) &link_qlat[0], 4*9*2,io_use, false, sizeof(RealD), false,9*2, read);
  io_use.io_close(file);
    
  //////double precision eigen vector in ps base
  ///rotate_gwu_vec_file(&prop_qlat[0], 12,noden, false);
  /////Do not rotate source, 
  ///gwu_to_cps_rotation_vec(&link_qlat[0], 12,noden, true, true,true);
  ///4 dire --> c0 -- > c1

  if(read == true)
  for (size_t index = 0; index < noden; ++index)
  {
    ColorMatrixT<Td>& res = gf.get_elem_offset(index*gf.multiplicity+0);

    for(Int dir=0;dir<4;dir++)
    {
      for(Int c0=0;c0<3;c0++)
      for(Int c1=0;c1<3;c1++)
      {
        Int dir0  = ((dir*3+c0)*3+c1)    ;
        Int dir1R = ((dir*3+c1)*3+c0)*2+0;
        Int dir1I = ((dir*3+c1)*3+c0)*2+1;
        res.p[dir0] = qlat::ComplexT<Td>(link_qlat[dir1R*noden + index], link_qlat[dir1I*noden + index]);
      }
    }
  }
}

template <class Td>
void save_gwu_link(const char *filename,GaugeFieldT<Td> &gf){
  load_gwu_link(filename, gf, false);
}

template <class Td>
void load_gwu_link(std::string &filename,GaugeFieldT<Td> &gf){
  // char tem[500];
  // ssprintf(tem,filename.c_str());
  // load_gwu_link(tem,gf);
  load_gwu_link(filename.c_str(),gf);
}

template<typename Ty>
void load_gwu_noies(const char *filename,std::vector<qlat::FieldM<Ty, 1> > &noises,bool read=true){

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ) and IsTypeComplex<Ty>() == 0){abort_r("Cannot understand the input format! \n");}
  if(noises.size() == 0 and read == false){return ;}
  if(noises.size() == 0 and read == true){abort_r("Initialize your vectors");}

  const Long Nnoi = noises.size();

  Qassert(noises[0].initialized);
  io_vec& io_use = get_io_vec_plan(noises[0].geo());

  FILE* file;
  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*2)*sizeof(float);

  std::vector<double > prop_noi;
  prop_noi.resize(noden*2);
  move_index mv_civ;

  for(Long iv=0;iv<Nnoi;iv++)
  {
    qlat::FieldM<Ty, 1>& noi = noises[iv];
    if(read==true){
      if(iv == 0){
        size_t sizen = get_file_size_MPI(filename);
        if(sizen != Nnoi*2*Fsize){abort_r("noise size wrong! \n");}
      }
      const Geometry& geo = io_use.geo();// geo.multiplicity=1;
      if(!noi.initialized)noi.init(geo);
    }
    if(read==false){
      ComplexT<double> *src = (ComplexT<double>*) &prop_noi[0];
      Ty* res = (Ty*) qlat::get_data(noi).data();
      #pragma omp parallel for
      for(size_t isp=0;isp<noden;isp++){
        src[isp] = ComplexT<double>(res[isp].real(),res[isp].imag());
      }
      mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(RealD), false);
    }

    read_kentucky_vector(file,(char*) &prop_noi[0], 2,io_use, true, sizeof(RealD), false, 1,read);

    /////Copy noise vectors
    if(read==true){
      //reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD));
      //io_use.mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD), false);
      mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD), false);
      ComplexT<double> *src = (ComplexT<double>*) &prop_noi[0];
      Ty* res = (Ty*) qlat::get_data(noi).data();
      #pragma omp parallel for
      for(size_t isp=0;isp<noden;isp++){
        res[isp]= Ty(src[isp].real(),src[isp].imag());
      }
    }
  }

  io_use.io_close(file);
}

template<typename Ty>
void save_gwu_noies(const char *filename,std::vector<qlat::FieldM<Ty,1>> &noi){
  load_gwu_noies(filename,noi,false);
}


template<class Fieldy>
void load_gwu_noiP(const char *filename, Fieldy& noi,bool read=true, bool GPU = false){
  TIMER("load_gwu_noi");
  Qassert(GetBasicDataType<Fieldy>::get_type_name() != std::string("unknown_type"));
  using D = typename GetBasicDataType<Fieldy>::ElementaryType;
  Qassert(IsBasicTypeReal<D>());
  using Ty = ComplexT<D >;

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ) and IsTypeComplex<Ty>() == 0){abort_r("Cannot understand the input format! \n");}

  Qassert(noi.initialized and noi.multiplicity == 1);
  io_vec& io_use = get_io_vec_plan(noi.geo());

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*2)*sizeof(float);

  std::vector<double > prop_noi;
  prop_noi.resize(noden*2);
  move_index mv_civ;

  if(read==true){
  size_t sizen = get_file_size_MPI(filename);
  if(sizen != 2*Fsize){abort_r("noise size wrong! \n");}
  }
  if(read==false){
    ComplexT<double> *src = (ComplexT<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    cpy_GPU(src, res, noden, 0, GPU);
    //#pragma omp parallel for
    //for(size_t isp=0;isp<noden;isp++)src[isp] = ComplexT<double>(res[isp].real(),res[isp].imag());
    //reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(RealD));
    //io_use.mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(RealD), false);
    mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(RealD), false);
  }

  FILE* file;

  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  read_kentucky_vector(file,(char*) &prop_noi[0], 2,io_use, true, sizeof(RealD), false, 1,read);
  io_use.io_close(file);

  /////Copy noise vectors
  if(read==true){
    //reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD));
    //io_use.mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD), false);
    mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(RealD), false);
    ComplexT<double> *src = (ComplexT<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    cpy_GPU(res, src, noden, GPU, 0);
    //#pragma omp parallel for
    //for(size_t isp=0;isp<noden;isp++)res[isp]= Ty(src[isp].real(),src[isp].imag());
  }
}

template<typename Ty>
void load_gwu_noiG(const char *filename,qlat::FieldG<Ty> &noi){
  Qassert(noi.initialized and noi.multiplicity == 1);
  load_gwu_noiP(filename,noi,true, 1);
}

template<typename Ty>
void save_gwu_noiG(const char *filename,qlat::FieldG<Ty> &noi){
  Qassert(noi.initialized and noi.multiplicity == 1);
  load_gwu_noiP(filename,noi,false, 1);
}

template<typename Ty>
void load_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi){
  load_gwu_noiP(filename,noi,true);
}

template<typename Ty>
void save_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi){
  load_gwu_noiP(filename,noi,false);
}

template <typename Td>
void save_gwu_noiP(const char *filename,Propagator4dT<Td>& prop){
  qlat::FieldM<qlat::ComplexD,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  
  const long noden = prop.geo().local_volume();
  qthread_for(index, noden, {
    qlat::WilsonMatrixT<Td>&  src =  prop.get_elem_offset(index);
    double sum = 0.0;
    for(Int d1=0;d1<12;d1++)
    for(Int d0=0;d0<12;d0++)
    {
      sum += std::fabs(src(d1,d0).real());
      sum += std::fabs(src(d1,d0).imag());
    }
    qlat::ComplexD phase = qlat::ComplexD(src(0,0).real(),src(0,0).imag());

    if(sum >1e-8){noi.get_elem_offset(index) = 1.0*phase;}
  });

  save_gwu_noi(filename,noi);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <typename Td>
void save_gwu_noiP(std::string &filename,Propagator4dT<Td>& prop){
  // char tem[500];
  // ssprintf(tem,filename.c_str());
  // save_gwu_noiP(tem,prop);
  save_gwu_noiP(filename.c_str(),prop);
}

template <typename Td>
void noi_to_propP(qlat::FieldM<qlat::ComplexD,1> &noi,Propagator4dT<Td>& prop, Int dir = 0){
  for (Long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrixT<Td>& res =  prop.get_elem_offset(index);
    if(dir==0)for(Int d0=0;d0<12;d0++){res(d0,d0) = noi.get_elem_offset(index);}
    if(dir==1)for(Int d0=0;d0<12;d0++){noi.get_elem_offset(index) = res(d0,d0);}
  }
}

template <typename Td>
void load_gwu_noiP(const char *filename,Propagator4dT<Td>& prop){
  ////io_vec& io_use = get_io_vec_plan(prop.geo());
  qlat::FieldM<qlat::ComplexD,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  load_gwu_noi(filename,noi);
  prop.init(noi.geo());

  noi_to_propP(noi, prop, 0);
  
}

template <typename Td>
void load_gwu_noiP(std::string &filename,Propagator4dT<Td>& prop){
  // char tem[500];
  // ssprintf(tem,filename.c_str());
  // load_gwu_noiP(tem,prop);
  load_gwu_noiP(filename.c_str(),prop);
}

/////================END of NO checksum read/write

//////Assume memory allocated already
template<class T, typename Ty>
void copy_noise_to_vec(T* noiP, Ty* buf, const Geometry& geo, const Long bfac, const Int dir=1, const QMEM Gmem = QMCPU)
{
  TIMERB("copy_noise_to_vec");
  const Long Nd = geo.local_volume() * bfac;
  if(dir == 1){cpy_GPU(buf, noiP, Nd, QMCPU, Gmem);}
  if(dir == 0){cpy_GPU(noiP, buf, Nd, Gmem, QMCPU);}

  //#pragma omp parallel for 
  //for (Long index = 0; index < geo.local_volume(); ++index)
  //{
  //  for(Int bi=0;bi<bfac;bi++){
  //    if(dir==1){buf[index*bfac+bi] = noiP[index*bfac+bi];}
  //    if(dir==0){noiP[index*bfac+bi]= buf[index*bfac+bi];}
  //  }
  //}
}

inline void open_file_qlat_noisesT(const char *filename, Int bfac, inputpara& in, bool read=true, bool single_file=true, Int N_noi=-1, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), bool rotate_bfac = true)
{
  in.bsize = sizeof(float );
  in.rotate_bfac = rotate_bfac;
  in.do_checksum = true;
  in.N_noi = N_noi; in.ncur = 0;
  ////Long N_noi = 1;
  if(bfac==1)in.rotate_bfac = false;
  in.bfac_write = bfac;
  if(in.rotate_bfac)in.bfac_write = 1;
  in.read = read;
  in.single_file = single_file;
  in.filename = std::string(filename);

  ////size_t off_file = 0;
  ////inputpara in;
  if(in.read == true){
    in.load_para(filename, false);
    if(in.VECS_TYPE != VECS_TYPE){
      qmessage("Noise type wrong, file %s, input %s, expect %s \n", filename, VECS_TYPE.c_str(), in.VECS_TYPE.c_str());
      abort_r("");
    }
    if(in.nvec <= 0){qmessage("%s \n", filename);abort_r("File noise vector size Wrong! \n");}
    if(in.bfac != 1){Qassert(in.bfac == bfac);}
    in.bfac_write = in.bfac;

    if(in.OBJECT != std::string("BEGIN_Vecs_HEAD")){abort_r("File head wrong");}
    Int type = get_save_type(in.save_type);
    if(type == 0){in.bsize=sizeof(RealD);in.single_file=false;}
    if(type == 1){in.bsize=sizeof(float) ;in.single_file=true; }

    //////Check file sizes
    size_t sizen = get_file_size_MPI(filename) - in.off_file;  //Qassert(sizen == string_to_size(in.total_size));
    if(sizen != string_to_size(in.total_size)){
      qmessage("size  %zu %zu .", sizen, string_to_size(in.total_size));
      abort_r("FILE size not match with head !\n");}

    size_t Vsize = size_t(in.nx)*in.ny*in.nz*in.nt*size_t(bfac*2);
    in.N_noi = in.nvec/(bfac/in.bfac_write);
    size_t Fsize = (in.N_noi + 0)*Vsize*in.bsize;  //Qassert(Fsize <= string_to_size(in.total_size));
    if(Fsize > string_to_size(in.total_size)){abort_r("FILE size too small for vectors read !\n");}
    if(Fsize != string_to_size(in.total_size)){in.do_checksum = false;}

    //////Check file sizes
    ////off_file = in.off_file + n0*Vsize*bsize;
    in.end_of_file = sizen + in.off_file;

  }
  if(in.read == false){
    in.VECS_TYPE = VECS_TYPE;
    in.INFO_LIST  = INFO_LIST;
    if(!IO_ENDIAN){in.FILE_ENDIAN = std::string("BIGENDIAN");}
    if( IO_ENDIAN){in.FILE_ENDIAN = std::string("LITTLEENDIAN");}
    if(in.N_noi <= 0){qmessage("write noises size zero \n");return;}
    in.nvec = in.N_noi*(bfac/in.bfac_write);

    //////in.nx = io_use.nx;in.ny = io_use.ny;in.nz = io_use.nz;in.nt = io_use.nt;
    if(in.nx == 0 or in.ny == 0 or in.nz == 0 or in.nt == 0){abort_r("Set up input dim first to write!\n");}
    in.bfac = in.bfac_write;in.checksum = 0;
    if(in.single_file==false){in.bsize=sizeof(RealD);in.save_type = std::string("Double");}
    if(in.single_file==true ){in.bsize=sizeof(float) ;in.save_type = std::string("Single");}

    size_t Fsize = size_t(in.N_noi) * in.nx*in.ny*in.nz*in.nt* size_t(bfac*2);
    Fsize = Fsize*in.bsize;in.total_size = print_size(Fsize);
    ////qmessage("size of file %zu \n", Fsize);
    vecs_head_write(in, filename, true);
    in.end_of_file = Fsize + in.off_file;

  }

  if(in.bfac_write == bfac){in.rotate_bfac = false;}
}

inline void close_file_qlat_noisesT(FILE* file, io_vec& io_use, inputpara& in)
{
  io_use.io_close(file);
  if(in.read==false){
    in.checksum = io_use.full_crc;
    vecs_head_write(in, in.filename.c_str(), false);
  }
  if(in.read==true and in.ncur == in.N_noi and io_use.do_checksum == true){
    if(in.checksum != io_use.full_crc){
      qmessage("File %s check sum wrong, %12X %12X ! \n ", in.filename.c_str(), in.checksum, io_use.full_crc);
      // abort only when input checksum is not zero
      Int bad_crc = io_use.full_crc;
      if(bad_crc != 0)
      {
        abort_r("");
      }    
    }    
  }
}

template <class Ty>
void load_qlat_noisesT_core(FILE* file, std::vector<Ty*  > &noises, const Int bfac, const QMEM Gmem, const Geometry& geo, io_vec& io_use, inputpara& in, Int n0=0, Int n1=-1)
{
  bool read        = in.read;
  bool rotate_bfac = in.rotate_bfac;
  Int bfac_write   = in.bfac_write;
  size_t bsize     = in.bsize;
  bool single_file = in.single_file;

  if(read == true ){if(n1 == -1){n1 = in.N_noi;}}
  if(read == false){
    if(n1 == -1){n1 = n0 + noises.size();}
    if(n1 != -1){
    if(n1 - n0 > int(noises.size())){
      qmessage("Give more noises %d, n0 %d, n1 %d !\n", int(noises.size()), n0, n1);
      abort_r();}
    }
  }
  if(n1 <= n0 or n0 < 0){qmessage("Need read more vectors, n0 %d, n1 %d !\n", n0, n1);abort_r();}
  if(n1 > in.N_noi){qmessage("Need set input more vectors, n1 %d, N_noi %d !\n", n1, in.N_noi);abort_r();}
  in.ncur = n1;
  Int nread = n1 - n0;

  if(read == true){
    Qassert(noises.size() == (LInt) nread);
  }

  /////if(read == false){geo = noises[0].geo();}
  size_t Vsize = size_t(in.nx)*in.ny*in.nz*in.nt*size_t(bfac*2);
  //size_t off_file = in.off_file + n0*Vsize*bsize;
  size_t off_file = size_t(n0)*Vsize*bsize;
  /////qmessage(" ionum off %zu, n0 %zu, n1 %zu, Vsize %zu, bsize %zu \n", off_file, size_t(n0), size_t(n1), Vsize, size_t(bsize));
  io_use.io_off(file, off_file, true);

  Int io_gn = IO_GN;
  if(in.nvec/io_gn < 1){io_gn = 1;}

  bool Rendian = IO_ENDIAN;
  if(in.FILE_ENDIAN == std::string("BIGENDIAN")){   Rendian = false;}
  if(in.FILE_ENDIAN == std::string("LITTLEENDIAN")){Rendian = true ;}

  /////int bufN = io_use.ionum;
  Int bufN = io_use.ionum;
  /////qmessage("ionum %d %d \n", io_use.ionum, N_noi);

  //void* buf;
  //buf = aligned_alloc_no_acc(bufN*bfac*io_use.noden * 2*bsize);
  qlat::vector<char > buf_vec;buf_vec.resize(bufN*bfac*io_use.noden * 2*bsize);
  void* buf = (void*) buf_vec.data();
  qlat::ComplexD*  bufD = (qlat::ComplexD* ) buf;
  qlat::ComplexF* bufF = (qlat::ComplexF*) buf;

  move_index mv_civ;

  //////false big endian, true small endian
  Int bi = 0;
  for(Int ni = 0; ni < nread; ni++)
  {
    if(read==false){
      if(!single_file)copy_noise_to_vec(noises[ni], &bufD[bi*bfac*io_use.noden], geo, bfac, 1, Gmem);
      if( single_file)copy_noise_to_vec(noises[ni], &bufF[bi*bfac*io_use.noden], geo, bfac, 1, Gmem);}
    bi = bi + 1;

    if(bi == bufN or ni == (nread - 1)){
    if(read==false)if(rotate_bfac)mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 1, 2*bsize, false);

    read_kentucky_vector(file,(char*) buf, bi*bfac/bfac_write, io_use, Rendian, bfac_write*bsize*2, single_file, io_gn , read);

    if(read==true)if(rotate_bfac)mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 0, 2*bsize, false);
    if(read==true)for(Int nbi=0; nbi < bi; nbi++){int na = ni - bi + 1;
    {
      if(!single_file)copy_noise_to_vec(noises[na + nbi], &bufD[nbi*bfac*io_use.noden], geo, bfac, 0, Gmem);
      if( single_file)copy_noise_to_vec(noises[na + nbi], &bufF[nbi*bfac*io_use.noden], geo, bfac, 0, Gmem);
    }}

    bi = 0;
    }
  }

  buf = NULL;bufD = NULL; bufF = NULL;
  //free(buf);

}

template <class Ty>
void initialize_Field(std::vector<qlat::FieldG<Ty> > &noises, const Int bfac, const Int nread, inputpara& in)
{
  Qassert(in.nx != 0 and in.ny != 0 and in.nz != 0 and in.nt != 0);
  const Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
  const Geometry& geo = get_geo_cache(total_site);
  if(noises.size() != (LInt) nread){
    noises.resize(0);
    noises.resize(nread);
    for(LInt ni=0;ni<noises.size();ni++){
      noises[ni].init(geo, bfac, QMGPU, QLAT_DEFAULT);
    }
  }
}

/*
  initialize noises with in, if in.read
  setup bufP from noises
  return n1 --> end of files or noises size
*/
template <class Ty, Int bfac>
inline Int load_qlat_noisesT_ini(std::vector<Ty* >& bufP, std::vector<qlat::FieldM<Ty, bfac> > &noises, const Int bc, inputpara& in, Int n0, Int n1=-1)
{
  Qassert(bc == bfac);
  if(in.read == true){
    if(n1 == -1){n1 = in.N_noi;}
  }
  if(in.read == false){
    if(n1 == -1){n1 = n0 + noises.size();}
    if(n1 != -1){
    if(n1 - n0 > int(noises.size())){
      qmessage("Give more noises %d, n0 %d, n1 %d !\n", int(noises.size()), n0, n1);
      abort_r();}
    }
  }
  if(n1 <= n0 or n0 < 0){qmessage("Need read more vectors, n0 %d, n1 %d !\n", n0, n1);abort_r();}
  if(n1 > in.N_noi){qmessage("Need set input more vectors, n1 %d, N_noi %d !\n", n1, in.N_noi);abort_r();}
  Int nread = n1 - n0;

  if(in.read == true){
    //Geometry geo;
    //Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
    //geo.init(total_site);
    QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
    const Geometry& geo = get_geo_cache(Coordinate(in.nx, in.ny, in.nz, in.nt));
    QLAT_DIAGNOSTIC_POP;

    if(noises.size() != (LInt) nread){
      noises.resize(0);
      noises.resize(nread);
      for(LInt ni=0;ni<noises.size();ni++)noises[ni].init(geo);
    }
  }

  bufP.resize(noises.size());
  for(unsigned int iv=0;iv<noises.size();iv++)
  {
    bufP[iv] = (Ty*) qlat::get_data(noises[iv]).data();
  }

  return n1;
}

/*
  Copy of the same function for FieldG for initialization
  May try to initialize fields on seperate functions ?
*/
template <class Ty>
inline Int load_qlat_noisesT_ini(std::vector<Ty* >& bufP, std::vector<qlat::FieldG<Ty> > &noises, const Int bfac, inputpara& in, Int n0, Int n1=-1)
{
  if(in.read == true){
    if(n1 == -1){n1 = in.N_noi;}
  }
  if(in.read == false){
    if(n1 == -1){n1 = n0 + noises.size();}
    if(n1 != -1){
    if(n1 - n0 > int(noises.size())){
      qmessage("Give more noises %d, n0 %d, n1 %d !\n", int(noises.size()), n0, n1);
      abort_r();}
    }
  }
  if(n1 <= n0 or n0 < 0){qmessage("Need read more vectors, n0 %d, n1 %d !\n", n0, n1);abort_r();}
  if(n1 > in.N_noi){qmessage("Need set input more vectors, n1 %d, N_noi %d !\n", n1, in.N_noi);abort_r();}
  Int nread = n1 - n0;

  if(in.read == true){
    //Geometry geo;
    //Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
    //geo.init(total_site);
    const Geometry& geo = get_geo_cache(Coordinate(in.nx, in.ny, in.nz, in.nt));

    if(noises.size() != (LInt) nread){
      noises.resize(0);
      noises.resize(nread);
      for(LInt ni=0;ni<noises.size();ni++){
        noises[ni].init(geo, bfac, QMGPU, QLAT_DEFAULT);
      }
    }
  }

  // dc inside vol
  if(bfac != 1){Qassert(noises[0].mem_order == QLAT_DEFAULT);}
  bufP.resize(noises.size());
  for(unsigned int iv=0;iv<noises.size();iv++)
  {
    bufP[iv] = (Ty*) qlat::get_data(noises[iv]).data();
  }

  return n1;
}


/*
  inputpara in should be initialized
*/
template <class Ty, Int bfac>
void load_qlat_noisesT(FILE* file, std::vector<qlat::FieldM<Ty, bfac> > &noises, io_vec& io_use, inputpara& in, Int n0=0, Int n1=-1)
{
  std::vector<Ty* > bufP;
  n1 = load_qlat_noisesT_ini(bufP, noises, bfac, in, n0, n1);
  load_qlat_noisesT_core(file, bufP, bfac, QMCPU, noises[0].geo(), io_use, in, n0, n1);
}

////initialize the instruct and end of file
template <class Ty>
Geometry load_qlat_noisesT_file_ini(const char *filename, const Int N_noi, const Int bfac, inputpara& in, const Geometry& geo, bool read=true, bool single_file=true, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), const bool rotate_bfac = true){
  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ) and IsTypeComplex<Ty>() == 0){
    abort_r("Cannot understand the input format! \n");}

  open_file_qlat_noisesT(filename, bfac, in, read, single_file, N_noi, VECS_TYPE, INFO_LIST, rotate_bfac);

  Geometry geo_copy = geo;
  if(read == true ){
    Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
    geo_copy.init(total_site);
  }

  IOvecKey fkey(geo_copy, IO_DEFAULT, in.do_checksum);
  io_vec& io_use = get_io_vec_plan(fkey);
  io_use.end_of_file = in.end_of_file;
  return geo_copy;
}

template <class Ty>
void load_qlat_noisesT(const char *filename, std::vector<Ty*  > &noises, const Int bfac, inputpara& in, const Geometry& geo, bool read=true, bool single_file=true, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), Int n0=0,Int n1=-1, bool rotate_bfac = true){
  TIMERB("load_qlat_noisesT kernel");

  Long N_noi = 0;
  if(read == false){
    ////if(n0 != 0 or n1 != -1){abort_r("Write mode shoude have n0 0, n1 -1 . ! \n");}
    N_noi = noises.size();
  }
  Geometry geo_copy = load_qlat_noisesT_file_ini<Ty>(filename, N_noi, bfac, in, geo, read, single_file, VECS_TYPE, INFO_LIST, rotate_bfac);

  if(in.read == true){
    if(n1 == -1){n1 = in.N_noi;}
  }
  if(in.read == false){
    if(n1 == -1){n1 = n0 + noises.size();}
  }

  //////io_vec io_use(geo, IO_DEFAULT, IO_THREAD, in.do_checksum);
  IOvecKey fkey(geo_copy, IO_DEFAULT, in.do_checksum);
  io_vec& io_use = get_io_vec_plan(fkey);
  io_use.end_of_file = in.end_of_file;

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");
  //size_t off_file = in.off_file + n0*Vsize*bsize;
  /////qmessage(" ionum off %zu, n0 %zu, n1 %zu, Vsize %zu, bsize %zu \n", off_file, size_t(n0), size_t(n1), Vsize, size_t(bsize));

  io_use.io_off(file, in.off_file, true);  ////shift file for the head
  load_qlat_noisesT_core(file, noises, bfac, QMCPU, geo_copy, io_use, in, n0, n1);

  close_file_qlat_noisesT(file, io_use, in);
}

template <class Ty, Int bfac>
void load_qlat_noisesT(const char *filename, std::vector<qlat::FieldM<Ty, bfac> > &noises, bool read=true, bool single_file=true, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), Int n0=0,Int n1=-1, bool rotate_bfac = true){

  Long N_noi = 0;
  inputpara in; Geometry geo;
  if(read == false){
    Qassert(noises[0].initialized);
    N_noi = noises.size();
    geo = noises[0].geo();
    in.read_geo(geo);
  }

  geo = load_qlat_noisesT_file_ini<Ty>(filename, N_noi, bfac, in, geo, read, single_file, VECS_TYPE, INFO_LIST, rotate_bfac);

  std::vector<Ty* > bufP;n1 = load_qlat_noisesT_ini(bufP, noises, bfac, in, n0, n1);

  IOvecKey fkey(geo, IO_DEFAULT, in.do_checksum);
  io_vec& io_use = get_io_vec_plan(fkey);
  io_use.end_of_file = in.end_of_file;

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");

  io_use.io_off(file, in.off_file, true);  ////shift file for the head
  load_qlat_noisesT_core(file, bufP, bfac, QMCPU, geo, io_use, in, n0, n1);

  close_file_qlat_noisesT(file, io_use, in);
}

/*
  need bfac
  QMGPU by default
*/
template <class Ty>
void load_qlat_noisesG(const char *filename, std::vector<qlat::FieldG<Ty> > &noises, const Int bfac, bool read=true, bool single_file=true, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), Int n0=0,Int n1=-1, bool rotate_bfac = true){

  Long N_noi = 0;
  inputpara in; Geometry geo;
  if(read == false){
    Qassert(noises[0].initialized);
    N_noi = noises.size();
    geo = noises[0].geo();
    in.read_geo(geo);
  }

  geo = load_qlat_noisesT_file_ini<Ty>(filename, N_noi, bfac, in, geo, read, single_file, VECS_TYPE, INFO_LIST, rotate_bfac);

  std::vector<Ty* > bufP;n1 = load_qlat_noisesT_ini(bufP, noises, bfac, in, n0, n1);

  IOvecKey fkey(geo, IO_DEFAULT, in.do_checksum);
  io_vec& io_use = get_io_vec_plan(fkey);
  io_use.end_of_file = in.end_of_file;

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");

  io_use.io_off(file, in.off_file, true);  ////shift file for the head
  load_qlat_noisesT_core(file, bufP, bfac, QMGPU, geo, io_use, in, n0, n1);

  close_file_qlat_noisesT(file, io_use, in);
}

template <class T, Int civ>
void load_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, civ> > &noises, bool read , bool single_file=true, const std::string& INFO_LIST = std::string("NONE"),Int n0 = 0, Int n1=-1, std::string info = "NONE")
{
  TIMERC("load/save qlat eigen");
  std::string VECS_TYPE;
  if( info == std::string("NONE") ){
    std::string tmp = ssprintf("Eigen_system_nvec.%d.tzyx.R/I", civ);
    VECS_TYPE = tmp;
  }else{VECS_TYPE = info;}
  
  //VECS_TYPE = std::string("Eigen_system_nvec.12.tzyx.R/I");
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST, n0, n1);
}

template <class T, Int civ>
void load_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, civ> > &noises, Int n0=0,Int n1=-1, std::string info = "NONE")
{
  bool read = true; bool single_file = true; std::string INFO_LIST = std::string("NONE");
  load_qlat_eigen(filename, noises, read, single_file, INFO_LIST, n0, n1, info);
}

template <class T, Int civ>
void save_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, civ> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), std::string info = "NONE"){
  load_qlat_eigen(filename, noises, false, single_file, INFO_LIST, 0, -1, info);
}

////===Check eigen system type
inline Int check_eigen_qlat(const char *filename, Int n1, inputpara& in)
{
  in.load_para(filename, false);
  if(in.VECS_TYPE != std::string("NONE"))
  {
    if(in.OBJECT != std::string("BEGIN_Vecs_HEAD")){
      qmessage("%s \n", filename);
      abort_r("File head wrong");
    }

    Qassert(in.nvec > 0);////Qassert(in.bfac == bfac_write);
    Int type = get_save_type(in.save_type);
    in.bsize = 8; in.single_file = true;
    if(type == 0){in.bsize=sizeof(RealD);in.single_file=false;}
    if(type == 1){in.bsize=sizeof(float) ;in.single_file=true; }
    //////Check file sizes
    size_t sizen = get_file_size_MPI(filename) - in.off_file;  //Qassert(sizen == string_to_size(in.total_size));
    if(sizen != string_to_size(in.total_size)){abort_r("FILE size not match with head !\n");}

    in.Vsize = size_t(in.nx)*in.ny*in.nz*size_t(in.bfac*2);
    size_t Fsize = size_t(n1)*(12)*in.Vsize*in.bsize;  //Qassert(Fsize <= string_to_size(in.total_size));
    if(Fsize  > string_to_size(in.total_size)){abort_r("FILE size too small for vectors read !\n");}

    ////string_to_size(in.total_size) + in.off_file;
    /////io_use.end_of_file = sizen + in.off_file;

    //////Check file sizes
    ////off_file = in.off_file + nini*Vsize*bsize;
    /////if(Fsize != string_to_size(in.total_size)){io_use.do_checksum = false;}

    if(!in.single_file){return 2;}
    if( in.single_file){return 3;}

    /////return 2;  ////qlat noises type
  }
  return  -1;
}

////return 0, gwu double, 1 gwu single, 2 qlat eigen with inputpar written
inline Int check_qlat_eigen_file_type(const char *filename, io_vec &io_use,Int n1, inputpara& in){
  Int type = check_eigen_qlat(filename, n1, in);
  if(type != -1){return type;}

  bool check = true;
  type = check_Eigen_file_type(filename, io_use, n1, check);
  return type;
}
////===Check eigen system type


template <class T, Int bfac>
void load_qlat_noises(const char *filename, std::vector<qlat::FieldM<T, bfac> > &noises, bool read=true, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), Int n0=0,Int n1=-1)
{
  TIMERB("load/save qlat noises");
  std::string VECS_TYPE = std::string("Noise_Vectors");
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST, n0, n1);
}

template <class T, Int bfac>
void save_qlat_noises(const char *filename, std::vector<qlat::FieldM<T, bfac> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_noises(filename, noises, false, single_file, INFO_LIST);
}

template <class T, Int bfac>
void save_qlat_noises(const std::string& filename, std::vector<qlat::FieldM<T, bfac> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_noises(filename.c_str(), noises, false, single_file, INFO_LIST);
}

template <class T, Int bfac>
void load_qlat_noise(const char *filename, qlat::FieldM<T, bfac> &noise, bool read=true, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  std::vector<qlat::FieldM<T, bfac> > noises;
  ///noises[0].init(noise.geo());
  if(read == false){noises.resize(1);noises[0] = noise;}
  load_qlat_noises(filename, noises, read, single_file, INFO_LIST, 0, 1);
  if(read == true ){noise = noises[0];}
}

template <class T, Int bfac>
void save_qlat_noise(const char *filename, qlat::FieldM<T, bfac> &noise, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_noise(filename, noise, false, single_file, INFO_LIST);
}

/* 
  Templates for FieldG loaders
  Need sort the code for more general purpose
  FieldG default bfac = 1
*/
template <class T>
void load_qlat_noises(const char *filename, std::vector<qlat::FieldG<T> > &noises, bool read=true, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), Int n0=0,Int n1=-1, const Int bfac = 1)
{
  TIMERB("load/save qlat noises");
  std::string VECS_TYPE = std::string("Noise_Vectors");
  load_qlat_noisesG(filename, noises, bfac, read, single_file, VECS_TYPE, INFO_LIST, n0, n1);
}

template <class T>
void save_qlat_noises(const char *filename, std::vector<qlat::FieldG<T> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), const Int bfac = 1){
  load_qlat_noises(filename, noises, false, single_file, INFO_LIST, 0, -1, bfac);
}

template <class T>
void save_qlat_noises(const std::string& filename, std::vector<qlat::FieldG<T> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), const Int bfac = 1){
  load_qlat_noises(filename.c_str(), noises, false, single_file, INFO_LIST, 0, -1, bfac);
}


//////Assume memory allocated already
template<class T, typename Td>
void copy_noise_to_prop(qlat::FieldM<T, 12*12>& noise, Propagator4dT<Td>& prop, Int dir=1)
{
  TIMERB("copy_noise_to_prop");
  if(dir == 1){prop.init(noise.geo());}
  if(dir == 0){noise.init(prop.geo());}
  T* noi = (T*) qlat::get_data(noise).data();
  #pragma omp parallel for 
  for (Long index = 0; index < prop.geo().local_volume(); ++index)
  {
    ///qlat::WilsonMatrixT<T>& src =  prop.get_elem_offset(index);
    T* src   = (T*) &noi[index*12*12];
    qlat::ComplexT<Td>* res  = &prop.get_elem_offset(index)(0,0);
    
    for(Int d0=0;d0<12;d0++)
    {
      for(Int d1=0;d1<12;d1++)
      {
        //////copy to prop
        if(dir==0){src[d1*12+d0] = res[d0*12+d1];}
        //////copy to buf
        if(dir==1){res[d0*12+d1] = src[d1*12+d0];}

      }

    }
  }
}

//////Assume memory allocated already
template<class T, typename Td>
void copy_noises_to_prop(std::vector<qlat::FieldM<T, 12*12> >& noises, Propagator4dT<Td>& prop, Int dir=1)
{
  TIMERB("copy_noises_to_prop");
  if(dir == 1){
    if(!prop.initialized or prop.geo() != noises[0].geo())
    {
      prop.init(noises[0].geo());
    }
  }
  if(dir == 0){
    if(noises.size() != 1 or !noises[0].initialized or noises[0].geo() != prop.geo())
    {
      noises.resize(0);noises.resize(1);noises[0].init(prop.geo());
    }
  }
  copy_noise_to_prop(noises[0], prop, dir);
}

template <typename Td>
void load_qlat_prop(const char *filename, Propagator4dT<Td>& prop, bool read=true, bool single_file=true){
  std::string VECS_TYPE = std::string("Propagator");
  std::string INFO_LIST  = std::string("src 12, sink 12, tzyx, R/I");
  std::vector<qlat::FieldM<qlat::ComplexT<Td >, 12*12> > noises;
  if(read == false){copy_noises_to_prop(noises, prop, 0);}
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST);
  if(read == true){copy_noises_to_prop(noises, prop, 1);}
  
}

template <typename Td>
void load_qlat_prop(const std::string& filename, Propagator4dT<Td>& prop, bool read=true, bool single_file=true){
  load_qlat_prop(filename.c_str(), prop, read, single_file);
}

template <typename Td >
void save_qlat_prop(const char *filename,Propagator4dT<Td >& prop, bool single_file=true){
  load_qlat_prop(filename, prop, false, single_file);
}

template <typename Td >
void save_qlat_prop(const std::string& filename,Propagator4dT<Td >& prop, bool single_file=true){
  load_qlat_prop(filename.c_str(), prop, false, single_file);
}

// need some cleaning .........
template <typename Td>
void load_qlat_prop_buf(const char *filename, Propagator4dT<Td>& prop, std::vector<qlat::FieldM<qlat::ComplexT<Td >, 12*12> >& noises, bool read=true, bool single_file=true){
  std::string VECS_TYPE = std::string("Propagator");
  std::string INFO_LIST  = std::string("src 12, sink 12, tzyx, R/I");
  if(read == false){copy_noises_to_prop(noises, prop, 0);}
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST);
  if(read == true){copy_noises_to_prop(noises, prop, 1);}
}

template <typename Td >
void save_qlat_prop_buf(const char *filename,Propagator4dT<Td >& prop, std::vector<qlat::FieldM<qlat::ComplexT<Td >, 12*12> >& noises, bool single_file=true){
  load_qlat_prop_buf(filename, prop, noises, false, single_file);
}

template <typename Td>
void load_qlat_link(const char *filename,GaugeFieldT<Td> &gf, bool read = true , bool single_file=false){
  std::string VECS_TYPE = std::string("Links");
  std::string INFO_LIST  = std::string("dir 4, cxc 9, tzyx, R/I");
  Qassert(gf.initialized);
  qlat::Geometry& geo = gf.geo();

  inputpara in;
  if(read == false){
    in.read_geo(geo);
  }

  std::vector<qlat::ComplexT<Td >* > src; src.resize(1);
  src[0] = (qlat::ComplexT<Td >*) qlat::get_data(gf).data();

  load_qlat_noisesT<qlat::ComplexT<Td >>(filename, src, 4*9, in, geo, read, single_file, VECS_TYPE, INFO_LIST);
}

template <typename Td>
void load_qlat_link(const std::string& filename,GaugeFieldT<Td> &gf, bool read = true , bool single_file=false){
  load_qlat_link(filename.c_str(), gf, read, single_file);
}

template <typename Td>
void save_qlat_link(const char *filename,GaugeFieldT<Td> &gf, bool single_file = false){
  load_qlat_link(filename, gf, false, single_file);
}

template <typename Td>
void save_qlat_link(const std::string& filename,GaugeFieldT<Td> &gf, bool single_file = false){
  load_qlat_link(filename.c_str(), gf, false, single_file);
}

/////nvec needed for checksum
inline FILE* open_eigensystem_file(const char *filename, Int nini, Int nvec, bool read, io_vec& io_use, inputpara& in, Int save_type = 3)
{
  in.single_file = true;
  in.read = read;
  in.filename = std::string(filename);

  if(read==true ){in.file_type = check_qlat_eigen_file_type(filename, io_use, nvec, in);}
  if(read == false){
    if(save_type != 0 and save_type != 1 and save_type != 2 and save_type != 3){abort_r("Save type unknown! \n ");}
    in.file_type = save_type;
  }

  if(in.file_type == 0 or in.file_type == 2){in.single_file = false;}
  if(in.file_type == 1 or in.file_type == 3){in.single_file = true ;}

  if(in.file_type == 0 or in.file_type == 1){
    io_use.do_checksum = false;
  }
  if(in.file_type == 2 or in.file_type == 3){
    if(read == false){in.read_geo(io_use.geo());}
    Int bfac_eigen = 12;
    std::string VECS_TYPE = std::string("Eigen_system_nvec.12.tzyx.R/I");
    std::string INFO_LIST = std::string("NONE");
    bool rotate_bfac = true;
    open_file_qlat_noisesT(filename, bfac_eigen, in, read, in.single_file, nvec, VECS_TYPE, INFO_LIST, rotate_bfac);

    if(io_use.end_of_file != 0){abort_r("Please close file first! \n");}

    io_use.end_of_file = in.end_of_file;
    io_use.do_checksum = in.do_checksum;
  }

  //qmessage("test sites %d %d %d %d ; %d %d %d %d \n", io_use.nx, io_use.ny, io_use.nz, io_use.nt, in.nx, in.ny, in.nz, in.nt);
  if(in.nx != 0){Qassert(io_use.nx == in.nx);}
  if(in.ny != 0){Qassert(io_use.ny == in.ny);}
  if(in.nz != 0){Qassert(io_use.nz == in.nz);}
  if(in.nt != 0){Qassert(io_use.nt == in.nt);}

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");

  //////shift the file to nini position
  {
    size_t bsize = sizeof(RealD);int bfac = 12;
    if( in.single_file){bsize=sizeof(float) ;}
    if(!in.single_file){bsize=sizeof(RealD);}
    size_t Vsize = size_t(io_use.nx)*io_use.ny*io_use.nz*io_use.nt*size_t(bfac*2);
    size_t off_file = size_t(nini)*Vsize*bsize;
    ////qmessage("off each %ld, nini %ld %ld %ld...\n", Long(off_file), Long(nini), Long(Vsize), Long(bsize));
    if(in.file_type == 2 or in.file_type == 3){off_file += in.off_file;}
    io_use.io_off(file, off_file, false);
  }

  return file;
}


template <class Ty>
void load_eigensystem_vecs(FILE* file, std::vector<qlat::FieldM<Ty, 12> > &noises, io_vec& io_use, inputpara& in, Int n0=0, Int n1=-1)
{

  if(in.file_type == 0 or in.file_type == 1)
  {
    if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
    Int n_vec = n1-n0;
    ////if(in.file_type == 0){read_single = false;}
    ////if(in.file_type == 1){read_single = true ;}

    bool check = true;
    ////load_gwu_eigen(file, noises, io_use, n0, n1, check, in.read, in.single_file );
    if(in.read == true){
    if(noises.size() < (LInt) n_vec)
    {
      noises.resize(0);
      noises.resize(n_vec);
      for(Int iv=0;iv<n_vec;iv++){noises[iv].init(io_use.geo());}
    }}

    if(sizeof(Ty) == 2*sizeof(RealD)){
      std::vector<double* > respD;respD.resize(n_vec);
      for(Int iv=0;iv<n_vec;iv++){respD[iv] = (double*) qlat::get_data(noises[iv]).data();}
      load_gwu_eigen(file, respD, io_use,n0,n1,check, in.read, in.single_file );
    }
    if(sizeof(Ty) == 2*sizeof(float) ){
      std::vector<RealF*  > respF;respF.resize(n_vec);
      for(Int iv=0;iv<n_vec;iv++){respF[iv] = (RealF*) qlat::get_data(noises[iv]).data();}
      load_gwu_eigen(file, respF, io_use,n0,n1,check, in.read, in.single_file );
    }

  }

  if(in.file_type == 2 or in.file_type == 3)
  {
    load_qlat_noisesT(file, noises, io_use, in, n0, n1);
  }

}

inline void close_eigensystem_file(FILE* file, io_vec& io_use, inputpara& in){

  if(in.file_type == 0 or in.file_type == 1){io_use.io_close(file);}

  if(in.file_type == 2 or in.file_type == 3)
  {
    close_file_qlat_noisesT(file, io_use, in);
  }

}

template <class Td>
void load_qlat_vecs(const char *filename, Td* prop, const Int nvec, io_vec &io_use, const bool single = false, bool Rendian=false)
{
  ////if(sizeof(Td) != sizeof(double ) and sizeof(Td) != sizeof(float )){abort_r("Cannot understand the input format! \n");}
  size_t noden = io_use.noden;
  size_t Fsize = (noden)*io_use.Nmpi*nvec*sizeof(float);

  //bool read = true;
  if(single == false){Fsize = Fsize * 2;}
  size_t sizen = get_file_size_MPI(filename);

  if(sizen != 2*Fsize and sizen != Fsize){qmessage("File %s \n",filename);abort_r("File size wrong! \n");}
  const Int ncomplex = sizen / Fsize;
  Int s_inner = sizeof(float);if(single == false){s_inner = sizeof(RealD);}
  s_inner = s_inner * ncomplex;

  FILE* file = io_use.io_read(filename,"rb");
  read_kentucky_vector(file, (char*) prop, nvec, io_use, Rendian, s_inner, single, 1);
  io_use.io_close(file);
}


//template<typename Ty, Int dir>
//void copy_eo_cs_to_fieldM(qlat::vector_gpu<Ty >& res, const Int civ, const Geometry& geo, vector_cs<Ty >& even, vector_cs<Ty >& odd,
//  Int e0, Int e1, Int o0, Int o1, qlat::vector<Long >& map, Int mode = 0)
//{
//  TIMER_FLOPS("copy_eo_cs_to_fieldM");
//  Qassert(even.initialized and odd.initialized);
//  ////fft_desc_basic& fd = get_fft_desc_basic_plan(res.geo());
//  const bool GPU = even.GPU;
//  const Long V = geo.local_volume();
//  const Long Vh= V/2;
//  if(dir == 1){if(Long(res.size()) != civ*V){res.resize(civ * V, GPU);}}
//  if(dir == 0){Qassert(Long(res.size()) == V * civ );}
//  Qassert(res.GPU == even.GPU);
//
//  const Int DIM = 3;
//  Qassert(civ % DIM == 0);
//  if(map.size() == 0){get_index_mappings_reverse(map, geo);}
//  Int nvec = civ / DIM;
//  Qassert(e1 > e0 and o1 > o0);
//  Int ne = e1 - e0;
//  Int no = o1 - o0;
//  Qassert(ne == no);
//  Qassert(ne <= nvec and no <= nvec and e1 <= even.nvec and o1 <= odd.nvec);
//  const Int b_size = even.b_size;
//  const Int btotal = even.btotal;
//
//  Qassert(btotal * b_size == DIM*V/2);
//  qlat::vector<Ty** > eP;eP.resize(2*ne);
//  //qlat::vector<Ty** > oP;oP.resize(no);
//  ////Ty** eP = even.get_pointers(ni)
//  for(Int ei=0;ei<ne;ei++){eP[ei]      = even.get_pointers(ei + e0);}
//  for(Int oi=0;oi<no;oi++){eP[ne + oi] =  odd.get_pointers(oi + o0);}
//  const Long* mapP = (Long*) qlat::get_data(map).data();
//
//  ////NtNzNyNx, DIM x nvec
//  Ty* r = (Ty*) qlat::get_data(res).data();
//  //if(mode == 1)
//  {
//  qGPU_for(qi, V/2, GPU,{
//    //const Int ni = ci / DIM;
//    //const Int c  = ci % DIM;
//    for(Int eo = 0; eo < 2;eo++)
//    {
//      const Long quda_idx = eo*Vh + qi;
//      const Long qlat_idx_4d = mapP[quda_idx];
//      Ty* rr = &r[qlat_idx_4d*civ];
//      for(Int c = 0; c < 3 ; c++)
//      {
//        Long bv =  qi*DIM + c;
//        ////const Long bv = qi*DIM + c; ////quda vectors this order
//        if(mode == 1){bv = c*Vh + DIM;} ////quda vectors this order
//        const Long bi = bv / b_size;
//        const Long bj = bv % b_size;
//        for(Int ni=0;ni<ne;ni++)
//        {
//          {
//          if(dir == 1){rr[c*nvec + ni] = eP[eo*ne + ni][bi][bj];}
//          if(dir == 0){eP[eo*ne + ni][bi][bj] = rr[c*nvec + ni];}
//          }
//        }
//      }
//    }
//  });}
//
//  timer.flops += double(V) * DIM * ne * sizeof(Ty);
//}

/////load the even vectors with zero mass
/////mode_c = 0, default even with vol -> color
template<typename Ty >
inline void load_eo_evecs(const char* filename, vector_cs<Ty >& even, qlat::vector<Ty >& evals, std::vector<double>& err,
  const Geometry& geo, const Int N0=0, const Int N1=-1,
  double mass = 0.0, Int mode_c = 0, const bool single_file = true, const bool read = true ,
  std::string VECS_TYPE = std::string("EO_Eigensystem"), const Int n_off_file = 0)
{
  TIMERB("load_eo_evecs");
  const fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
  if(even.nvec == 0){return ;}
  io_vec& io_use = get_io_vec_plan_with_checksum(geo);
  io_use.do_checksum = true; // turn on checksum and initialize files
  io_use.ini_crc(true);
  const Int DIM = 3;
  const bool rotate_bfac = true; ////default rotate color to outside to save memory

  //char fileE[600];
  std::string fileE = ssprintf("%s.evals", filename);
  Int nini = N0;
  Int Nmax = even.nvec - nini;
  if(N1 != -1){
    Nmax = N1 - nini;
  }
  Qassert(Nmax > 0);
  /////qmessage("save Nmax %3d \n", Nmax);

  Int nvec  = 0;
  Int nhalf = 0;
  double mass_file = 0.0;
  if(read == true){
    inputpara in;
    in.load_para(filename, false);
    nhalf = in.nvec / DIM;Qassert(rotate_bfac == true );
    std::vector<std::string > mL = stringtolist(in.INFO_LIST);
    double mre = stringtodouble(mL[1]);
    Qassert(in.nx == fd.nx and in.ny == fd.ny and in.nz == fd.nz and in.nt == fd.nt);
    ///if(nhalf <= 0){return ;}

    std::vector<double > values, errors;
    ////load only the values exists
    if(get_file_size_MPI(fileE.c_str(), true) > 0){
      load_txt_eigenvalues(values, errors, fileE.c_str());
      if(2 * nhalf >= Long(values.size()) ){nvec = values.size();}////DIM 3 may be needed here for check
      else{nvec = 2* nhalf;}
      Qassert(nvec <= Long(values.size()) );///nhalf*2 could be larger than values.size()
      if(nvec >  Nmax){
        nvec  = Nmax;
        nhalf = (nvec + 1)/2;
      }
      Qassert(nvec == Nmax);

      evals.resize(nvec);err.resize(nvec);
      for(Int n=  0;n<nvec;n++){
        Qassert(qlat::qnorm( values[n*2 + 1] ) < 1e-10);
        evals[n] = Ty(values[n*2+0] - 4.0*mre*mre, 0.0);
        err[n]   = errors[n];
      }
    }else{
      ////half is not correct due to bfac rotations
      nvec = 2 * nhalf;
      if(nvec >  Nmax or nvec == 0){
        nvec  = Nmax;
        nhalf = (nvec + 1)/2;
      }
      ////qmessage("nhalf %5d, nvec %5d, Nmax %5d \n", nhalf, nvec, Nmax);
      Qassert(nvec == Nmax);
    }
  }
  if(read == false){
    nvec = Nmax;
    nhalf = (Nmax + 1)/2;////one more vectors if not devided by 2
    if(nvec <= 0){return ;}

    ////save only it's none zero
    if(evals.size() > 0){
      std::vector<double > values;values.resize(nvec*2);
      mass_file = 0.0;
      for(Int n=  0;n<nvec;n++){
        values[n*2 + 0] = evals[n].real() - 4.0 * mass * mass; /// subtract to zero mass
        values[n*2 + 1] = 0.0;
      }
      save_txt_eigenvalues(values, err, fileE.c_str(), "Fermions EO");
    }
  }
  qmessage("nvec %5d, nhalf %5d, ionum %5d \n", nvec, nhalf, io_use.ionum);
  print_mem_info();

  /////std::string VECS_TYPE("EO_Eigensystem");
  Int Ngroup = io_use.ionum;
  if(Ngroup > nhalf){Ngroup = nhalf;}
  if(Ngroup <= 0){abort_r("ionum wrong!");}
  //char infoL[500];
  std::string infoL = ssprintf("mass %.8f", mass_file);
  std::string INFO_LIST(infoL.c_str());
  const Long V = geo.local_volume();
  qlat::vector<Long > map;
  std::vector<Ty*  > noises;noises.resize(Ngroup);
  std::vector<qlat::vector_gpu<Ty > > eig;eig.resize(Ngroup);
  for(Int iv=0;iv<Ngroup;iv++){
    eig[iv].resize(V * DIM, QMSYNC);///default on SYNC
  }
  for(Int iv=0;iv<Ngroup;iv++){noises[iv] = (Ty*) qlat::get_data(eig[iv]).data();}
  vector_cs<Ty > tmp_end;tmp_end.resize(1, QMSYNC, even);
  qmessage("after io mem allocate");print_mem_info();

  ////load evecs
  {
  inputpara in;
  /////const Int ntotal = nhalf;
  in.nx = io_use.nx;in.ny = io_use.ny;in.nz = io_use.nz;in.nt = io_use.nt;
  ////FILE* file_read  = open_eigensystem_file(filename, nini, ntotal, true , io_use , in_read_eigen , 2);
  ////close_eigensystem_file(file_read , io_use , in_read_eigen );

  Geometry geo_copy = load_qlat_noisesT_file_ini<Ty>(filename, nhalf, DIM, in, geo, read, single_file, VECS_TYPE, INFO_LIST, rotate_bfac);
  Qassert(geo_copy == geo);

  io_use.end_of_file = in.end_of_file;

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");

  io_use.io_off(file, in.off_file, true);  ////shift file for the head

  ////const Int mode_c = 0;
  Int nini_off = 0;////will only affect read == true
  ////only supports offset multiple of 2 due to even-odd writtings together
  if(read == true){Qassert(n_off_file % 2 == 0);}

  std::vector<Long > jobA = job_create(nhalf, Ngroup);
  for(LInt jobi=0;jobi < jobA.size()/2; jobi++)
  {
    if(read == true){if(jobi == 0){nini_off = n_off_file / 2;}}
    if(jobi != 0){nini_off = 0;}///default zero offset, by default it will always read next rounds
    const Long n0   = jobA[jobi*2 + 0];
    const Long ncut = jobA[jobi*2 + 1];
    noises.resize(ncut);for(Int iv=0;iv<ncut;iv++){noises[iv] = (Ty*) qlat::get_data(eig[iv]).data();}
    qmessage("load %5d, dN %5d, N %5d ", int(n0), int(ncut), int(nhalf));
    //qmessage("load %5d, dN %5d, N %5d, inioff %5d ", int(n0), int(ncut), int(nhalf), int(nini_off));

    if(read == false)
      for (Int iv = 0; iv < ncut; iv++) {
        const Int ne = (n0 + iv) * 2 + nini;
        if (ne + 2 <= even.nvec) {
          copy_eo_cs_to_fieldM(eig[iv], 3, geo_copy, even, even, ne, ne + 1,
                               ne + 1, ne + 2, map, mode_c);
        } else {
          copy_eo_cs_to_fieldM(eig[iv], 3, geo_copy, even, tmp_end, ne, ne + 1,
                               0, 1, map, mode_c);
        }
      }

    /////will shift file from current position
    load_qlat_noisesT_core(file, noises, DIM, QMCPU, geo_copy, io_use, in, nini_off, ncut + nini_off);
    ////print_mem_info();

    if(read == true)
    for(Int iv=0;iv<ncut;iv++){
      ////qmessage("iv %8d, ncut %8d \n", iv, int(ncut));
      const Int ne = (n0+iv) * 2 + nini;
      if(ne+2 <= even.nvec){
        copy_fieldM_to_eo_cs(even,   even ,eig[iv], 3, geo_copy, ne, ne+1, ne+1, ne+2, map, mode_c);
      }else{                                               
        copy_fieldM_to_eo_cs(even, tmp_end,eig[iv], 3, geo_copy, ne, ne+1, 0, 1, map, mode_c);
      }
    }
  }

  close_file_qlat_noisesT(file, io_use, in);
  }
  io_use.clear_buf();
}

template<typename Ty >
inline void load_eo_evecs(const std::string& filename, vector_cs<Ty >& even, qlat::vector<Ty >& evals, std::vector<double>& err,
  const Geometry& geo, const Int N0=0, const Int N1=-1,
  double mass = 0.0, Int mode_c = 0, const bool single_file = true, const bool read = true ,
  std::string VECS_TYPE = std::string("EO_Eigensystem"), const Int n_off_file = 0)
{
  load_eo_evecs(filename.c_str(), even, evals, err, geo, N0, N1, mass, mode_c, single_file, read, VECS_TYPE, n_off_file);
}

template<typename Ty>
inline void save_eo_evecs(const char* filename, vector_cs<Ty >& even, qlat::vector<Ty >& evals, std::vector<double>& err,
  const Geometry& geo, const Int N0=0, const Int N1=-1, double mass = 0.0, Int mode_c = 0, const bool single_file = true,
  std::string VECS_TYPE = std::string("EO_Eigensystem"))
{
  load_eo_evecs(filename, even, evals, err, geo, N0, N1, mass, mode_c, single_file, false, VECS_TYPE);
}

template<typename Ty>
inline void save_eo_evecs(const std::string& filename, vector_cs<Ty >& even, qlat::vector<Ty >& evals, std::vector<double>& err,
  const Geometry& geo, const Int N0=0, const Int N1=-1, double mass = 0.0, Int mode_c = 0, const bool single_file = true,
  std::string VECS_TYPE = std::string("EO_Eigensystem"))
{
  load_eo_evecs(filename.c_str(), even, evals, err, geo, N0, N1, mass, mode_c, single_file, false, VECS_TYPE);
}

}

#endif
