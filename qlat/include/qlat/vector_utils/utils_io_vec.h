// utils_io_vec.h
// Gen Wang
// Jan. 2021

#ifndef UTILS_IO_VEC_H
#define UTILS_IO_VEC_H
#pragma once

#include "general_funs.h"
#define IO_DEFAULT  0
#define IO_ENDIAN false
#define IO_GN 2
#define IO_THREAD -1

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

  std::vector<int > node_ioL;
  int ionum;
  int nx,ny,nz,nt;
  int Nx,Ny,Nz,Nt;
  size_t vol,noden;
  int rank,Nmpi;
  std::vector<int > nv,Nv;
  qlat::Geometry geop;
  int threadio;
  std::vector<FILE* > file_omp;
  std::vector<size_t > currsend,currrecv,currspls,currrpls;
  int MPI_size_c;
  char* tmp;
  bool do_checksum;
  std::vector<crc32_t > io_crc;
  crc32_t full_crc;

  size_t end_of_file;
  move_index mv_civ;

  char* buf;size_t size_buf;

  //////only used with eigen system

  //////each noden has the memory order to t,z,y,x
  std::vector<std::vector<LInt > > map_Nitoi;
  ///void get_pos(int i);
  //FILE* io_read(const char *filename,const char* mode);
  //void io_close(FILE *file);

  FILE* io_read(const char *filename,const char* mode){
    int curr_threadio = threadio;
    int Nv = omp_get_max_threads();
    if(curr_threadio < 1 or curr_threadio > Nv){curr_threadio = Nv;}
    if(curr_threadio * ionum > 512 ){curr_threadio = 1;}
    threadio = curr_threadio;

    if(do_checksum)ini_crc(true);
    if(node_ioL[rank]>=0){
      /////Currently can only open one file with openmp
      qassert(file_omp.size() == 0);
      if(threadio!=1)for(int i=0;i<threadio;i++){file_omp.push_back(fopen(filename, mode));}
      return fopen(filename, mode);
    }else{return NULL;}

  }
  
  void io_close(FILE *file){
    if(do_checksum){sum_crc();end_of_file = 0;do_checksum=false;}
    if(node_ioL[rank]>=0){
      if(threadio!=1)for(int i=0;i<threadio;i++){fclose(file_omp[i]);}
      file_omp.resize(0);
      fclose(file);
    }
    file = NULL;
    if(buf != NULL){free(buf);size_buf = 0;buf = NULL;}
  }

  /////default initial file pos
  void io_off(FILE *file,size_t off, bool default_pos = true){
    if( default_pos){
      if(node_ioL[rank]>=0){
      fseek(file , off , SEEK_CUR );}
    }
    if(!default_pos){
      if(node_ioL[rank]>=0){
      fseek(file , off , SEEK_SET );}
    }
  }

  void ini_MPI(int size_c0)
  {
    if(MPI_size_c == size_c0)return;
    MPI_size_c = size_c0;
    size_t size_c = size_c0 * noden;
    if(node_ioL[rank]>=0){
      ////if(tmp != NULL){delete []tmp;tmp=NULL;}
      ////tmp = new char[size_c0 * vol];
      if(tmp != NULL){free(tmp);tmp=NULL;}
      tmp = (char*) aligned_alloc_no_acc(size_c0 * vol);
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
      for(int n=0;n<Nmpi;n++){
        currsend[n] = size_c;
        currspls[n] = size_c*n;
      }
    }
    for(int n=0;n<Nmpi;n++){
      if(node_ioL[n]>=0){
        currrecv[n] = size_c;
        currrpls[n] = size_c*node_ioL[n];
      }
    }
  }

  void ini_crc(bool new_file){
    io_crc.resize(Nmpi);
    for(LInt i=0;i<io_crc.size();i++){io_crc[i] = 0;}
    if(new_file)full_crc = 0;
  }
  crc32_t sum_crc(){
    sum_all_size(&io_crc[0], io_crc.size());
    for(LInt i=0;i<io_crc.size();i++){full_crc ^= io_crc[i];}
    return full_crc;
  }

  ////fd = &fds;
  io_vec(const qlat::Geometry& geo,int ionum_or,bool checkorder=true, int threadio_set = IO_THREAD, bool do_checksum_set=false){
  /////x,y,z,t
  geop = geo;
  threadio = threadio_set;
  MPI_size_c = 0;tmp = NULL;
  do_checksum = do_checksum_set;
  end_of_file = 0;

  //int Nv = omp_get_max_threads();
  //if(threadio < 1 or threadio > Nv){threadio = Nv;}

  nv.resize(4);Nv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];

  vol   =  nx*ny*nz*nt;
  noden =  Nx*Ny*Nz*Nt;

  Nmpi  = qlat::get_num_node();
  rank  = qlat::get_id_node();

  map_Nitoi.resize(Nmpi);
  for(int ri=0;ri<Nmpi;ri++){
    map_Nitoi[ri].resize(noden);
  }

  for(size_t isp=0;isp<size_t(noden);isp++){
    qlat::Coordinate ts = geo.coordinate_from_index(isp);
    qlat::Coordinate gs = geo.coordinate_g_from_l(ts);
    LInt offv = ((gs[3]*nz+gs[2])*ny+gs[1])*nx+gs[0];
    map_Nitoi[rank][isp] = offv;
  }

  for(int ri=0;ri<Nmpi;ri++)MPI_Bcast(&map_Nitoi[ri][0], noden*sizeof(LInt), MPI_CHAR, ri, get_comm());

  if(checkorder){
    int flag = 0;
    for(int i=0;i<Nmpi;i++)
    for(size_t isp=0;isp<size_t(noden/Nx);isp++){
      LInt inix = map_Nitoi[i][isp*Nx+ 0];
      for(int xi=0;xi<Nx;xi++)if(map_Nitoi[i][isp*Nx+xi] != inix + xi){flag=1;}
    }
    sum_all_size(&flag,1);
    if(flag>0){abort_r("Layout not continuous in x! \n");}
  }

  node_ioL.resize(Nmpi);
  if(ionum_or > 0 and ionum_or < Nmpi){ionum=ionum_or;}else{
    std::string val = get_env(std::string("q_io_vec_ionum"));
    if(val == ""){ionum =32;}else{ionum = stringtonum(val);}
    if(ionum < 0 or ionum > Nmpi){ionum = Nmpi;}
    ////print0("==io number %d \n", ionum);
  }

  for(int ni=0;ni<Nmpi;ni++){node_ioL[ni]=-1;}
  int off = Nmpi/ionum;
  int countN = 0;
  for(int ni=0;ni<Nmpi;ni++){
    if(ni%off==0 and countN < ionum){node_ioL[ni]=countN;countN++;}
  }

  ///io_crc.resize(ionum);
  ini_crc(true);

  buf = NULL;size_buf = 0;
  }


  ~io_vec(){
    //if(node_ioL[rank]>=0){if(tmp != NULL){delete []tmp;tmp=NULL;}}
    if(node_ioL[rank]>=0){if(tmp != NULL){free(tmp);tmp=NULL;}}
    map_Nitoi.resize(0);
    node_ioL.resize(0);
    currsend.resize(0);
    currrecv.resize(0);
    currspls.resize(0);
    currrpls.resize(0);
    io_crc.resize(0);
    ////geop = NULL;
    if(buf != NULL){free(buf);size_buf = 0;buf = NULL;}
  }


};

inline void send_vec_kentucky(char* src,char* res,int dsize,int gN, io_vec& io, bool read=true)
{
  const std::vector<int >& node_ioL = io.node_ioL;
  int rank = io.rank;size_t vol = io.vol;size_t noden = io.noden;
  ///int Nmpi=io.Nmpi;

  int Nx = io.Nx;

  int size_c0 = gN*dsize;
  ////int size_c = gN*noden*dsize;

  ////char* tmp=NULL;
  io.ini_MPI(size_c0);
  MPI_Barrier(get_comm());

  ///position p;
  if(read==true)
  {
  TIMERA("IO sort mem time");
  if(node_ioL[rank]>=0)
  for(int ri=0;ri<io.Nmpi;ri++)for(int gi=0;gi<gN;gi++){
      #pragma omp parallel for
      for(size_t isp=0;isp<size_t(noden/Nx);isp++){
        char* pres=&io.tmp[(ri*gN+gi)*noden*dsize];
        char* psrc=&src[gi*vol*dsize];
        size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
        memcpy(&pres[isp*Nx*dsize],&psrc[offv*dsize],Nx*dsize);
      }
  }
  }

  {
  TIMERA("IO MPI time");
  /////for(int io=0;io<node_ioL.size();io++)if(node_ioL[io]>=0)
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
     (char*) io.tmp, &io.currsend[0], &io.currspls[0], 
     (char*)    res, &io.currrecv[0], &io.currrpls[0], tem_comm);}

  if(read==false)
  {MPI_Alltoallv_Send_Recv(
     (char*)    res, &io.currrecv[0], &io.currrpls[0],
     (char*) io.tmp, &io.currsend[0], &io.currspls[0], tem_comm);}

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
  for(int ri=0;ri<io.Nmpi;ri++)for(int gi=0;gi<gN;gi++){
    #pragma omp parallel for
    for(size_t isp=0;isp<size_t(noden/Nx);isp++){
      char* pres=&io.tmp[(ri*gN+gi)*noden*dsize];
      char* psrc=&src[gi*vol*dsize];
      size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
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
inline void read_kentucky_vector(FILE *file,char* props,int Nvec,io_vec& io,bool Rendian=false, int dsize=8, bool single_file=false,  int gN=1, bool read = true)
{
  TIMER("IO VECS");
  timeval tm0,tm1,tm2,tm3;
  gettimeofday(&tm0, NULL);
  double mpi_t = 0.0;

  const std::vector<int >& node_ioL = io.node_ioL;
  int rank = io.rank;size_t vol = io.vol;int ionum=io.ionum;size_t noden=io.noden;

  if(io.do_checksum){if(io.end_of_file == 0){abort_r("io_vec need end of file for check sum! \n ");}}

  ////char* buf=NULL;if(node_ioL[rank]>=0){buf = new char[gN*vol*dsize];}
  char* buf=NULL;size_t size_buf = gN*vol*dsize;
  if(node_ioL[rank]>=0){
    if(io.size_buf != size_buf){
      if(io.buf != NULL){free(io.buf);io.size_buf = 0;io.buf = NULL;}
      io.buf = (char *)aligned_alloc_no_acc(size_buf);io.size_buf = size_buf;}
  }
  buf = io.buf;
  char res[ionum*noden*gN*dsize];
  size_t sizec = vol*dsize;

  size_t off_file = 0;////size_t sem=0;
  std::vector<size_t > semL;semL.resize(io.threadio);for(LInt i=0;i<semL.size();i++){semL[i]=0;}
  if(node_ioL[rank]>=0){off_file = ftell(file);}

  int curr_v = 0;
  for(int readi=0;readi<Nvec;readi++){

    /////From res to buf
    if(read==false)
    {
    TIMERA("IO copy mem time");
    for(int iou=0;iou<ionum;iou++)
    for(int gi=0;gi<gN;gi++)if(curr_v + iou*gN + gi<Nvec){
      int offv = curr_v + iou*gN + gi;
      //memcpy(&res[(iou*gN+gi)*noden*dsize + 0],&props[(offv)*noden*dsize + 0],noden*dsize);

      char* pres=&res[(iou*gN+gi)*noden*dsize + 0];
      char* psrc=&props[(offv)*noden*dsize + 0];
      #pragma omp parallel for
      for(size_t isp=0;isp<size_t(noden*dsize);isp++){pres[isp] = psrc[isp];}
    }
      gettimeofday(&tm2, NULL);
      send_vec_kentucky((char*) &buf[0],(char*) &res[0], dsize,gN, io, read);
      gettimeofday(&tm3, NULL);
      mpi_t += tm3.tv_sec - tm2.tv_sec;
      mpi_t += (tm3.tv_usec - tm2.tv_usec)/1000000.0;
    }


    {
    TIMERA("IO disk");
    if(node_ioL[rank]>=0)
    if(curr_v + node_ioL[rank]*gN < Nvec)
    {
      size_t offv = curr_v + node_ioL[rank]*gN;
      fseek ( file , off_file + offv*sizec , SEEK_SET );
      size_t pos_file_cur = off_file + offv*sizec;

      int rN = gN;if(offv+rN>=size_t(Nvec)){rN = Nvec-offv;}

      /////Switch endian of the file write
      if(read==false){
        int dsize_single = sizeof(double);int fac = dsize/sizeof(double);
        if(single_file == true){dsize_single = sizeof(float);fac = dsize/sizeof(float);}
        //////Write Prop and link
        if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
        //////Write source noise
        if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
      }

      if(rN>0){
        int curr_threadio = io.threadio;
        if(curr_threadio==1){
          //if(read==true ) fread(&buf[0], 1, rN*sizec, file);
          //if(read==false)fwrite(&buf[0], 1, rN*sizec, file);
          if(read==true ){semL[0] =  fread(&buf[0], rN*sizec, 1, file);}
          if(read==false){semL[0] = fwrite(&buf[0], rN*sizec, 1, file);}
          if(semL[0] != 1){printf("Reading/Writing error %zu 1 \n", semL[0] );}
        }


        ////omp reading
        if(curr_threadio>1){
          size_t Nvec = rN*sizec;
          size_t Group = (Nvec-1)/curr_threadio + 1;
          size_t off_file_local = ftell(file);
          #pragma omp parallel for
          for(int tid=0;tid<curr_threadio;tid++)
          {
            //FILE *file_omp = file;
            /////FILE *fp2 = fdopen (dup (fileno (fp)), "r");
            size_t currN = Group; if((tid+1)*Group > Nvec){currN = Nvec - tid*Group;}
            if(currN > 0){
              size_t iniN  = tid*Group;//size_t endN = tid*Group + currN;
              //FILE *file_omp = fdopen (dup (fileno (file)), "r");
              fseek(io.file_omp[tid] , off_file_local + iniN, SEEK_SET );

              //if(read==true ) fread(&buf[iniN], 1, currN, io.file_omp[tid]);
              //if(read==false)fwrite(&buf[iniN], 1, currN, io.file_omp[tid]);

              if(read==true ){semL[tid] =  fread(&buf[iniN], currN, 1, io.file_omp[tid]);}
              if(read==false){semL[tid] = fwrite(&buf[iniN], currN, 1, io.file_omp[tid]);}
            }
          }
          for(int tid=0;tid<curr_threadio;tid++){if(semL[tid] != 1){printf("Reading/Writing error %zu 1 \n", semL[tid] );}}
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
        int dsize_single = sizeof(double);int fac = dsize/sizeof(double);
        if(single_file == true){dsize_single = sizeof(float);fac = dsize/sizeof(float);}
        //////Read Prop and link
        if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
        //////Read source noise
        if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
      }

    }
    }

    ////////Reordersrc buf
    ///if(gN!=1){if(node_ioL[rank]>=0)reorder_civ((char*)&buf[0],(char*)&buf[0],1,gN,vol,0 ,sizeof(double));}
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
    for(int iou=0;iou<ionum;iou++)
    for(int gi=0;gi<gN;gi++)if(curr_v + iou*gN + gi<Nvec){
      int offv = curr_v + iou*gN + gi;
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

  if(read==true) print0("READ_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);
  if(read==false)print0("WRIT_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);

  MPI_Barrier(get_comm());
  fflush(stdout);
}

inline void save_gwu_eigenvalues(std::vector<double > &values,std::vector<double > &errors,const char* filename, const char* sDescription)
{
  if(qlat::get_id_node() == 0)
  {
    int nvec = values.size()/2;
    FILE* filew = fopen(filename, "w");

    fprintf(filew, "Eigenvalues and eigenvectors for the %s\n", sDescription);
    fprintf(filew, "Each eigenvector is preceded by a line describing the eigenvalue.\n");
    fprintf(filew, "The residue is defined as norm(mat.vec-lambda.vec).\n");
    fprintf(filew, "The format is: a tag EIGV, the real and imaginary part of the eigenvalue and the residue.\n");

    for(int iv=0;iv<nvec;iv++)fprintf(filew, "EIGV %+.15le\t%+.15le\t%.10le\n", values[iv*2+0],values[iv*2+1],errors[iv]);

    fclose(filew);
  }
}

inline void load_gwu_eigenvalues(std::vector<double > &values,std::vector<double > &errors,const char* filename)
{
  int nvec = 0;
  values.resize(0);errors.resize(0);
  if(qlat::get_id_node() == 0)
  {
    FILE* filer = fopen(filename, "r");
    char sTemp[300],tem[300];
    char* ftem = NULL;
    for(int i=0;i<4;i++){
      ftem = fgets(sTemp, 300, filer);
      if(ftem == NULL){print0("Read eigenvalues error!\n");}
    }

    //while(!feof(filer))
    while(fgets(tem, 300, filer) != NULL){
      //fgets(tem, 300, filer);
      std::string re = std::string(tem);
      int size = re.size();
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

//
//order C, iDir, C_col, C_row, RealIm, t,z,y,x
  
template<typename Ty>
void rotate_gwu_vec_file(Ty* src,int n_vec,size_t noden,bool single_file,bool read=true){
  TIMERB("Rotate gwu file vec");

  if(single_file == true){
    return;
  }

  //size_t noden = (*props[0]).desc->sites_on_node;
  size_t Np=noden*12*2;
  for(int dc1=0;dc1<n_vec;dc1++){
    std::vector<Ty > p;Ty *q;
    p.resize(Np);
    //memcpy(&p[0],&src[dc1*Np],sizeof(Ty)*p.size());
    #pragma omp parallel for
    for(size_t isp=0;isp<p.size();isp++){p[isp] = src[dc1*Np + isp];}

    q = &src[dc1*Np];
    if(read==true){
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++){
      for(int dc0=0;dc0<12;dc0++)
      {   
        q[isp*12*2 + dc0*2 + 0] = p[(0*12 + dc0)*noden + isp];
        q[isp*12*2 + dc0*2 + 1] = p[(1*12 + dc0)*noden + isp];
      }
    }}

    if(read==false){
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++){
      for(int dc0=0;dc0<12;dc0++)
      {   
        q[(0*12 + dc0)*noden + isp] = p[isp*12*2 + dc0*2 + 0];
        q[(1*12 + dc0)*noden + isp] = p[isp*12*2 + dc0*2 + 1];
      }
    }}
  }
}

///////Ty is real as the transformation have only real numbers
template<typename Ty>
void gwu_to_cps_rotation_vec(Ty* src,int n_vec,size_t noden,bool source=false,bool PS0=true,bool PS1=true, bool read=true){
  TIMERB("Rotate gwu file vec");
  ///prop
  if(source == true)if(n_vec%12!=0){print0("n_vec %8d.\n",n_vec);abort_r("source vector size wrong.\n");}
  //size_t noden = geo.node_site[0]*geo.node_site[1]*geo.node_site[2]*geo.node_site[3];

  int signT = 1;if(read==false){signT=-1;}
  //std::vector<Typ0 > q;q.resize(12*2)
  const double sqrt2=std::sqrt(2.0);
  //////Rotate source dirac
  if(source == true){
    int dr,d0,d1;Ty *q;
    int Nprop = n_vec/12;
    
    size_t Np=noden*12*2;
    for(int ip=0;ip<Nprop;ip++){
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
  for(int dc1=0;dc1<n_vec;dc1++){
  #pragma omp parallel for
  for(size_t isp=0;isp<noden;isp++){
    int dr,d0,d1;Ty *q;
    std::vector<Ty > p;p.resize(12*2);
    q = &src[(dc1*noden + isp)*12*2 + 0];

    if(PS1 == true){
      memcpy(&p[0],&src[(dc1*noden + isp)*12*2 + 0],sizeof(Ty)*p.size());
      dr=0;d0=1;d1=3;for(int c=0;c<6;c++)q[dr*6+c] = (-p[d0*6+c]*signT + p[d1*6+c])/sqrt2;
      dr=1;d0=0;d1=2;for(int c=0;c<6;c++)q[dr*6+c] = ( p[d0*6+c]*signT - p[d1*6+c])/sqrt2;
      dr=2;d0=1;d1=3;for(int c=0;c<6;c++)q[dr*6+c] = (-p[d0*6+c] - p[d1*6+c]*signT)/sqrt2;
      dr=3;d0=0;d1=2;for(int c=0;c<6;c++)q[dr*6+c] = (+p[d0*6+c] + p[d1*6+c]*signT)/sqrt2;
    }else{
      dr = 2;for(int c=0;c<6;c++)q[dr*6+c] = -1.0*q[dr*6+c];
      dr = 3;for(int c=0;c<6;c++)q[dr*6+c] = -1.0*q[dr*6+c];
    }
  }}

}

template<typename Ty>
Ty get_norm_vec(Ty *src,size_t noden){
  Ty res = 0.0;

  std::complex<Ty > *p = (std::complex<Ty >*) src;
  /////need sum 12 first to reduce float sum error
  #pragma omp parallel for reduction(+: res)
  for(size_t isp=0;isp<noden;isp++)
  {
    Ty a = 0.0;
    for(unsigned int dc=0;dc<12;dc++){
       a += std::norm(p[isp*12+dc]);
    }
    res += a;
  }

  //print0("==omp_get_max_threads %8d \n ",omp_get_max_threads());
  sum_all_size(&res,1);
  return res;

}

inline int test_single(const char *filename,io_vec &io_use,int iv=0){
  
  double normd=0.0;double normf=0.0;
  int n_vec = 1;
  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  FILE* file;
  double err = 1e-3;

  {
  std::vector<double > prop_E;
  prop_E.resize(n_vec*12*noden*2);

  file = io_use.io_read(filename,"rb");
  io_use.io_off(file,iv*Fsize*2, false);
  read_kentucky_vector(file,(char*) &prop_E[0], n_vec*12*2,io_use, false, sizeof(double), false, 12*2);
  io_use.io_close(file);

  rotate_gwu_vec_file(&prop_E[0],n_vec,noden, false);
  gwu_to_cps_rotation_vec(&prop_E[0],n_vec,noden, false, true, true);

  normd = get_norm_vec(&prop_E[0],noden);
  }

  if(fabs(normd - 1.0) < err)return 0;

  {
  std::vector<float > prop_E;
  prop_E.resize(n_vec*12*noden*2);
  file = io_use.io_read(filename,"rb");
  io_use.io_off(file,iv*Fsize, false);

  read_kentucky_vector(file,(char*) &prop_E[0], n_vec,io_use, true, 6*2*sizeof(double), true);
  io_use.io_close(file);

  rotate_gwu_vec_file(&prop_E[0],n_vec,noden, true);
  gwu_to_cps_rotation_vec(&prop_E[0],n_vec,noden, false, false, false);

  normf = get_norm_vec(&prop_E[0],noden);
  }


  if(fabs(normf - 1.0) < err)return 1;

  print0("Norm of vector double %.3e, %.3e.\n",normd,normd-1.0);
  print0("Norm of vector single %.3e, %.3e.\n",normf,normf-1.0);
  return -1;
}
/////================END of NO checksum read/write

inline int check_Eigen_file_type(const char *filename, io_vec &io_use,int n1,bool check){

  size_t noden = io_use.noden;
  size_t sizen = get_file_size_MPI(filename);

  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  if(sizen <  Fsize){print0("%s \n",filename);abort_r("Eigen file size too small. \n");}
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

  int tmp = test_single(filename,io_use,0);
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
void load_gwu_eigen(FILE* file,std::vector<Ty* > resp,io_vec &io_use,int n0,int n1,
  bool check=true, bool read=true, bool read_single=true)
{
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(resp.size() < size_t(n1-n0)){abort_r("Final point size wrong!\n");}
  if(sizeof(Ty) != sizeof(double) and sizeof(Ty) != sizeof(float)){abort_r("Need double or float pointer! \n");}
  ////if(n0!=0 and read==false){abort_r("Write vectors with off not zero!");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  bool single = read_single;
  //bool single = true;
  //////read==false, write only single vectors
  //if(read==true){
  //  int type = check_Eigen_file_type(filename,io_use,n1,check);
  //  if(type == 0){single = false;}
  //  if(type == 1){single = true ;}
  //  if(type !=0 and type != 1){abort_r("Eigen system not in gwu format! \n");}
  //}
  ////FILE* file;

  if(single == true){
    /////Load single precision
      
    int count = 0;int off = io_use.ionum;
    std::vector<float > prop_E;
    int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    ////if(read==true )file = io_use.io_read(filename,"rb");
    ////if(read==false)file = io_use.io_read(filename,"wb");
    io_use.io_off(file, n0*Fsize, true);

    for(int iv=0;iv<n_vec;iv++){
      int ri = off;if(count + off > n_vec){ri = n_vec - count;}
      //int offE = count*12*noden*2;

      if(read==false){
        for(int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
           prop_E[ic*noden*12*2 + isp] = resp[count + ic][isp];
        }
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, false, false, read);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, true, read);
      }

      /////print0("WRITE!!!!!!\n");
      read_kentucky_vector(file,(char*) &prop_E[0], ri,io_use, true, 6*2*sizeof(double), true, 1, read);
      fflush_MPI();

      if(read==true){
        //rotate_qlat_to_gwu(prop_E, &src_new.vec[0],geo, true);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, true, read);
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, false, false, read);
        for(int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
          resp[count + ic][isp] = prop_E[ic*noden*12*2 + isp];}
      }

      count += ri;if(count >= n_vec){break;}
    }
    /////io_use.io_close(file);

    /////rotate_qlat_to_gwu(prop_E,&src_new.vec[0],geo);
    /////////From ky to milc
    /////for(int iv=0;iv<n_vec;iv++)ga.ga[4][0].multiply(*src_new.vec[iv],*src_new.vec[iv]);
  }else{
    /////Load double precision

    int count = 0;int off = io_use.ionum;
    std::vector<double > prop_E;
    int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    //////file = io_use.io_read(filename,"rb");
    io_use.io_off(file,n0*Fsize*2, true);
      
    for(int iv=0;iv<n_vec;iv++){
      int ri = off;if(count + off > n_vec){ri = n_vec - count;}
      //int offE = count*12*noden*2;

      if(read==false){
        for(int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
          prop_E[ic*noden*12*2 + isp] = resp[count + ic][isp];
        }
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, true, true, read);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, false, read);
      }

      read_kentucky_vector(file,(char*) &prop_E[0], ri*12*2,io_use, false, sizeof(double), false, 12*2, read);
      fflush_MPI();

      if(read==true){
        /////double precision eigen vector in ps base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, false, read);
        ////Do not rotate source, rotate sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, true, true, read);

        for(int ic=0;ic<ri;ic++)
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
//void load_gwu_eigen(FILE* file, std::vector<qlat::FieldM<Ty , 12> > res,io_vec &io_use,int n0,int n1,
//  bool check=true, bool read=true, bool read_single=true)
//{
//  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
//  int n_vec = n1-n0;
//  size_t noden = io_use.noden;
//  if(read == true)
//  if(res.size() != n_vec)
//  {
//    res.resize(0);
//    res.resize(n_vec);
//    for(int iv=0;iv<n_vec;iv++){res[iv].init(io_use.geop);}
//  }
//  std::vector<Ty* > resp;resp.resize(n_vec);
//  for(int iv=0;iv<n_vec;iv++){resp[iv] = (Ty*) qlat::get_data(res[iv]).data();}
//  load_gwu_eigen(file, resp, io_use,n0,n1,check, read, read_single);
//}

inline FILE* open_gwu_eigen(const char *filename,io_vec &io_use, bool read=true)
{
  if(read==true ){return io_use.io_read(filename,"rb");}
  if(read==false){return io_use.io_read(filename,"wb");}
  return NULL;
}

template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty* > resp,io_vec &io_use,int n0,int n1,
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
void load_gwu_eigen(FILE *file,std::vector<Ty > &res,io_vec &io_use,int n0,int n1,bool check=true,bool read=true, bool read_single=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  int n_vec = n1-n0;
  size_t noden = io_use.noden;
  if(read == true)res.resize(n_vec*noden*12*2);

  std::vector<Ty* > resp;resp.resize(n_vec);
  for(int iv=0;iv<n_vec;iv++){resp[iv] = &res[iv*noden*12*2];}
  load_gwu_eigen(file, resp, io_use,n0,n1,check,read, read_single);
}

//////Ty should be float or double
template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty > &res,io_vec &io_use,int n0,int n1,bool check=true,bool read=true){
  /////FILE* file = NULL;
  bool read_single = true;
  if(read==true){read_single = check_Eigen_file_type(filename, io_use, n1, check);}
  FILE* file = open_gwu_eigen(filename, io_use, read);
  load_gwu_eigen(file, res, io_use,n0,n1,check,read, read_single);
  io_use.io_close(file);
}


template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &eigen,io_vec &io_use,int n0,int n1,bool check=true,bool read=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  int n_vec = n1-n0;
  if(read==true){
    eigen.resize(0);
    eigen.resize(n_vec);
    for(int iv=0;iv<n_vec;iv++)eigen[iv].init(io_use.geop);
  }

  if(sizeof(Ty) == 2*sizeof(double )){
    std::vector<double* > resp;resp.resize(n_vec);
    for(int iv=0;iv<n_vec;iv++){
      resp[iv]=(double*)(qlat::get_data(eigen[iv]).data());
    }
    load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
  }

  if(sizeof(Ty) == 2*sizeof(float )){
    std::vector<float* > resp;resp.resize(n_vec);
    for(int iv=0;iv<n_vec;iv++){
      resp[iv]=(float*)(qlat::get_data(eigen[iv]).data());
    }
    load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
  }

  ////return eigen;
}

template<typename Ty>
void save_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &eigen,io_vec &io_use,int n0,int n1,bool check=true){
  load_gwu_eigen(filename,eigen,io_use,n0,n1,check,false);
}

inline void load_gwu_eigen(const char *filename,EigenM &Mvec,io_vec &io_use,int n0,int n1,bool check=true,bool read=true)
{
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  int n_vec = n1-n0;
  if(n_vec > int(Mvec.size())){abort_r("Read number of eigen larger than memory. \n");}

  long Nsize = io_use.noden*12;

  std::vector<Ftype* > resp;resp.resize(n_vec);
  for(int iv=0;iv<n_vec;iv++){
    if(Mvec[iv].size() != Nsize){
      if(read==false)abort_r("Mvec not initialized! \n");
      if(read==true ){Mvec[iv].resize(Nsize);set_zero(Mvec[iv]);}
    }
    resp[iv]= (Ftype*)(&Mvec[iv][0]);
  }
  load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);
}

inline void save_gwu_eigen(const char *filename,EigenM &Mvec,io_vec &io_use,int n0,int n1,bool check=true)
{
  load_gwu_eigen(filename, Mvec,io_use,n0,n1,check,false);
}


template<typename Ty>
void load_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &prop,io_vec &io_use,bool read=true){

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);

  FILE* file;bool single = true;
  if(read==true){
  size_t sizen = get_file_size_MPI(filename);

  if(sizen != 2*Fsize*12 and sizen != Fsize*12){print0("File %s \n",filename);abort_r("prop size wrong! \n");}
  prop.resize(0);
  prop.resize(12);
  for(int iv=0;iv<12;iv++)prop[iv].init(io_use.geop);
  if(sizen == 2*Fsize*12){single=false;}
  }
    
  ////Can only write with single!
  if(read==false){single = true;}
  
  if(single == true)
  {
    ////Single vector read
    std::vector<float > prop_qlat;
    prop_qlat.resize(12*noden*12*2);

    if(read==false){
    for(int iv=0;iv<12;iv++){
      Ty* res   = (Ty*) qlat::get_data(prop[iv]).data();
      std::complex<float> *src = (std::complex<float>*) &prop_qlat[iv*noden*12*2];
      #pragma omp parallel for
      for(size_t isp=0;isp<noden*12;isp++)src[isp] = std::complex<float>(res[isp].real(),res[isp].imag());
    }

    ////Do not rotate source, in ps/ky base
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true, true, read);
    /////double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, true, read);
    }

    if(read==true )file = io_use.io_read(filename,"rb");
    if(read==false)file = io_use.io_read(filename,"wb");
    read_kentucky_vector(file,(char*) &prop_qlat[0], 12,io_use, true, 6*2*sizeof(double), true, 1 , read);
    io_use.io_close(file);
      
    if(read==true){
    /////double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, true);
    ////Do not rotate source, in ps/ky base
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true, true);

    for(int iv=0;iv<12;iv++){
      Ty* res   = (Ty*) qlat::get_data(prop[iv]).data();
      std::complex<float> *src = (std::complex<float>*) &prop_qlat[iv*noden*12*2];
      #pragma omp parallel for
      for(size_t isp=0;isp<noden*12;isp++)res[isp]= Ty(src[isp].real(),src[isp].imag());
    }
    }
  }

  if(single == false)
  {
    std::vector<double > prop_qlat;
    prop_qlat.resize(12*noden*12*2);

    file = io_use.io_read(filename,"rb");
    read_kentucky_vector(file,(char*) &prop_qlat[0], 12*2*12,io_use, false, sizeof(double), false, 1);
    io_use.io_close(file);
      
    ///double precision eigen vector in ps base
    rotate_gwu_vec_file(&prop_qlat[0], 12,noden, false);
    //Do not rotate source, 
    gwu_to_cps_rotation_vec(&prop_qlat[0], 12,noden, true, true,true);

    for(int iv=0;iv<12;iv++){
      Ty* res = (Ty*) qlat::get_data(prop[iv]).data();
      std::complex<double> *src = (std::complex<double>*) &prop_qlat[iv*noden*12*2];
      for(size_t isp=0;isp<noden*12;isp++)res[isp]= Ty(src[isp].real(),src[isp].imag());
    }

  }
}

template<typename Ty>
void save_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &prop,io_vec &io_use){
  load_gwu_prop(filename,prop,io_use,false);
}

//////final result 12*12 --> Nt*Nxyz
template<typename Ty>
void load_gwu_prop(const char *filename, qpropT& res,io_vec &io_use,bool read=true){
  if(read == true ){res.init(io_use.geop);}
  if(read == false){abort_r("Not supported! \n");}

  long sizeF = io_use.geop.local_volume();

  std::vector<qlat::FermionField4dT<Ty> > prop;
  load_gwu_prop(filename, prop, io_use, read);

  move_index mv_civ;
  Ty* p0; Ty* p1;Ty* pt;
  pt = (Ty*) qlat::get_data(res).data();

  for(int iv=0;iv<12;iv++){
    p0 = (Ty*) qlat::get_data(prop[iv]).data();
    p1 = (Ty*) &pt[iv*12*sizeF + 0];
    mv_civ.dojob(p0, p1, 1, 12, sizeF, 1, 1, false);
  }
}

/////V -- 12a x 12b   to   12b x 12a -- V
template<class T, typename Ty>
void prop4d_to_qprop(qpropT& res, Propagator4dT<T>& src, int dir = 1){
  TIMERA("prop4d_to_qprop");
  if(dir == 1){qassert(src.initialized);res.init(src.geo());}
  if(dir == 0){qassert(res.initialized);src.init(res.geo());}

  long sizeF = src.geo().local_volume();

  move_index mv_civ;
  T* ps; Ty* pt;
  ps = (T* ) qlat::get_data(src).data();
  pt = (Ty*) qlat::get_data(res).data();

  ////V x 12 a x 12 b to 12b x 12a x V
  if(dir == 1){
    qthread_for(isp, long(sizeF),{
      T buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = ps[isp*12*12 + i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        pt[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
    mv_civ.move_civ_out(pt, pt, 1, sizeF, 12*12, 1, false);
  }
  ////12 a x 12 b x V to V x 12b x 12a
  if(dir == 0){
    mv_civ.move_civ_in(pt, pt, 1, 12*12, sizeF, 1, false);
    qthread_for(isp, long(sizeF),{
      T buf[12*12];for(unsigned int i=0;i<12*12;i++){buf[i] = pt[isp*12*12 + i];}
      for(unsigned int d0=0;d0<12;d0++)
      for(unsigned int d1=0;d1<12;d1++)
      {
        ps[(isp*12+d0)*12+d1] = buf[d1*12+d0];
      }
    });
  }
}

template<class T, typename Ty>
void qprop_to_prop4d(Propagator4dT<T>& res, qpropT& src){
  prop4d_to_qprop(src, res, 0);
}

template<class T, typename Ty>
void prop4d_to_Fermion(Propagator4dT<T>& prop,std::vector<qlat::FermionField4dT<Ty > > &buf, int dir=1){

  ////if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}
  if(dir==1){buf.resize(0);buf.resize(12);for(int iv=0;iv<12;iv++){buf[iv].init(prop.geo());}}
  if(dir==0){qassert(buf.size() == 12);prop.init(buf[0].geo());}

  #pragma omp parallel for
  for (long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrixT<T>& src =  prop.get_elem(index);
    for(int d0=0;d0<12;d0++)
    {
      ////v0(s*3 + c0, ga.ind[d0]*3 + c1)
      Ty* res = (Ty*)&(buf[d0].get_elem(index));
      for(int d1=0;d1<12;d1++)
      {
        if(dir==1){res[d1] = src(d1, d0);}
        if(dir==0){src(d1, d0) = res[d1];}
      }

    }
  }

}

template<class T, typename Ty>
void Fermion_to_prop4d(Propagator4dT<T>& prop, std::vector<qlat::FermionField4dT<Ty > > &buf){
  qassert(buf.size() == 12);
  prop4d_to_Fermion(prop, buf, 0);
}


template <class T>
void save_gwu_prop(const char *filename,Propagator4dT<T>& prop){
  qassert(prop.initialized);
  io_vec io_use(prop.geo(),IO_DEFAULT);
  std::vector<qlat::FermionField4dT<qlat::ComplexF> > prop_qlat;
  prop4d_to_Fermion(prop,prop_qlat, 1);
  save_gwu_prop(filename,prop_qlat,io_use);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <class T>
void save_gwu_prop(std::string &filename,Propagator4dT<T>& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  save_gwu_prop(tem,prop);
}

template <class T>
void load_gwu_prop(const char *filename,Propagator4dT<T>& prop){
  qassert(prop.initialized);
  io_vec io_use(prop.geo(),IO_DEFAULT);
  std::vector<qlat::FermionField4dT<qlat::Complex > > prop_qlat;
  load_gwu_prop(filename,prop_qlat,io_use);
  prop4d_to_Fermion(prop,prop_qlat, 0);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <class T>
void load_gwu_prop(std::string &filename,Propagator4dT<T>& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_prop(tem,prop);
}



template <class T>
void load_gwu_link(const char *filename,GaugeFieldT<T> &gf, bool read = true){
  io_vec io_use(gf.geo(),8);
  //if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ))
  //{abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(4*9*noden*2)*sizeof(double);

  std::vector<double > link_qlat;
  link_qlat.resize(4*9*noden*2);

  if(read==true)
  {
    size_t sizen = get_file_size_MPI(filename);
    if(sizen != Fsize){abort_r("Link size wrong! \n");}
    /////gf.init(*io_use.geop);
  }


  if(read == false)
  for (size_t index = 0; index < noden; ++index)
  {
    ColorMatrixT<T>& res = gf.get_elem(index*gf.geo().multiplicity+0);

    for(int dir=0;dir<4;dir++)
    {
      for(int c0=0;c0<3;c0++)
      for(int c1=0;c1<3;c1++)
      {
        int dir0  = ((dir*3+c0)*3+c1)    ;
        int dir1R = ((dir*3+c1)*3+c0)*2+0;
        int dir1I = ((dir*3+c1)*3+c0)*2+1;
        link_qlat[dir1R*noden + index] = res.p[dir0].real();
        link_qlat[dir1I*noden + index] = res.p[dir0].imag();
      }
    }
  }

  FILE* file;
  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  read_kentucky_vector(file,(char*) &link_qlat[0], 4*9*2,io_use, false, sizeof(double), false,9*2, read);
  io_use.io_close(file);
    
  //////double precision eigen vector in ps base
  ///rotate_gwu_vec_file(&prop_qlat[0], 12,noden, false);
  /////Do not rotate source, 
  ///gwu_to_cps_rotation_vec(&link_qlat[0], 12,noden, true, true,true);
  ///4 dire --> c0 -- > c1

  if(read == true)
  for (size_t index = 0; index < noden; ++index)
  {
    ColorMatrixT<T>& res = gf.get_elem(index*gf.geo().multiplicity+0);

    for(int dir=0;dir<4;dir++)
    {
      for(int c0=0;c0<3;c0++)
      for(int c1=0;c1<3;c1++)
      {
        int dir0  = ((dir*3+c0)*3+c1)    ;
        int dir1R = ((dir*3+c1)*3+c0)*2+0;
        int dir1I = ((dir*3+c1)*3+c0)*2+1;
        res.p[dir0] = T(link_qlat[dir1R*noden + index], link_qlat[dir1I*noden + index]);
      }
    }
  }
}

template <class T>
void save_gwu_link(const char *filename,GaugeFieldT<T> &gf){
  load_gwu_link(filename, gf, false);
}

template <class T>
void load_gwu_link(std::string &filename,GaugeFieldT<T> &gf){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_link(tem,gf);
}


template<typename Ty>
void load_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi,io_vec &io_use,bool read=true){

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*2)*sizeof(float);

  std::vector<double > prop_noi;
  prop_noi.resize(noden*2);

  if(read==true){
  size_t sizen = get_file_size_MPI(filename);
  if(sizen != 2*Fsize){abort_r("noise size wrong! \n");}

  qlat::Geometry geo = io_use.geop;geo.multiplicity=1;
  noi.init(geo);
  }
  if(read==false){
    std::complex<double> *src = (std::complex<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++)src[isp] = std::complex<double>(res[isp].real(),res[isp].imag());
    //reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(double));
    io_use.mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(double), false);
  }

  FILE* file;

  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  read_kentucky_vector(file,(char*) &prop_noi[0], 2,io_use, true, sizeof(double), false, 1,read);
  io_use.io_close(file);

  /////Copy noise vectors
  if(read==true){
    //reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(double));
    io_use.mv_civ.dojob((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(double), false);
    std::complex<double> *src = (std::complex<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    #pragma omp parallel for
    for(size_t isp=0;isp<noden;isp++)res[isp]= Ty(src[isp].real(),src[isp].imag());
  }

}

template<typename Ty>
void save_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi,io_vec &io_use){
  load_gwu_noi(filename,noi,io_use,false);
}

template <class T>
void save_gwu_noiP(const char *filename,Propagator4dT<T>& prop){
  io_vec io_use(prop.geo(),IO_DEFAULT);
  qlat::FieldM<qlat::Complex,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  
  size_t noden = io_use.noden;
  for (size_t index = 0; index < noden; ++index)
  {
    qlat::WilsonMatrixT<T>&  src =  prop.get_elem(index);
    double sum = 0.0;
    for(int d1=0;d1<12;d1++)
    for(int d0=0;d0<12;d0++)
    {
      sum += std::fabs(src(d1,d0).real());
      sum += std::fabs(src(d1,d0).imag());
    }
    qlat::Complex phase = qlat::Complex(src(0,0).real(),src(0,0).imag());

    if(sum >1e-8)noi.get_elem(index) = 1.0*phase;
  }

  save_gwu_noi(filename,noi ,io_use);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

template <class T>
void save_gwu_noiP(std::string &filename,Propagator4dT<T>& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  save_gwu_noiP(tem,prop);
}

template <class T>
void noi_to_propP(qlat::FieldM<qlat::Complex,1> &noi,Propagator4dT<T>& prop, int dir = 0){
  for (long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrixT<T>& res =  prop.get_elem(index);
    if(dir==0)for(int d0=0;d0<12;d0++){res(d0,d0) = noi.get_elem(index);}
    if(dir==1)for(int d0=0;d0<12;d0++){noi.get_elem(index) = res(d0,d0);}
  }
}

template <class T>
void load_gwu_noiP(const char *filename,Propagator4dT<T>& prop){
  io_vec io_use(prop.geo(),IO_DEFAULT);
  qlat::FieldM<qlat::Complex,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  load_gwu_noi(filename,noi ,io_use);
  prop.init(noi.geo());

  noi_to_propP(noi, prop, 0);
  
}

template <class T>
void load_gwu_noiP(std::string &filename,Propagator4dT<T>& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_noiP(tem,prop);
}

/////================END of NO checksum read/write

//////Assume memory allocated already
template<class T, typename Ty, int bfac>
void copy_noise_to_vec(qlat::FieldM<T, bfac>& noi, Ty* buf, int dir=1)
{
  TIMERB("copy_noise_to_vec");
  T* noiP = (T*) qlat::get_data(noi).data();
  #pragma omp parallel for 
  for (long index = 0; index < noi.geo().local_volume(); ++index)
  {
    for(int bi=0;bi<bfac;bi++){
      if(dir==1){buf[index*bfac+bi] = noiP[index*bfac+bi];}
      if(dir==0){noiP[index*bfac+bi]= buf[index*bfac+bi];}
    }
    //////T* src  = &noiP[];
    //Ty* res = (Ty*) (&(buf[index]));
    //if(dir==0){src[0] = res[0];}
    ////////copy to buf
    //if(dir==1){res[0] = src[0];}
  }
}

inline void open_file_qlat_noisesT(const char *filename, int bfac, inputpara& in, bool read=true, bool single_file=true, int N_noi=-1, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), bool rotate_bfac = true)
{
  in.bsize = sizeof(float );
  in.rotate_bfac = rotate_bfac;
  in.do_checksum = true;
  in.N_noi = N_noi; in.ncur = 0;
  ////long N_noi = 1;
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
    if(in.VECS_TYPE != VECS_TYPE){print0("Noise type wrong, file %s \n", filename);abort_r("");}
    if(in.nvec <= 0){print0("%s \n", filename);abort_r("File noise vector size Wrong! \n");}
    if(in.bfac != 1){qassert(in.bfac == bfac);}
    in.bfac_write = in.bfac;

    if(in.OBJECT != std::string("BEGIN_Vecs_HEAD")){abort_r("File head wrong");}
    int type = get_save_type(in.save_type);
    if(type == 0){in.bsize=sizeof(double);in.single_file=false;}
    if(type == 1){in.bsize=sizeof(float) ;in.single_file=true; }

    //////Check file sizes
    size_t sizen = get_file_size_MPI(filename) - in.off_file;  //qassert(sizen == string_to_size(in.total_size));
    if(sizen != string_to_size(in.total_size)){
      print0("size  %zu %zu .", sizen, string_to_size(in.total_size));
      abort_r("FILE size not match with head !\n");}

    size_t Vsize = size_t(in.nx)*in.ny*in.nz*in.nt*size_t(bfac*2);
    in.N_noi = in.nvec/(bfac/in.bfac_write);
    size_t Fsize = (in.N_noi + 0)*Vsize*in.bsize;  //qassert(Fsize <= string_to_size(in.total_size));
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
    if(in.N_noi <= 0){print0("write noises size zero \n");return;}
    in.nvec = in.N_noi*(bfac/in.bfac_write);

    //////in.nx = io_use.nx;in.ny = io_use.ny;in.nz = io_use.nz;in.nt = io_use.nt;
    if(in.nx == 0 or in.ny == 0 or in.nz == 0 or in.nt == 0){abort_r("Set up input dim first to write!\n");}
    in.bfac = in.bfac_write;in.checksum = 0;
    if(in.single_file==false){in.bsize=sizeof(double);in.save_type = std::string("Double");}
    if(in.single_file==true ){in.bsize=sizeof(float) ;in.save_type = std::string("Single");}

    size_t Fsize = size_t(in.N_noi) * in.nx*in.ny*in.nz*in.nt* size_t(bfac*2);
    Fsize = Fsize*in.bsize;in.total_size = print_size(Fsize);
    ////print0("size of file %zu \n", Fsize);
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
    if(in.checksum != io_use.full_crc){print0("File %s check sum wrong ! \n ", in.filename.c_str());abort_r("");}
  }
}

template <class Ty, int bfac>
void load_qlat_noisesT(FILE* file, std::vector<qlat::FieldM<Ty, bfac> > &noises, io_vec& io_use, inputpara& in, int n0=0, int n1=-1)
{
  bool read        = in.read;
  bool rotate_bfac = in.rotate_bfac;
  int bfac_write   = in.bfac_write;
  size_t bsize     = in.bsize;
  bool single_file = in.single_file;

  if(read == true ){if(n1 == -1){n1 = in.N_noi;}}
  if(read == false){
    if(n1 == -1){n1 = n0 + noises.size();}
    if(n1 != -1){
    if(n1 - n0 > int(noises.size())){
      print0("Give more noises %d, n0 %d, n1 %d !\n", int(noises.size()), n0, n1);
      abort_r();}
    }
  }
  if(n1 <= n0 or n0 < 0){print0("Need read more vectors, n0 %d, n1 %d !\n", n0, n1);abort_r();}
  if(n1 > in.N_noi){print0("Need set input more vectors, n1 %d, N_noi %d !\n", n1, in.N_noi);abort_r();}
  in.ncur = n1;
  int nread = n1 - n0;

  if(read == true){
    Geometry geo;
    Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
    geo.init(total_site, 1);

    if(noises.size() != (LInt) nread){
      noises.resize(0);
      noises.resize(nread);
      for(LInt ni=0;ni<noises.size();ni++)noises[ni].init(geo);
    }

    //////in.N_noi = noises.size();
  }

  /////if(read == false){geo = noises[0].geo();}
  size_t Vsize = size_t(in.nx)*in.ny*in.nz*in.nt*size_t(bfac*2);
  //size_t off_file = in.off_file + n0*Vsize*bsize;
  size_t off_file = size_t(n0)*Vsize*bsize;
  /////print0(" ionum off %zu, n0 %zu, n1 %zu, Vsize %zu, bsize %zu \n", off_file, size_t(n0), size_t(n1), Vsize, size_t(bsize));
  io_use.io_off(file, off_file, true);

  int io_gn = IO_GN;
  if(in.nvec/io_gn < 1){io_gn = 1;}


  bool Rendian = IO_ENDIAN;
  if(in.FILE_ENDIAN == std::string("BIGENDIAN")){   Rendian = false;}
  if(in.FILE_ENDIAN == std::string("LITTLEENDIAN")){Rendian = true ;}

  /////int bufN = io_use.ionum;
  int bufN = io_use.ionum;
  /////print0("ionum %d %d \n", io_use.ionum, N_noi);

  void* buf;
  buf = aligned_alloc_no_acc(bufN*bfac*io_use.noden * 2*bsize);
  qlat::Complex*  bufD = (qlat::Complex* ) buf;
  qlat::ComplexF* bufF = (qlat::ComplexF*) buf;

  //////false big endian, true small endian
  int bi = 0;
  for(int ni = 0; ni < nread; ni++)
  {
    if(read==false){
      if(!single_file)copy_noise_to_vec(noises[ni], &bufD[bi*bfac*io_use.noden], 1);
      if( single_file)copy_noise_to_vec(noises[ni], &bufF[bi*bfac*io_use.noden], 1);}
    bi = bi + 1;

    if(bi == bufN or ni == (nread - 1)){
    if(read==false)if(rotate_bfac)io_use.mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 1, 2*bsize, false);

    read_kentucky_vector(file,(char*) buf, bi*bfac/bfac_write, io_use, Rendian, bfac_write*bsize*2, single_file, io_gn , read);

    if(read==true)if(rotate_bfac)io_use.mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 0, 2*bsize, false);
    if(read==true)for(int nbi=0; nbi < bi; nbi++){int na = ni - bi + 1;
    {
      if(!single_file)copy_noise_to_vec(noises[na + nbi], &bufD[nbi*bfac*io_use.noden], 0);
      if( single_file)copy_noise_to_vec(noises[na + nbi], &bufF[nbi*bfac*io_use.noden], 0);
    }}

    bi = 0;
    }

  }

  bufD = NULL; bufF = NULL;
  free(buf);

}

//{
//  io_vec io_use(geo, IO_DEFAULT, true, IO_THREAD, do_checksum);
//}

template <class Ty, int bfac>
void load_qlat_noisesT(const char *filename, std::vector<qlat::FieldM<Ty, bfac> > &noises, bool read=true, bool single_file=true, const std::string& VECS_TYPE = std::string("NONE"), const std::string& INFO_LIST = std::string("NONE"), int n0=0,int n1=-1, bool rotate_bfac = true){
  TIMERB("load_qlat_noisesT kernel");
  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){
    abort_r("Cannot understand the input format! \n");}

  long N_noi = 0;
  inputpara in; Geometry geo;
  if(read == false){
    ////if(n0 != 0 or n1 != -1){abort_r("Write mode shoude have n0 0, n1 -1 . ! \n");}
    N_noi = noises.size();
    geo = noises[0].geo();
    in.read_geo(geo);
  }

  open_file_qlat_noisesT(filename, bfac, in, read, single_file, N_noi, VECS_TYPE, INFO_LIST, rotate_bfac);

  if(read == true ){
    Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
    geo.init(total_site, 1);
  }

  io_vec io_use(geo, IO_DEFAULT, true, IO_THREAD, in.do_checksum);
  io_use.end_of_file = in.end_of_file;

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");
  //size_t off_file = in.off_file + n0*Vsize*bsize;
  /////print0(" ionum off %zu, n0 %zu, n1 %zu, Vsize %zu, bsize %zu \n", off_file, size_t(n0), size_t(n1), Vsize, size_t(bsize));

  io_use.io_off(file, in.off_file, true);  ////shift file for the head
  load_qlat_noisesT(file, noises, io_use, in, n0, n1);

  close_file_qlat_noisesT(file, io_use, in);

  //int bsize = sizeof(float );

  //bool do_checksum = true;
  //int N_noi = 1;
  //if(bfac==1)rotate_bfac = false;
  //int bfac_write = bfac;
  //if(rotate_bfac)bfac_write = 1;
  //bool Rendian = IO_ENDIAN;
  //////int Nmpi = qlat::get_num_node();

  //size_t off_file = 0;
  //Geometry geo;
  //if(read == true){
  //  in.load_para(filename, false);

  //  if(in.VECS_TYPE != VECS_TYPE){print0("Noise type wrong, file %s \n", filename);abort_r("");}
  //  Coordinate total_site = Coordinate(in.nx, in.ny, in.nz, in.nt);
  //  geo.init(total_site, 1);
  //  qassert(in.nvec > 0);////qassert(in.bfac == bfac_write);

  //  if(in.bfac != 1){qassert(in.bfac == bfac);}
  //  bfac_write = in.bfac;

  //  if(noises.size() != 0){noises.resize(0);}
  //  int Ntotal = in.nvec/(bfac/bfac_write);
  //  if(nread == -1)noises.resize(Ntotal);
  //  else{
  //    if(nread + nini > Ntotal){print0("Load noise number large, file %s \n", filename);abort_r("");}
  //    noises.resize(nread);
  //  }
  //  for(LInt ni=0;ni<noises.size();ni++)noises[ni].init(geo);
  //  N_noi = noises.size();
  //  if(in.FILE_ENDIAN == std::string("BIGENDIAN")){   Rendian = false;}
  //  if(in.FILE_ENDIAN == std::string("LITTLEENDIAN")){Rendian = true ;}
  //}
  //if(read == false){
  //  in.VECS_TYPE = VECS_TYPE;
  //  in.INFO_LIST  = INFO_LIST;
  //  if(!IO_ENDIAN){in.FILE_ENDIAN = std::string("BIGENDIAN");}
  //  if( IO_ENDIAN){in.FILE_ENDIAN = std::string("LITTLEENDIAN");}
  //  N_noi = noises.size();
  //  in.nvec = N_noi*(bfac/bfac_write);
  //  if(N_noi < 1){print0("write noises size zero \n");return;}
  //  geo = noises[0].geo();

  //}
  //if(bfac_write == bfac)rotate_bfac = false;

  //io_vec io_use(geo, IO_DEFAULT, true, IO_THREAD, do_checksum);
  //int io_gn = IO_GN;
  //if(in.nvec/io_gn < 1){io_gn = 1;}

  //if(read == true){
  //  if(in.OBJECT != std::string("BEGIN_Vecs_HEAD")){abort_r("File head wrong");}
  //  int type = get_save_type(in.save_type);
  //  if(type == 0){bsize=sizeof(double);single_file=false;}
  //  if(type == 1){bsize=sizeof(float) ;single_file=true; }

  //  //////Check file sizes
  //  size_t sizen = get_file_size_MPI(filename) - in.off_file;  //qassert(sizen == string_to_size(in.total_size));
  //  if(sizen != string_to_size(in.total_size)){abort_r("FILE size not match with head !\n");}

  //  size_t Vsize = io_use.Nmpi*size_t(io_use.noden*bfac*2);
  //  size_t Fsize = (N_noi + nini)*Vsize*bsize;  //qassert(Fsize <= string_to_size(in.total_size));
  //  if(Fsize > string_to_size(in.total_size)){abort_r("FILE size too small for vectors read !\n");}
  //  if(Fsize != string_to_size(in.total_size)){io_use.do_checksum = false;}

  //  //////Check file sizes

  //  off_file = in.off_file + nini*Vsize*bsize;
  //  io_use.end_of_file = sizen + in.off_file;
  //}

  //if(read == false){
  //  in.nx = io_use.nx;in.ny = io_use.ny;in.nz = io_use.nz;in.nt = io_use.nt;
  //  in.bfac = bfac_write;in.checksum = 0;
  //  if(single_file==false){bsize=sizeof(double);in.save_type = std::string("Double");}
  //  if(single_file==true ){bsize=sizeof(float) ;in.save_type = std::string("Single");}

  //  size_t Fsize = N_noi* io_use.Nmpi*size_t(io_use.noden*bfac*2);
  //  Fsize = Fsize*bsize;in.total_size = print_size(Fsize);
  //  ////print0("size of file %zu \n", Fsize);
  //  off_file = vecs_head_write(in, filename, true);
  //  io_use.end_of_file = Fsize + off_file;
  //}

  //FILE* file=NULL;
  //if(read==true )file = io_use.io_read(filename,"rb");
  //if(read==false)file = io_use.io_read(filename,"wb");
  //if(io_use.node_ioL[io_use.rank]>=0){fseek(file , off_file, SEEK_SET );} 

  //int bufN = io_use.ionum;
  ///////print0("ionum %d %d \n", io_use.ionum, N_noi);

  //void* buf;
  //buf = aligned_alloc_no_acc(bufN*bfac*io_use.noden * 2*bsize);
  //qlat::Complex*  bufD = (qlat::Complex* ) buf;
  //qlat::ComplexF* bufF = (qlat::ComplexF*) buf;


  ////qlat::vector<qlat::Complex  > bufD;// buf.resize(bufN*bfac*io_use.noden );
  ////qlat::vector<qlat::ComplexF > bufF;// buf.resize(bufN*bfac*io_use.noden );
  ////if(!single_file)bufD.resize(bufN*bfac*io_use.noden);
  ////if( single_file)bufF.resize(bufN*bfac*io_use.noden);

  ////////false big endian, true small endian
  //int bi = 0;
  //for(int ni = 0; ni < N_noi; ni++)
  //{
  //  if(read==false){
  //    if(!single_file)copy_noise_to_vec(noises[ni], &bufD[bi*bfac*io_use.noden], 1);
  //    if( single_file)copy_noise_to_vec(noises[ni], &bufF[bi*bfac*io_use.noden], 1);}
  //  bi = bi + 1;

  //  if(bi == bufN or ni == (N_noi - 1)){
  //  //if(read==false)if(rotate_bfac)reorder_civ((char*) buf,(char*) buf, bi, bfac, io_use.noden, 1, 2*bsize);
  //  if(read==false)if(rotate_bfac)io_use.mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 1, 2*bsize);

  //  ////load_qlat_vecs(filename, &buf[0], io_use, bi*bfac/bfac_write, bfac_write, in, read, single_file);
  //  read_kentucky_vector(file,(char*) buf, bi*bfac/bfac_write, io_use, Rendian, bfac_write*bsize*2, single_file, io_gn , read);

  //  //if(read==true)if(rotate_bfac)reorder_civ((char*) buf,(char*) buf, bi, bfac, io_use.noden, 0, 2*bsize);
  //  if(read==true)if(rotate_bfac)io_use.mv_civ.dojob((char*) buf,(char*) buf, bi, bfac, io_use.noden, 0, 2*bsize);
  //  if(read==true)for(int nbi=0; nbi < bi; nbi++){int n0 = ni - bi + 1;
  //  {
  //    if(!single_file)copy_noise_to_vec(noises[n0 + nbi], &bufD[nbi*bfac*io_use.noden], 0);
  //    if( single_file)copy_noise_to_vec(noises[n0 + nbi], &bufF[nbi*bfac*io_use.noden], 0);
  //  }}

  //  bi = 0;
  //  }

  //}

  //io_use.io_close(file);

  //if(read==false){in.checksum = io_use.full_crc;off_file = vecs_head_write(in, filename, false);}
  //if(read==true and io_use.do_checksum==true){
  //  if(in.checksum != io_use.full_crc){print0("File %s check sum wrong ! \n ", filename);abort_r("");}}

  //bufD = NULL; bufF = NULL;
  //free(buf);

  //if(read==false){
  //  for(int ni=0;ni<N_noi;ni++)copy_noise_to_vec(noises[ni], &buf[ni*bfac*io_use.noden], 1);
  //}
  ///////TODO, need split the jobs if N_noi is too large
  //load_qlat_vecs(filename, &buf[0], io_use, N_noi, bfac_write, in, read, single_file, nini);
  //if(read==true){
  //  for(int ni=0;ni<N_noi;ni++)copy_noise_to_vec(noises[ni], &buf[ni*bfac*io_use.noden], 0);
  //}

}

template <class T>
void load_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, 12> > &noises, bool read , bool single_file=true, const std::string& INFO_LIST = std::string("NONE"),int n0 = 0, int n1=-1)
{
  TIMERC("load/save qlat eigen");
  std::string VECS_TYPE = std::string("Eigen_system_nvec.12.tzyx.R/I");
  ////int nread = -1;if(n1 > 0){qassert(n0 < n1); nread = n1 - n0;}
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST, n0, n1);
}

template <class T>
void load_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, 12> > &noises, int n0=0,int n1=-1)
{
  bool read = true; bool single_file = true; std::string INFO_LIST = std::string("NONE");
  load_qlat_eigen(filename, noises, read, single_file, INFO_LIST, n0, n1);
}

template <class T>
void save_qlat_eigen(const char *filename, std::vector<qlat::FieldM<T, 12> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_eigen(filename, noises, false, single_file, INFO_LIST);
}


////===Check eigen system type
inline int check_eigen_qlat(const char *filename, int n1, inputpara& in)
{
  in.load_para(filename, false);
  if(in.VECS_TYPE != std::string("NONE"))
  {
    if(in.OBJECT != std::string("BEGIN_Vecs_HEAD")){
      print0("%s \n", filename);
      abort_r("File head wrong");
    }

    qassert(in.nvec > 0);////qassert(in.bfac == bfac_write);
    int type = get_save_type(in.save_type);
    in.bsize = 8; in.single_file = true;
    if(type == 0){in.bsize=sizeof(double);in.single_file=false;}
    if(type == 1){in.bsize=sizeof(float) ;in.single_file=true; }
    //////Check file sizes
    size_t sizen = get_file_size_MPI(filename) - in.off_file;  //qassert(sizen == string_to_size(in.total_size));
    if(sizen != string_to_size(in.total_size)){abort_r("FILE size not match with head !\n");}

    in.Vsize = size_t(in.nx)*in.ny*in.nz*size_t(in.bfac*2);
    size_t Fsize = size_t(n1)*(12)*in.Vsize*in.bsize;  //qassert(Fsize <= string_to_size(in.total_size));
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
inline int check_qlat_eigen_file_type(const char *filename, io_vec &io_use,int n1, inputpara& in){
  int type = check_eigen_qlat(filename, n1, in);
  if(type != -1){return type;}

  bool check = true;
  type = check_Eigen_file_type(filename, io_use, n1, check);
  return type;
}
////===Check eigen system type


template <class T, int bfac>
void load_qlat_noises(const char *filename, std::vector<qlat::FieldM<T, bfac> > &noises, bool read=true, bool single_file=true, const std::string& INFO_LIST = std::string("NONE"), int n0=0,int n1=-1)
{
  TIMERB("load/save qlat noises");
  std::string VECS_TYPE = std::string("Noise_Vectors");
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST, n0, n1);
}

template <class T, int bfac>
void save_qlat_noises(const char *filename, std::vector<qlat::FieldM<T, bfac> > &noises, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_noises(filename, noises, false, single_file, INFO_LIST);
}

template <class T, int bfac>
void load_qlat_noise(const char *filename, qlat::FieldM<T, bfac> &noise, bool read=true, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  std::vector<qlat::FieldM<T, bfac> > noises;
  ///noises[0].init(noise.geo());
  if(read == false){noises.resize(1);noises[0] = noise;}
  load_qlat_noises(filename, noises, read, single_file, INFO_LIST, 0, 1);
  if(read == true ){noise = noises[0];}
}

template <class T, int bfac>
void save_qlat_noise(const char *filename, qlat::FieldM<T, bfac> &noise, bool single_file=true, const std::string& INFO_LIST = std::string("NONE")){
  load_qlat_noise(filename, noise, false, single_file, INFO_LIST);
}

//////Assume memory allocated already
template<class T, typename Ty>
void copy_noise_to_prop(qlat::FieldM<T, 12*12>& noise, Propagator4dT<Ty>& prop, int dir=1)
{
  TIMERB("copy_noise_to_prop");
  if(dir == 1){prop.init(noise.geo());}
  if(dir == 0){noise.init(prop.geo());}
  T* noi = (T*) qlat::get_data(noise).data();
  #pragma omp parallel for 
  for (long index = 0; index < prop.geo().local_volume(); ++index)
  {
    ///qlat::WilsonMatrixT<T>& src =  prop.get_elem(index);
    T* src   = (T*) &noi[index*12*12];
    Ty* res  = &prop.get_elem(index)(0,0);
    
    for(int d0=0;d0<12;d0++)
    {
      for(int d1=0;d1<12;d1++)
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
template<class T, typename Ty>
void copy_noises_to_prop(std::vector<qlat::FieldM<T, 12*12> >& noises, Propagator4dT<Ty>& prop, int dir=1)
{
  TIMERB("copy_noises_to_prop");
  if(dir == 1){prop.init(noises[0].geo());}
  if(dir == 0){noises.resize(0);noises.resize(1);noises[0].init(prop.geo());}
  copy_noise_to_prop(noises[0], prop, dir);

  //T* noi = (T*) qlat::get_data(noises[0]).data();
  //#pragma omp parallel for 
  //for (long index = 0; index < prop.geo().local_volume(); ++index)
  //{
  //  ///qlat::WilsonMatrixT<T>& src =  prop.get_elem(index);
  //  T* src   = (T*) &noi[index*12*12];
  //  Ty* res  = &prop.get_elem(index)(0,0);
  //  
  //  for(int d0=0;d0<12;d0++)
  //  {
  //    for(int d1=0;d1<12;d1++)
  //    {
  //      //////copy to prop
  //      if(dir==0){src[d1*12+d0] = res[d0*12+d1];}
  //      //////copy to buf
  //      if(dir==1){res[d0*12+d1] = src[d1*12+d0];}

  //    }

  //  }
  //}
}

template <class Ty>
void load_qlat_prop(const char *filename, Propagator4dT<Ty>& prop, bool read=true, bool single_file=true){
  std::string VECS_TYPE = std::string("Propagator");
  std::string INFO_LIST  = std::string("src 12, sink 12, zyxt, R/I");
  std::vector<qlat::FieldM<Ty, 12*12> > noises;
  if(read == false){copy_noises_to_prop(noises, prop, 0);}
  load_qlat_noisesT(filename, noises, read, single_file, VECS_TYPE, INFO_LIST);
  if(read == true){copy_noises_to_prop(noises, prop, 1);}
  
}

template <class Ty >
void save_qlat_prop(const char *filename,Propagator4dT<Ty >& prop, bool single_file=true){
  load_qlat_prop(filename, prop, false, single_file);
}

/////nvec needed for checksum
inline FILE* open_eigensystem_file(const char *filename, int nini, int nvec, bool read, io_vec& io_use, inputpara& in, int save_type = 3)
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

  if(in.file_type == 2 or in.file_type == 3){
    if(read == false){in.read_geo(io_use.geop);}
    int bfac_eigen = 12;
    std::string VECS_TYPE = std::string("Eigen_system_nvec.12.tzyx.R/I");
    std::string INFO_LIST = std::string("NONE");
    bool rotate_bfac = true;
    open_file_qlat_noisesT(filename, bfac_eigen, in, read, in.single_file, nvec, VECS_TYPE, INFO_LIST, rotate_bfac);

    if(io_use.end_of_file != 0){abort_r("Please close file first! \n");}

    io_use.end_of_file = in.end_of_file;
    io_use.do_checksum = in.do_checksum;
  }

  FILE* file=NULL;
  if(read==true )file = io_use.io_read(in.filename.c_str(),"rb");
  if(read==false)file = io_use.io_read(in.filename.c_str(),"wb");

  //////shift the file to nini position
  {
    size_t bsize = sizeof(double);int bfac = 12;
    if(!in.single_file){bsize=sizeof(float) ;}
    if( in.single_file){bsize=sizeof(double);}
    size_t Vsize = size_t(in.nx)*in.ny*in.nz*in.nt*size_t(bfac*2);
    size_t off_file = size_t(nini)*Vsize*bsize;
    if(in.file_type == 2 or in.file_type == 3){off_file += in.off_file;}
    io_use.io_off(file, off_file, false);
  }

  return file;
}


template <class Ty>
void load_eigensystem_vecs(FILE* file, std::vector<qlat::FieldM<Ty, 12> > &noises, io_vec& io_use, inputpara& in, int n0=0, int n1=-1)
{

  if(in.file_type == 0 or in.file_type == 1)
  {
    if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
    int n_vec = n1-n0;
    ////if(in.file_type == 0){read_single = false;}
    ////if(in.file_type == 1){read_single = true ;}

    bool check = true;
    ////load_gwu_eigen(file, noises, io_use, n0, n1, check, in.read, in.single_file );
    if(in.read == true){
    if(noises.size() < (LInt) n_vec)
    {
      noises.resize(0);
      noises.resize(n_vec);
      for(int iv=0;iv<n_vec;iv++){noises[iv].init(io_use.geop);}
    }}

    if(sizeof(Ty) == 2*sizeof(double)){
      std::vector<double* > respD;respD.resize(n_vec);
      for(int iv=0;iv<n_vec;iv++){respD[iv] = (double*) qlat::get_data(noises[iv]).data();}
      load_gwu_eigen(file, respD, io_use,n0,n1,check, in.read, in.single_file );
    }
    if(sizeof(Ty) == 2*sizeof(float) ){
      std::vector<float*  > respF;respF.resize(n_vec);
      for(int iv=0;iv<n_vec;iv++){respF[iv] = (float*) qlat::get_data(noises[iv]).data();}
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



}

#endif
