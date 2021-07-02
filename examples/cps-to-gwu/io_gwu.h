// io_gwu.h
// Gen Wang
// Jan. 2021

#ifndef IO_GWU_H
#define IO_GWU_H
#pragma once

#include <qlat/qlat.h>
#include "general_funs.h"

//#include "utils_float_type.h"

//#include "../genwang/general.h"
//using namespace qcd;

////Target, read double prop, double eigensystem
////read single prop, single eigensystem
////all in cps base,
////prop --> 12 *vol * d * c * complex, eigensytem, n * volv * d * c * complex
////Define gammas in PS, chiral, cps, base , ga
//

namespace qlat
{

inline bool is_big_endian_gwu(void){
  int num = 1;
  if(*(char *)&num == 1)
  {
    return false;
    //printf("\nLittle-Endian\n");
  }
  else
  {
    return true;
    //printf("Big-Endian\n");
  }
}

inline void swapbytes(void *_object, size_t _size)
{
   unsigned char *start, *end;
   for ( start = (unsigned char *)_object, end = start + _size - 1; start < end; ++start, --end )
   {   
      unsigned char swap = *start;
      *start = *end;
      *end = swap;
   }   
}

inline void switchendian(void *buffer,size_t length,int dsize)
{
  #pragma omp parallel for 
  for(size_t i=0;i<length;i++)
  {
    char* pos = ((char*)buffer) + i*dsize;
    swapbytes(pos,dsize);
  }

}

class io_gwu
{
public:

  ///fft_desc_basic* fd_basic;
  std::vector<int > node_ioL;
  int ionum;
  int nx,ny,nz,nt;
  int Nx,Ny,Nz,Nt;
  size_t vol;
  int rank,noden,Nmpi;
  std::vector<int > nv,Nv;
  const qlat::Geometry* geop;


  //////each noden has the memory order to t,z,y,x
  std::vector<std::vector<LInt > > map_Nitoi;
  ///void get_pos(int i);
  //FILE* io_read(const char *filename,const char* mode);
  //void io_close(FILE *file);
  ////void io_off(FILE *file,size_t off);

  FILE* io_read(const char *filename,const char* mode){
    if(node_ioL[rank]>=0){
      return fopen(filename, mode);
    }else{return NULL;}
  }
  
  void io_close(FILE *file){
    if(node_ioL[rank]>=0){fclose(file);}
    file = NULL;
  }

  void io_off(FILE *file,size_t off){
    if(node_ioL[rank]>=0){
    fseek(file , off , SEEK_CUR );
    }
  }


  ////fd = &fds;
  io_gwu(const qlat::Geometry& geo,int ionum_or,bool checkorder=true){
  /////x,y,z,t
  geop = &geo;
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
  if(ionum_or > 0 and ionum_or < Nmpi){ionum=ionum_or;}else{ionum = Nmpi;}

  for(int ni=0;ni<Nmpi;ni++){node_ioL[ni]=-1;}
  int off = Nmpi/ionum;
  int countN = 0;
  for(int ni=0;ni<Nmpi;ni++){
    if(ni%off==0 and countN < ionum){node_ioL[ni]=countN;countN++;}
  }
  }

  ~io_gwu(){
  map_Nitoi.resize(0);
  node_ioL.resize(0);
  }


};

inline void send_vec_kentucky(char* src,char* res,int dsize,int gN,const io_gwu& io, bool read=true)
{
  const std::vector<int >& node_ioL = io.node_ioL;
  int rank = io.rank;size_t vol = io.vol;int noden = io.noden;int Nmpi=io.Nmpi;

  int Nx = io.Nx;

  int size_c = gN*noden*dsize;

  char* tmp=NULL;
  if(node_ioL[rank]>=0){tmp = new char[gN*io.vol*dsize];}

  ///position p;
  if(read==true)
  {
  //TIMER("IO sort mem time");
  if(node_ioL[rank]>=0)
  for(int ri=0;ri<io.Nmpi;ri++)for(int gi=0;gi<gN;gi++){
      #pragma omp parallel for
      for(size_t isp=0;isp<size_t(noden/Nx);isp++){
        char* pres=&tmp[(ri*gN+gi)*noden*dsize];
        char* psrc=&src[gi*vol*dsize];
        size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
        memcpy(&pres[isp*Nx*dsize],&psrc[offv*dsize],Nx*dsize);
      }
  }
  }

  {
  /////TIMER("IO MPI time");
  /////for(int io=0;io<node_ioL.size();io++)if(node_ioL[io]>=0)
  /////{
  /////  MPI_Scatter(tmp,size_c,MPI_CHAR,&res[node_ioL[io]*size_c],size_c,MPI_CHAR,io,get_comm());
  /////  ////if(rank==io){memcpy(&res[node_ioL[io]*size_c],tmp,size_c);}
  /////}

  std::vector<int > currsend,currrecv,currspls,currrpls;
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
      ////if(rank == n)for(int ni=0;ni<Nmpi;ni++)currsend[ni] = size_c;
      currrecv[n] = size_c;
      currrpls[n] = size_c*node_ioL[n];
    }
  }

  //if(qlat::get_num_node() != 1)
  //{
  if(read==true)
  {MPI_Alltoallv(tmp,(int*) &currsend[0],(int*) &currspls[0], MPI_CHAR,
                 res,(int*) &currrecv[0],(int*) &currrpls[0], MPI_CHAR, get_comm());}
  if(read==false)
  {MPI_Alltoallv(res,(int*) &currrecv[0],(int*) &currrpls[0], MPI_CHAR,
                 tmp,(int*) &currsend[0],(int*) &currspls[0], MPI_CHAR, get_comm());}
  //}
  //else{
  //  //memcpy(res,tmp, size_c);
  //  #pragma omp parallel for
  //  for(size_t isp=0;isp<size_c;isp++){res[isp] = tmp[isp];}
  //}

  }

  if(read == false)
  {
  //TIMER("IO sort mem time");
  if(node_ioL[rank]>=0)
  for(int ri=0;ri<io.Nmpi;ri++)for(int gi=0;gi<gN;gi++){
    #pragma omp parallel for
    for(size_t isp=0;isp<size_t(noden/Nx);isp++){
      char* pres=&tmp[(ri*gN+gi)*noden*dsize];
      char* psrc=&src[gi*vol*dsize];
      size_t offv = io.map_Nitoi[ri][isp*Nx+ 0];
      memcpy(&psrc[offv*dsize],&pres[isp*Nx*dsize],Nx*dsize);
    }
  }
  }

  if(node_ioL[rank]>=0){delete []tmp;}

}

inline void read_kentucky_vector(FILE *file,char* props,int Nvec,const io_gwu& io,bool Rendian=false, int dsize=8, bool single_file=false,  int gN=1, bool read = true)
{
  //TIMER("new io load");
  timeval tm0,tm1,tm2,tm3;
  gettimeofday(&tm0, NULL);
  double mpi_t = 0.0;

  const std::vector<int >& node_ioL = io.node_ioL;
  int rank = io.rank;size_t vol = io.vol;int ionum=io.ionum;int noden=io.noden;

  char* buf=NULL;if(node_ioL[rank]>=0)buf = new char[gN*vol*dsize];
  char res[ionum*noden*gN*dsize];
  size_t sizec = vol*dsize;

  size_t off_file = 0;
  if(node_ioL[rank]>=0){off_file = ftell(file);}

  int curr_v = 0;
  for(int readi=0;readi<Nvec;readi++){

    /////From res to buf
    if(read==false)
    {
    //TIMER("IO copy mem time");
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
    //TIMER("IO disk");
    if(node_ioL[rank]>=0)
    if(curr_v + node_ioL[rank]*gN < Nvec)
    {
      int offv = curr_v + node_ioL[rank]*gN;
      fseek ( file , off_file + offv*sizec , SEEK_SET );

      int rN = gN;if(offv+rN>=Nvec){rN = Nvec-offv;}

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
        if(read==true ) fread(&buf[0], 1, rN*sizec, file);
        if(read==false)fwrite(&buf[0], 1, rN*sizec, file);
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
    //TIMER("IO copy mem time");
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
  if(node_ioL[rank]>=0){delete []buf;}

  gettimeofday(&tm1, NULL);
  double time0 = tm1.tv_sec - tm0.tv_sec;
  time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
  double total_R = vol*Nvec*dsize/(1024*1024*1024.0);

  if(read==true) print0("READ_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);
  if(read==false)print0("WRIT_TIMING: total %.3e s , MPI %.3e s, disk %.3f GB/s, eff %.3f GB/s \n",time0,mpi_t,total_R/(time0-mpi_t),total_R/time0);

}

//inline void save_kentucky_vector(FILE *file,char* props,int Nvec,const io_gwu& io,bool Rendian=false, int dsize=8, bool single_file=false,  int gN=1)
//{
//  //TIMER("new io load");
//  timeval tm0,tm1;gettimeofday(&tm0, NULL);
//
//  const std::vector<int >& node_ioL = io.node_ioL;
//  int rank = io.rank;size_t vol = io.vol;int ionum=io.ionum;int noden=io.noden;
//
//  char* buf=NULL;if(node_ioL[rank]>=0)buf = new char[gN*vol*dsize];
//  char res[ionum*noden*gN*dsize];
//  size_t sizec = vol*dsize;
//
//  size_t off_file = 0;
//  if(node_ioL[rank]>=0){off_file = ftell(file);}
//
//  int curr_v = 0;
//  for(int readi=0;readi<Nvec;readi++){
//    {
//    //TIMER("IO disk");
//    if(node_ioL[rank]>=0)
//    if(curr_v + node_ioL[rank]*gN < Nvec)
//    {
//      int offv = curr_v + node_ioL[rank]*gN;
//      fseek ( file , off_file + offv*sizec , SEEK_SET );
//
//      int rN = gN;if(offv+rN>=Nvec){rN = Nvec-offv;}
//      if(rN>0){fread(&buf[0], 1, rN*sizec, file);}
//      //fwrite(fread(&buf[0], 1, rN*sizec, file));
//
//      int dsize_single = sizeof(double);int fac = dsize/sizeof(double);
//      if(single_file == true){dsize_single = sizeof(float);fac = dsize/sizeof(float);}
//
//      //////Read Prop and link
//      if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
//      //////Read source noise
//      if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], gN*vol*fac,dsize_single);
//    }
//    }
//
//    ////////Reordersrc buf
//    ///if(gN!=1){if(node_ioL[rank]>=0)reorder_civ((char*)&buf[0],(char*)&buf[0],1,gN,vol,0 ,sizeof(double));}
//    /////
//      
//    send_vec_kentucky((char*) &buf[0],(char*) &res[0], dsize,gN, io);
//
//    {
//    //TIMER("IO copy mem time");
//    for(int io=0;io<ionum;io++)
//    for(int gi=0;gi<gN;gi++)if(curr_v + io*gN + gi<Nvec){
//      int offv = curr_v + io*gN + gi;
//
//      memcpy(&props[(offv)*noden*dsize + 0],&res[(io*gN+gi)*noden*dsize + 0],noden*dsize);
//    }
//    }
//    curr_v += ionum*gN;
//    if(curr_v >= Nvec)break;
//
//  //fseek ( pFile , 9 , SEEK_SET );
//  //fseek ( pFile , 9 , SEEK_CUR );
//  //fseek ( pFile , 9 , SEEK_END );
//  }
//
//  
//
//  if(node_ioL[rank]>=0){fseek ( file , off_file + Nvec*sizec , SEEK_SET );}
//  if(node_ioL[rank]>=0){delete []buf;}
//
//  gettimeofday(&tm1, NULL);
//  double time0 = tm1.tv_sec - tm0.tv_sec;
//  time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;
//  double total_R = vol*Nvec*dsize/(1024*1024*1024.0);
//  print0("READ_TIMING: total %.3e s  eff %.3f GB/s \n",time0,total_R/time0);
//
//}

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
    for(int i=0;i<4;i++)fgets(sTemp, 300, filer);

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
  //TIMER("Rotate gwu file vec");
  ////resize res vector
  ///res.resize(src.size());
  //size_t noden = geo.node_site[0]*geo.node_site[1]*geo.node_site[2]*geo.node_site[3];

  ///std::vector<int > nv,Nv;nv.resize(4);Nv.resize(4);
  ///for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  ///Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];
  ///int n_vec = src.size()/(noden*12*2);

  ////if(single_file == false and read=false){abort_r("Cannot do double file rotation now!\n");}

  if(single_file == true){
    //#pragma omp parallel for
    //for(size_t isp=0;isp<src.size();isp++)res[isp] = src[isp];
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

        ///res[(dc1*noden + isp)*12*2 + dc0*2 + 0] = src[(dc1*12*2 + 0*12 + dc0)*noden + isp];
        ///res[(dc1*noden + isp)*12*2 + dc0*2 + 1] = src[(dc1*12*2 + 1*12 + dc0)*noden + isp];
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
  //TIMER("Rotate gwu file vec");
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

///////TO DO, Need use omp here
////double get_norm_vec_qlat(qlat::Complex *src, size_t noden){
////  double res = 0.0;
////  /////need sum 12 first to reduce float sum error
////  #pragma omp parallel for reduction(+: res)
////  for(size_t isp=0;isp<noden;isp++)
////  {
////    double a = 0.0;
////    for(unsigned int dc=0;dc<12;dc++){
////       a += qnorm(src[isp*12+dc]);
////    }
////    res += a;
////  }
////  sum_all_size(&res,1);
////  return res;
////}

//template<typename Ty>
//double get_norm_vec_qlat(Ty *src, size_t noden){
//  double res = 0.0;
//  ///////need sum 12 first to reduce float sum error
//  //#pragma omp parallel for reduction(+: res)
//  for(size_t isp=0;isp<noden;isp++)
//  {
//    double a = 0.0;
//    for(unsigned int dc=0;dc<12;dc++){
//       a += qnorm(src[isp*12+dc]);
//    }
//    res += a;
//  }
//  sum_all_size(&res,1);
//  return res;
//}

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

inline int test_single(const char *filename,io_gwu &io_use,int iv=0){
  
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
  io_use.io_off(file,iv*Fsize*2);
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
  io_use.io_off(file,iv*Fsize);

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

inline bool check_Eigen_file_size(const char *filename,io_gwu &io_use,int n1,bool check){

  size_t noden = io_use.noden;
  size_t sizen = get_file_size_MPI(filename);

  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  if(sizen <  Fsize){print0(filename);abort_r("Eigen file size too small. \n");}
  if(sizen <2*Fsize){
    if(n1 != 1){abort_r("Eigen file size too small 1. \n");}
    /////read with single
    return true;
  }

  ////Default use single
  if(check == false){
    if(sizen < n1*Fsize){abort_r("Eigen file size too small 2. \n");}
    return true;
  }

  int tmp = test_single(filename,io_use,0);
  if(tmp == -1){abort_r("Eigen system file norm not correct.\n");}

  if(tmp == 1){if(sizen <  n1*Fsize){abort_r("Eigen file size too small 2. \n");}return true ;}
  if(tmp == 0){if(sizen <2*n1*Fsize){abort_r("Eigen file size too small 3. \n");}return false;}

  return false;
}

/////Read the data into the point resp which have memory allocated already
/////check --> if false, abort if file is double,
/////read --> flag for writing and reading, can only write single currently
template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty* > resp,io_gwu &io_use,int n0,int n1,
  bool check=true, bool read=true)
{
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(resp.size()!= size_t(n1-n0)){abort_r("Final point size wrong!\n");}
  ////if(n0!=0 and read==false){abort_r("Write vectors with off not zero!");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);
  bool single = true;
  ////read==false, write only single vectors
  if(read==true)single = check_Eigen_file_size(filename,io_use,n1,check);
  FILE* file;

  if(single == true){
    /////Load single precision
      
    int count = 0;int off = io_use.ionum;
    std::vector<float > prop_E;
    int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    if(read==true )file = io_use.io_read(filename,"rb");
    if(read==false)file = io_use.io_read(filename,"wb");
    io_use.io_off(file,n0*Fsize);

    for(int iv=0;iv<n_vec;iv++){
      //TIMER("New load Eigen");
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

      read_kentucky_vector(file,(char*) &prop_E[0], ri,io_use, true, 6*2*sizeof(double), true, 1, read);
      fflush_MPI();

      if(read==true){
        //rotate_qlat_to_gwu(prop_E, &src_new.vec[0],geo, true);
        /////single precision eigen vector in milc base
        rotate_gwu_vec_file(&prop_E[0],ri,noden, true);
        ////Do not rotate source and sink
        gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, false, false);
        for(int ic=0;ic<ri;ic++)
        #pragma omp parallel for
        for(size_t isp=0;isp<noden*12*2;isp++){
          resp[count + ic][isp] = prop_E[ic*noden*12*2 + isp];}
      }

      count += ri;if(count >= n_vec){break;}
    }
    io_use.io_close(file);

    /////rotate_qlat_to_gwu(prop_E,&src_new.vec[0],geo);
    /////////From ky to milc
    /////for(int iv=0;iv<n_vec;iv++)ga.ga[4][0].multiply(*src_new.vec[iv],*src_new.vec[iv]);
  }else{
    /////Load double precision

    int count = 0;int off = io_use.ionum;
    std::vector<double > prop_E;
    int n_vec = n1-n0;
    prop_E.resize(off*12*noden*2);

    file = io_use.io_read(filename,"rb");
    io_use.io_off(file,n0*Fsize*2);
      
    for(int iv=0;iv<n_vec;iv++){
      //TIMER("New load Eigen");
      int ri = off;if(count + off > n_vec){ri = n_vec - count;}
      //int offE = count*12*noden*2;
      read_kentucky_vector(file,(char*) &prop_E[0], ri*12*2,io_use, false, sizeof(double), false, 12*2);
      fflush_MPI();

      /////double precision eigen vector in ps base
      rotate_gwu_vec_file(&prop_E[0],ri,noden, false);
      ////Do not rotate source, rotate sink
      gwu_to_cps_rotation_vec(&prop_E[0],ri,noden, false, true, true);

      for(int ic=0;ic<ri;ic++)
      #pragma omp parallel for
      for(size_t isp=0;isp<noden*12*2;isp++){
        resp[count + ic][isp] = prop_E[ic*noden*12*2 + isp];
      }

      count += ri;if(count >= n_vec){break;}
    }

    io_use.io_close(file);

    /////////Copy back to gwu mem order
    /////rotate_qlat_to_gwu(prop_E,&src_new.vec[0],geo);
    /////////From cps to milc
    /////for(int iv=0;iv<n_vec;iv++)ga.ga[4][0].multiply(*src_new.vec[iv],*src_new.vec[iv]);
  }



}

//////Ty should be float or double
template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<Ty > &res,io_gwu &io_use,int n0,int n1,bool check=true,bool read=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  int n_vec = n1-n0;
  size_t noden = io_use.noden;
  res.resize(n_vec*noden*12*2);

  std::vector<Ty* > resp;resp.resize(n_vec);
  for(int iv=0;iv<n_vec;iv++){resp[iv] = &res[iv*noden*12*2];}

  load_gwu_eigen(filename,resp,io_use,n0,n1,check,read);

}


template<typename Ty>
void load_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &eigen,io_gwu &io_use,int n0,int n1,bool check=true,bool read=true){
  if(n1<n0 or n0<0){abort_r("Read number of eigen should be larger than 1. \n");}
  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  int n_vec = n1-n0;
  //size_t noden = io_use.noden;
  //res.resize(n_vec*noden*12*2);

  //std::vector<float*  > resf;resf.resize(n_vec);
  ////for(int iv=0;iv<n_vec;iv++){resp[iv] = &res[iv*noden*12*2];}
  //if(sizeof(Ty) == sizeof(float))
  //if(sizeof(Ty) == sizeof(float))

  //std::vector<qlat::FermionField4dT<Ty> > eigen;
  if(read==true){
    eigen.resize(0);
    eigen.resize(n_vec);
    for(int iv=0;iv<n_vec;iv++)eigen[iv].init(*io_use.geop);
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
void save_gwu_eigen(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &eigen,io_gwu &io_use,int n0,int n1,bool check=true){
  load_gwu_eigen(filename,eigen,io_use,n0,n1,check,false);
}

inline void load_gwu_eigen(const char *filename,EigenM &Mvec,io_gwu &io_use,int n0,int n1,bool check=true,bool read=true)
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

inline void save_gwu_eigen(const char *filename,EigenM &Mvec,io_gwu &io_use,int n0,int n1,bool check=true)
{
  load_gwu_eigen(filename, Mvec,io_use,n0,n1,check,false);
}


template<typename Ty>
void load_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &prop,io_gwu &io_use,bool read=true){

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*12*2)*sizeof(float);

  FILE* file;bool single = true;
  if(read==true){
  size_t sizen = get_file_size_MPI(filename);

  if(sizen != 2*Fsize*12 and sizen != Fsize*12){abort_r("prop size wrong! \n");}
  prop.resize(0);
  prop.resize(12);
  for(int iv=0;iv<12;iv++)prop[iv].init(*io_use.geop);
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
    read_kentucky_vector(file,(char*) &prop_qlat[0], 12*2*12,io_use, false, sizeof(double), false,12*2);
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
void save_gwu_prop(const char *filename,std::vector<qlat::FermionField4dT<Ty> > &prop,io_gwu &io_use){
  load_gwu_prop(filename,prop,io_use,false);
}

template<typename Ty>
void prop4d_to_Fermion(Propagator4d& prop,std::vector<qlat::FermionField4dT<Ty > > &buf, int dir=1){

  ////if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}
  if(dir==1){buf.resize(0);buf.resize(12);for(int iv=0;iv<12;iv++){buf[iv].init(prop.geo());qlat::set_zero(buf[iv]);}}
  if(dir==0){prop.init(buf[0].geo());}

  for (long index = 0; index < prop.geo().local_volume(); ++index)
  {
    qlat::WilsonMatrix& src =  prop.get_elem(index);
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

inline void save_gwu_prop(const char *filename,Propagator4d& prop){
  io_gwu io_use(prop.geo(),32);
  std::vector<qlat::FermionField4dT<qlat::ComplexF> > prop_qlat;
  prop4d_to_Fermion(prop,prop_qlat, 1);
  save_gwu_prop(filename,prop_qlat,io_use);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

inline void save_gwu_prop(std::string &filename,Propagator4d& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  save_gwu_prop(tem,prop);
}

inline void load_gwu_prop(const char *filename,Propagator4d& prop){
  io_gwu io_use(prop.geo(),32);
  std::vector<qlat::FermionField4dT<qlat::Complex > > prop_qlat;
  load_gwu_prop(filename,prop_qlat,io_use);
  prop4d_to_Fermion(prop,prop_qlat, 0);
  ///////load_gwu_prop(filename,prop,io_use,false);
}

inline void load_gwu_prop(std::string &filename,Propagator4d& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_prop(tem,prop);
}



////void save_gwu_link(const char *filename,GaugeField &gf,io_gwu &io_use,bool read=true){
////  load_gwu_link(filename,gf,io_use,false);
////}

inline void load_gwu_link(const char *filename,GaugeField &gf){
  io_gwu io_use(gf.geo(),8);
  //if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float ))
  //{abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(4*9*noden*2)*sizeof(double);

  FILE* file;
  //if(read==true)
  {
  size_t sizen = get_file_size_MPI(filename);
  if(sizen != Fsize){abort_r("Link size wrong! \n");}
  /////gf.init(*io_use.geop);
  }

  std::vector<double > link_qlat;
  link_qlat.resize(4*9*noden*2);

  file = io_use.io_read(filename,"rb");
  read_kentucky_vector(file,(char*) &link_qlat[0], 4*9*2,io_use, false, sizeof(double), false,9*2);
  io_use.io_close(file);
    
  //////double precision eigen vector in ps base
  ///rotate_gwu_vec_file(&prop_qlat[0], 12,noden, false);
  /////Do not rotate source, 
  ///gwu_to_cps_rotation_vec(&link_qlat[0], 12,noden, true, true,true);
  ///4 dire --> c0 -- > c1

  ////Geometry &geo = *io_use.geop;
  ////May need to check gf is double prec
  for (size_t index = 0; index < noden; ++index)
  {
    double* res = (double*) qlat::get_data(gf.get_elems(index)).data();

    ////for(int dir=0;dir<4*9*2;dir++)res[dir] = link_qlat[dir*noden + index];
    for(int dir=0;dir<4;dir++)
    for(int c0=0;c0<3;c0++)
    for(int c1=0;c1<3;c1++)
    for(int im=0;im<2;im++){
      int dir0 = ((dir*3+c0)*3+c1)*2+im;
      int dir1 = ((dir*3+c1)*3+c0)*2+im;
      res[dir0] = link_qlat[dir1*noden + index];
    }
  }
}

inline void load_gwu_link(std::string &filename,GaugeField &gf){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_link(tem,gf);
}


template<typename Ty>
void load_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi,io_gwu &io_use,bool read=true){

  if(sizeof(Ty) != 2*sizeof(double ) and sizeof(Ty) != 2*sizeof(float )){abort_r("Cannot understand the input format! \n");}

  size_t noden = io_use.noden;
  size_t Fsize = io_use.Nmpi*(noden*2)*sizeof(float);

  std::vector<double > prop_noi;
  prop_noi.resize(noden*2);

  if(read==true){
  size_t sizen = get_file_size_MPI(filename);
  if(sizen != 2*Fsize){abort_r("noise size wrong! \n");}

  qlat::Geometry geo = *io_use.geop;geo.multiplicity=1;
  noi.init(geo);
  }
  if(read==false){
    std::complex<double> *src = (std::complex<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    for(size_t isp=0;isp<noden;isp++)src[isp] = std::complex<double>(res[isp].real(),res[isp].imag());
    reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 1,sizeof(double));
  }

  FILE* file;

  if(read==true )file = io_use.io_read(filename,"rb");
  if(read==false)file = io_use.io_read(filename,"wb");
  read_kentucky_vector(file,(char*) &prop_noi[0], 2,io_use, true, sizeof(double), false, 1,read);
  io_use.io_close(file);

  /////Move complex to inner side
  if(read==true){
    reorder_civ((char*) &prop_noi[0],(char*) &prop_noi[0], 1, 2, noden, 0,sizeof(double));
    std::complex<double> *src = (std::complex<double>*) &prop_noi[0];
    Ty* res = (Ty*) qlat::get_data(noi).data();
    for(size_t isp=0;isp<noden;isp++)res[isp]= Ty(src[isp].real(),src[isp].imag());
  }

}

template<typename Ty>
void save_gwu_noi(const char *filename,qlat::FieldM<Ty,1> &noi,io_gwu &io_use){
  load_gwu_noi(filename,noi,io_use,false);
}

inline void save_gwu_noiP(const char *filename,Propagator4d& prop){
  io_gwu io_use(prop.geo(),32);
  qlat::FieldM<qlat::Complex,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  
  size_t noden = io_use.noden;
  for (size_t index = 0; index < noden; ++index)
  {
    qlat::WilsonMatrix& src =  prop.get_elem(index);
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

inline void save_gwu_noiP(std::string &filename,Propagator4d& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  save_gwu_noiP(tem,prop);
}

inline void load_gwu_noiP(const char *filename,Propagator4d& prop){
  io_gwu io_use(prop.geo(),32);
  qlat::FieldM<qlat::Complex,1> noi;
  noi.init(prop.geo());
  qlat::set_zero(noi);
  load_gwu_noi(filename,noi ,io_use);
  prop.init(noi.geo());
  
  size_t noden = io_use.noden;
  for (size_t index = 0; index < noden; ++index)
  {
    ////Assumed double
    ////double* res = (double*) qlat::get_data(prop.get_elems(index)).data();
    qlat::WilsonMatrix& res =  prop.get_elem(index);
    for(int d0=0;d0<12;d0++){res(d0,d0) = noi.get_elem(index);}

    ////src = (double*) qlat::get_data(noi.get_elems(index)).data();

    /////for(int d1=0;d1<12;d1++)
    /////for(int d0=0;d0<12;d0++)
    /////{
    /////  sum += std::fabs(src[(d1*12 + d0)*2+0]);
    /////  sum += std::fabs(src[(d1*12 + d0)*2+1]);
    /////}
    //{
    //  res(d0,d0) = noi.get_elems(index);
    //  //res[(d0*12+d0)*2+0] = src[0];
    //  //res[(d0*12+d0)*2+1] = src[1];
    //}
    ///qlat::Complex phase = qlat::Complex(src[0],src[1]);
    ///if(sum >1e-8)noi.get_elem(index) = 1.0*phase;
  }
  ///////load_gwu_prop(filename,prop,io_use,false);
}

inline void load_gwu_noiP(std::string &filename,Propagator4d& prop){
  char tem[500];
  sprintf(tem,filename.c_str());
  load_gwu_noiP(tem,prop);
}

inline void write_data(std::vector<double > dat,const char *filename, bool read=false){
  bool Rendian = true;
  int size = 0;

  if(read==false)size = dat.size();
  if(read==true){size_t sizeF = get_file_size_o(filename);size = sizeF/sizeof(double);dat.resize(size);}

  ////Set buf with size
  char* buf=NULL;
  buf = new char[size*sizeof(double)];

  ////Open file
  FILE* file = NULL;
  if(read==false)file = fopen(filename, "wb");
  if(read==true )file = fopen(filename, "rb");

  /////Switch endian of the file write
  if(read==false){
    memcpy(&buf[0],&dat[0],size*sizeof(double));
    if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size ,sizeof(double));
  }

  if(read==false)fwrite(&buf[0], 1, size*sizeof(double), file);
  if(read==true ) fread(&buf[0], 1, size*sizeof(double), file);

  /////Switch endian of the file write
  if(read==true ){
    if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size ,sizeof(double));
    memcpy(&dat[0],&buf[0],size*sizeof(double));
  }

  fclose(file);
  delete []buf;

}

inline void read_data(std::vector<double > dat,const char *filename){
  write_data(dat, filename, true);
}



}

#endif
