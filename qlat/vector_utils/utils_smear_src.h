// utils_smear_src.h
// Gen Wang
// Jul. 2021

#ifndef UTILS_SMEAR_SRC_H
#define UTILS_SMEAR_SRC_H
#pragma once

#include "general_funs.h"
#include <qlat/qcd-utils.h>
#include <qlat/qcd-prop.h>
//#include <openacc.h>
//#include "mpi-ext.h" /* Needed for CUDA-aware check */
 
namespace qlat{
template <class T>
void rotate_prop(Propagator4dT<T>& prop, int dir = 0)
{
  TIMER("Rotate color prop");
  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();

  qacc_for(isp, long(Nvol), {
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem(isp);
    qlat::WilsonMatrixT<T>  v1 =  prop.get_elem(isp);

    for(int c0 = 0;c0 < 3; c0++)
    for(int d0 = 0;d0 < 4; d0++)
    for(int c1 = 0;c1 < 3; c1++)
    for(int d1 = 0;d1 < 4; d1++)
    {
      //int a0 = (d0*3+c0)*12+ d1*3+c1;
      //LInt off = ((d0*4+d1)*3 + c0)*3 + c1;
      LInt off = (c0*3+c1)*16+d0*4 + d1;
      int a0 = off/12;
      int a1 = off%12;
      if(dir == 0)v0(a0, a1) = v1(d0*3+c0, d1*3+c1);
      if(dir == 1)v0(d0*3+c0, d1*3+c1) = v1(a0, a1);

      //if(dir==0)pE[off*Nvol + isp] = v0(d0*3 + c0, d1*3 + c1);
      //if(dir==1)v0(d0*3 + c0, d1*3 + c1) = pE[off*Nvol + isp];
    }
  });
}

template <class T>
__global__ void smear_global3(T* pres, T* psrc, const T* gf, double bw, double norm, long Nvol, const Geometry& geo,const Geometry& geo_v)
{

  __shared__ T ls[6*3*3 + 1];
  __shared__ T ps[6*3*3*16 + 1];
  /////__shared__ T ws[  3*3*16 + 1];

  //__shared__ T ls[3*3*6 + 1];
  //__shared__ T ps[16][3*6*3 + 1];
  //__shared__ T ws[16][3*3 + 1];

  ////////======c, mxn, a, mxw, b, nxw
  unsigned long  tid   =  threadIdx.x;
  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  ////unsigned long  Nvol  =  gridDim.x;
  const int dir_limit = 3;
  if(index < Nvol){

  ///////const Geometry& geo = prop.geo();
  const Coordinate xl = geo.coordinate_from_index(index);

  //const T* wm1 = &psrc[src_off*12*12];

  for (int dir = 0; dir < 2*dir_limit; ++dir){
    //if(tid < 9)ls[dir*9+tid] = gf[dir*Nvol*9 + index*9 + tid];
    //if(tid < 9)ls[tid*6 + dir] = gf[dir*Nvol*9 + index*9 + tid];
    if(tid < 9)ls[(tid/3)*6*3 + dir*3 + tid%3] = gf[dir*Nvol*9 + index*9 + tid];
    //if(tid < 9)ls[(tid/3)][dir*3 + tid%3] = gf[dir*Nvol*9 + index*9 + tid];
  }

  for (int dir = -dir_limit; dir < dir_limit; ++dir){
    const Coordinate xl1 = coordinate_shifts(xl, dir);
    long src_off = geo_v.offset_from_coordinate(xl1);
    ////const WilsonMatrixT<T >& wm1 = prop1.get_elem(xl1);
    //const T* wm1 = &psrc[src_off*12*12];

    //int off = (dir+3)*9*16;
    //for(int r=0;r<9;r++)ps[off+r*16+tid] = wm1[r*16+tid];
    //for(int r=0;r<9;r++)ps[ (dir*3 + (r/3))*3*16 + (r%3)*16+tid] = wm1[r*16+tid];
    //for(int r=0;r<9;r++)ps[ tid*6*3*3 + (dir*3 + (r/3))*3 + (r%3)] = wm1[r*16+tid];
    //for(int r=0;r<9;r++)ps[ tid*6*3*3 + (dir*3 + (r/3))*3 + (r%3)] = psrc[src_off*12*12 + r*16+tid];
    //for(int r=0;r<9;r++)ps[ tid][((dir+3)*3 + (r/3))*3 + (r%3)] = psrc[src_off*12*12 + r*16+tid];

    //for(int r=0;r<9;r++)ps[ tid*6*3*3 + ((dir+3)*3 + (r/3))*3 + (r%3)] = psrc[src_off*12*12 + r*16+tid];
    for(int r=0;r<9;r++)ps[ tid*6*3*3 + (r%3)*18+ (dir+3)*3 + (r/3)] = psrc[src_off*12*12 + r*16+tid];

    //const Coordinate xl = geo.coordinate_from_index(tid);
    //WilsonMatrixT<T>& wm = prop.get_elem(xl);
  }
  //for(int r=0;r<9;r++)ws[ r*16 + tid] = 0;

  __syncthreads();

  //Eigen::Map< Eigen::Matrix< T, 3   , 6*3,Eigen::RowMajor> > lE(&ls[0]);
  //Eigen::Map< Eigen::Matrix< T, 6*3 , 3  ,Eigen::RowMajor> > pE(&ps[tid][0]);
  //Eigen::Map< Eigen::Matrix< T, 3   , 3  ,Eigen::RowMajor> > wE(&ws[tid][0]);

  //Eigen::Map< Eigen::Matrix< T, 3   , 6*3,Eigen::RowMajor> > lE(&ls[0]);
  //Eigen::Map< Eigen::Matrix< T, 6*3 , 3  ,Eigen::RowMajor> > pE(&ps[tid*6*3*3]);
  //Eigen::Map< Eigen::Matrix< T, 3   , 3  ,Eigen::RowMajor> > wE(&ws[tid*3*3]);
  //wE  = lE * pE;

  T* lp = &ls[0];
  T* pp = &ps[tid*6*3*3];
  ////T* wp = &ws[tid*3*3];
  long res_off = geo.offset_from_coordinate(xl);
  T* wm = &pres[res_off*12*12];

  ///for(int i=0;i<9;i++){wp[i] = 0;}

  T tmp = 0.0;
  for(int c0=0;c0<3;c0++)
  for(int c1=0;c1<3;c1++)
  {
    tmp = 0.0;
    for(int di=0;di<18;di++)
    {
      //tmp += lp[c0*18 + di] * pp[di*3 + c1];
      tmp += lp[c0*18 + di] * pp[c1*18 + di];
    }

    //wp[c0*3 + c1 ] = tmp;
    wm[ (c0*3+c1)*16 + tid] = norm*wm[ (c0*3+c1)*16 + tid] + norm*bw*tmp;
  }

  //long src_off = geo_v.offset_from_coordinate(xl);
  //T* pcpy = &psrc[src_off*12*12];
  //for(int r=0;r<9;r++){pcpy[ tid*9 + r] = wm[ tid*9 + r];}

  //for(int r=0;r<9;r++){wm[ r*16+tid] = norm*wm[ r*16+tid] + norm*bw*ws[tid][r];}
  //for(int r=0;r<9;r++){wm[ r*16+tid] = norm*wm[ r*16+tid] + norm*bw*ws[tid*9 + r];}

  ////for(int r=0;r<9;r++){wm[ r*16+tid] += bw*ws[tid*9 + r];wm[ r*16+tid] *= norm;}

  //const Coordinate xl = geo.coordinate_from_index(tid);
  //WilsonMatrixT<T>& wm = prop.get_elem(xl);
  }
}

template <class T>
__global__ void cpy_smear_global3(T* p0, T* p1, long Nvol, const Geometry& geo,const Geometry& geo_v)
{
  ////////======c, mxn, a, mxw, b, nxw
  unsigned long  tid   =  threadIdx.x;
  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  __shared__ T ps[9 * 16];
  ////unsigned long  Nvol  =  gridDim.x;
  if(index < Nvol){

  const Coordinate xl = geo.coordinate_from_index(index);
  long off0 = geo.offset_from_coordinate(xl);
  long off1 = geo_v.offset_from_coordinate(xl);
  T* psm  = &p0[off0*12*12];
  T* pcpy = &p1[off1*12*12];

  //for(int r=0;r<9;r++){ps[r*16 + tid] = psm[r*16 + tid];}
  for(int r=0;r<9;r++){ps[tid*9 + r] = psm[r*16 + tid];}
  __syncthreads();

  for(int r=0;r<9;r++){pcpy[r*16 + tid] = ps[ tid*9 + r];}

  }
}

void get_mapvq(const std::vector<CommPackInfo> &pack_infos, qlat::vector<long >& mapvq, int dir=0)
{
  std::vector<long > mapv;//mapv.resize(4*pack_infos.size());
  for (long i = 0; i < (long)pack_infos.size(); ++i) {
    const CommPackInfo& cpi = pack_infos[i];
    for(long j=0; j< cpi.size ; j++)
    {
      long bufi = cpi.buffer_idx + j;
      long fi   = cpi.offset     + j;
      if(dir == 0){mapv.push_back(bufi);mapv.push_back(  fi);}
      if(dir == 1){mapv.push_back(  fi);mapv.push_back(bufi);}
    }
  }
  //////qlat::vector<long > mapvq;
  mapvq.resize(mapv.size());
  for(size_t mi=0;mi< mapv.size();mi++){mapvq[mi] = mapv[mi];}
}

template <class M>
__global__ void shift_copy_global(M* pres, M* psrc, long* map, long long Nvol)
{
  //unsigned int   nt    = 12;
  //long  index =  blockIdx.x*blockDim.x + threadIdx.x;
  //__shared__ char buf[32];
  ////unsigned long  Nvol  =  gridDim.x;
  //int ng = 16;
  ///long Nv = Mend/ng;
  /////A type with size 16

  //qlat::Complex* r = (qlat::Complex*) &pres[map[index*2+0]];
  //qlat::Complex* s = (qlat::Complex*) &psrc[map[index*2+1]];
  //unsigned int Mend = sizeof(M)/sizeof(qlat::Complex);
  //if(index < Nvol)

  //r[tid] = s[tid];


  //pres[map[index*2+0]] = psrc[map[index*2+1]];
  //long off=0;
  //while(off < Mend){for(unsigned int i=0;i<ng;i++){r[tid*ng+i] = s[tid*ng + i];}off+=nt*ng;}

  //const Coordinate xl = geo.coordinate_from_index(index);
  //long off0 = geo.offset_from_coordinate(xl);
  //long off1 = geo_v.offset_from_coordinate(xl);
  //T* psm  = &p0[off0*12*12];
  //T* pcpy = &p1[off1*12*12];

  ////for(int r=0;r<9;r++){ps[r*16 + tid] = psm[r*16 + tid];}
  //for(int r=0;r<9;r++){ps[tid*9 + r] = psm[r*16 + tid];}
  //__syncthreads();

  //for(int r=0;r<9;r++){pcpy[r*16 + tid] = ps[ tid*9 + r];}

  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int  tid    =  threadIdx.x;
  unsigned int   nt    = blockDim.x;

  __shared__ long m[2];
  if(tid < 2){
  m[tid] = map[index*2+tid];
  m[tid] = map[index*2+tid];
  }
  __syncthreads();

  double* r = (double*) &pres[m[0]];
  double* s = (double*) &psrc[m[1]];
  //double* r = (double*) &pres[map[index*2+0]];
  //double* s = (double*) &psrc[map[index*2+1]];

  unsigned int Mend = sizeof(M)/sizeof(double);
  unsigned int off = tid;
  while(off < Mend){r[off] = s[off];off += nt;}

}



template <class M>
void shift_copy(M* res, M* src,qlat::vector<long > mapvq)
{

  TIMER("shift copy");
  long Mend = sizeof(M);qassert(Mend%sizeof(double) == 0);
  //long nt = Mend/16; qassert(nt < 256)
  long nt = 24;
  long long Nvol = mapvq.size()/2;
  /////print0("==Nvol %d \n", Nvol);
  ///long long sn   = (Nvol+nt-1)/nt;

  dim3 dimBlock(  nt,   1, 1);
  dim3 dimGrid( Nvol,  1, 1);
  //if(dir == 0){shift_copy_global<M ><<< dimGrid, dimBlock >>>(buf, &f.field[0], &mapvq[0],Nvol);}
  //if(dir == 1){shift_copy_global<M ><<< dimGrid, dimBlock >>>(&f.field[0], buf, &mapvq[0],Nvol);}
  //if(dir == 1){shift_copy_global<M ><<< dimGrid, dimBlock >>>(&f.field[0], buf, &mapvq[0],Nvol);}
  //shift_copy_global<M ><<< dimGrid, dimBlock >>>(buf, f, &mapvq[0], Nvol);
  shift_copy_global<M ><<< dimGrid, dimBlock >>>(res, src, &mapvq[0], Nvol);

  qacc_barrier(dummy);

  //#pragma omp parallel for
  //for (long i = 0; i < (long)pack_infos.size(); ++i) {
  //  const CommPackInfo& cpi = pack_infos[i];
  //  if(dir == 0){
  //  cudaMemcpyAsync(&buf[cpi.buffer_idx], &f.get_elem(cpi.offset),
  //          cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);}
  //  if(dir == 1){
  //  cudaMemcpyAsync(&f.get_elem(cpi.offset), &buf[cpi.buffer_idx], 
  //          cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);}

  //  //Memcpy(&send_buffer[cpi.buffer_idx], &f.get_elem(cpi.offset),
  //  //       cpi.size * sizeof(M));
  //}
  //qacc_barrier(dummy);

  /////int count = 0;


}

template <class M>
void refresh_expanded_GPU(Field<M>& f, const CommPlan& plan, M* send_buffer, M* recv_buffer, qlat::vector<long > &mapvq_send, qlat::vector<long > &mapvq_recv)
{
  const long total_bytes =
      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded");
  timer.flops += total_bytes / 2;
  //vector<M> send_buffer(plan.total_send_size);
  //vector<M> recv_buffer(plan.total_recv_size);
    
  //{
  //TIMER("shift copy");
  //#pragma omp parallel for
  //for (long i = 0; i < (long)plan.send_pack_infos.size(); ++i) {
  //  const CommPackInfo& cpi = plan.send_pack_infos[i];
  //  cudaMemcpyAsync(&send_buffer[cpi.buffer_idx], &f.get_elem(cpi.offset),
  //         cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);
  //  //Memcpy(&send_buffer[cpi.buffer_idx], &f.get_elem(cpi.offset),
  //  //       cpi.size * sizeof(M));
  //}
  //qacc_barrier(dummy);
  //}

  //#pragma omp parallel for
  //for (long i = 0; i < mapvq_send.size()/2; ++i) {
  //  cudaMemcpyAsync(&send_buffer[mapvq_send[i*2+0]], &f.get_elem(mapvq_send[i*2+1]),
  //         sizeof(M), cudaMemcpyDeviceToDevice);
  //}
  //qacc_barrier(dummy);


  ////shift_copy(send_buffer, f , plan.send_pack_infos , 0);
  shift_copy(send_buffer, &f.field[0] , mapvq_send);

  vector<MPI_Request> send_reqs(plan.send_msg_infos.size());
  vector<MPI_Request> recv_reqs(plan.recv_msg_infos.size());

  {
    //sync_node();
    TIMER_FLOPS("refresh_expanded-comm");
    timer.flops +=
        (plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;
    {
      TIMER("refresh_expanded-comm-init");
      const int mpi_tag = 10;
      for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = plan.recv_msg_infos[i];
        MPI_Irecv(&recv_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);
      }
      for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = plan.send_msg_infos[i];
        MPI_Isend(&send_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);
      }
    }
    MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
    //sync_node();
  }

  //{
  //TIMER_FLOPS("refresh_expanded-comm");
  //timer.flops +=
  //    (plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;

  //int local_MPI = qlat::get_num_node();
  //qlat::vector<int > send;
  //qlat::vector<int > spls;
  //qlat::vector<int > recv;
  //qlat::vector<int > rpls;

  //{
  //TIMER("refresh_expanded-comm-init");
  //send.resize(local_MPI);
  //recv.resize(local_MPI);
  //spls.resize(local_MPI);
  //rpls.resize(local_MPI);
  //for(int n=0;n<local_MPI;n++){
  //  send[n] = 0;
  //  recv[n] = 0;
  //  spls[n] = 0;
  //  rpls[n] = 0;
  //}

  //for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
  //  const CommMsgInfo& cmi = plan.recv_msg_infos[i];
  //  recv[cmi.id_node] = cmi.size * sizeof(M) ;
  //  rpls[cmi.id_node] = cmi.buffer_idx * sizeof(M);
  //}
  //for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
  //  const CommMsgInfo& cmi = plan.send_msg_infos[i];
  //  send[cmi.id_node] = cmi.size * sizeof(M) ;
  //  spls[cmi.id_node] = cmi.buffer_idx * sizeof(M);
  //}
  //}

  //{
  //  TIMER("refresh_expanded-comm-mpi");
  //  MPI_Alltoallv(&send_buffer[0],(int*) &send[0],(int*) &spls[0], MPI_BYTE,
  //                &recv_buffer[0],(int*) &recv[0],(int*) &rpls[0], MPI_BYTE, get_comm());
  //}
  //
  //
  //}

  //shift_copy(recv_buffer, f ,plan.recv_pack_infos , 1);
  //shift_copy(recv_buffer, f , mapvq_recv , 1);
  //shift_copy(recv_buffer, &f.field[0] , mapvq_recv);

  shift_copy(&f.field[0], recv_buffer , mapvq_recv);

  //#pragma omp parallel for
  //for (long i = 0; i < mapvq_recv.size()/2; ++i) {
  //  cudaMemcpyAsync(&f.get_elem(mapvq_recv[i*2+0]), &recv_buffer[mapvq_recv[i*2+1]],
  //         sizeof(M), cudaMemcpyDeviceToDevice);
  //}
  //qacc_barrier(dummy);


  //{
  //TIMER("shift copy");
  //#pragma omp parallel for
  //for (long i = 0; i < (long)plan.recv_pack_infos.size(); ++i) {
  //  const CommPackInfo& cpi = plan.recv_pack_infos[i];
  //  cudaMemcpyAsync(&f.get_elem(cpi.offset), &recv_buffer[cpi.buffer_idx],
  //         cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);
  //  //memcpy(&f.get_elem(cpi.offset), &recv_buffer[cpi.buffer_idx],
  //  //       cpi.size * sizeof(M));
  //}
  //qacc_barrier(dummy);
  //}

}    


template <class T>
void smear_propagator3(Propagator4dT<T>& prop, const qlat::vector<T >& gf,
                      const double width, const int step, Propagator4dT<T>& prop_buf)
{
  const double aw   = 3.0*width*width/(2*step);
  const double bw = width*width/(4.0*step - 6.0*width*width);

  const Geometry& geo = prop.geo();
  long Nvol = geo.local_volume();
  //const Geometry geo1 = geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));

  int nt = 16;
  long sn = long(std::sqrt(Nvol))+1;
  dim3 dimBlock(nt,   1, 1);
  dim3 dimGrid(  sn, sn, 1);

  const CommPlan& plan = get_comm_plan(set_marks_field_1, "", prop_buf.geo());
  //vector< WilsonMatrixT<T> > send_buffer(plan.total_send_size);
  //vector< WilsonMatrixT<T> > recv_buffer(plan.total_recv_size);
  WilsonMatrixT<T> * send_buffer;WilsonMatrixT<T> * recv_buffer;
  cudaMalloc((void **)&send_buffer, plan.total_send_size*sizeof(WilsonMatrixT<T> ));
  cudaMalloc((void **)&recv_buffer, plan.total_recv_size*sizeof(WilsonMatrixT<T> ));

  qlat::vector<long > mapvq_send;
  qlat::vector<long > mapvq_recv;
  get_mapvq(plan.send_pack_infos, mapvq_send, 0);
  get_mapvq(plan.recv_pack_infos, mapvq_recv, 1);

  //    printf("Compile time check:\n");
  //#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  //    printf("This MPI library has CUDA-aware support.\n", MPIX_CUDA_AWARE_SUPPORT);
  //#elif defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
  //    printf("This MPI library does not have CUDA-aware support.\n");
  //#else
  //    printf("This MPI library cannot determine if there is CUDA-aware support.\n");
  //#endif /* MPIX_CUDA_AWARE_SUPPORT */
  // 
  //    printf("Run time check:\n");
  //#if defined(MPIX_CUDA_AWARE_SUPPORT)
  //    if (1 == MPIX_Query_cuda_support()) {
  //        printf("This MPI library has CUDA-aware support.\n");
  //    } else {
  //        printf("This MPI library does not have CUDA-aware support.\n");
  //    }   
  //#else /* !defined(MPIX_CUDA_AWARE_SUPPORT) */
  //    printf("This MPI library cannot determine if there is CUDA-aware support.\n");
  //#endif /* MPIX_CUDA_AWARE_SUPPORT */

  {
  TIMER_FLOPS("==smear in");
  long long Tfloat = 0;
  double mem       = 0.0;
  {long long Lat = prop.geo().local_volume();
  int nsrc = 12;
  long long vGb = Lat *nsrc*4;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
  timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

  for (int i = 0; i < step; ++i) {
    //print0("step %d \n",i);
    TIMER("Matrix multiply");
    //fflush_MPI();
    T* pres = &prop.get_elem(0)(0,0);
    T* psrc = &prop_buf.get_elem(0)(0,0);

    {
    TIMER("Copy prop");
    //qacc_for(index,  geo.local_volume(),{
    //  const Coordinate xl = geo.coordinate_from_index(index);
    //  WilsonMatrixT<T>& v1 = prop_buf.get_elem(xl);
    //  WilsonMatrixT<T>& v0 = prop.get_elem(xl);
    //  v1 = v0;
    //});
    //#pragma omp parallel for
    //for(long index=0;index<Nvol;index++)
    //{
    //  const Coordinate xl = geo.coordinate_from_index(index);
    //  WilsonMatrixT<T>& v1 = prop_buf.get_elem(xl);
    //  WilsonMatrixT<T>& v0 = prop.get_elem(xl);
    //  cudaMemcpyAsync(&v1(0,0), &v0(0,0), 12*12*sizeof(T),cudaMemcpyDeviceToDevice);
    //}
    //qacc_barrier(dummy);

    //cudaMemcpy(psrc, pres, Nvol*12*12*sizeof(T),cudaMemcpyDeviceToDevice);

    cpy_smear_global3<T ><<< dimGrid, dimBlock >>>(pres, psrc, Nvol, geo,  prop_buf.geo());
    qacc_barrier(dummy);

    }

    ////refresh_expanded_1(prop1);
    //refresh_expanded_1_GPU(prop_buf);
    refresh_expanded_GPU(prop_buf, plan, send_buffer, recv_buffer, mapvq_send, mapvq_recv);

    //refresh_expanded_GPU1(prop_buf, plan, &send_buffer[0], &recv_buffer[0]);

    {
    TIMER("kernal");
    smear_global3<T ><<< dimGrid, dimBlock >>>(pres, psrc, &gf[0], bw, 1-aw, Nvol, geo,  prop_buf.geo());
    qacc_barrier(dummy);
    }
  }
  }

  //qacc_for(index,  geo.local_volume(),{
  //  WilsonMatrixT<T>& wm = prop.get_elem(index);
  //  wm *= std::pow(1-aw, step);
  //});

  cudaFree(send_buffer);
  cudaFree(recv_buffer);


}


//__global__ void smear_global2(Propagator4dT<T>& prop, const Propagator4dT<T> &prop1, const T* gf, const double bw, long Nvol,const Geometry& geo,const Geometry& geo_v)
//__global__ void smear_global2(T* pres, Propagator4dT<T> &prop1, const T* gf, const double bw, long Nvol,const Geometry& geo)
template <class T>
__global__ void smear_global2(T* pres, const T* psrc, const T* gf, const double bw, long Nvol,const Geometry& geo,const Geometry& geo_v)
{

  long  index =  blockIdx.x*blockDim.x + threadIdx.x;
  const int dir_limit = 3;

  //const Geometry geo = prop.geo();
  if(index < Nvol){
  const Coordinate xl = geo.coordinate_from_index(index);
  long res_off = geo.offset_from_coordinate(xl);

  //WilsonMatrixT<T>& wm = prop.get_elem(xl);
  //T* wm = &pres[res_off*12*12];
  T* wm = &pres[res_off*12*12];

  //for(int c0=0;c0<3;c0++)
  //for(int c1=0;c1<3;c1++)
  //for(int c2=0;c2<3;c2++)
  //for(int j=0;j<16;j++)
  //{
  //  //wm[(c0*3+c2)*16+j] += (bw * lp[c0*3+c1] * wm1[(c1*3+c2)*16 + j ]);
  //  wm[(c0*3+c2)*16+j] += 0.0;
  //}

  //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  &wm(0,0));
  //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  wm.p);
  Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  wm);

  //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor> tmp_w;
  //for(int i=0;i<3;i++)for(int j=0;j<3*16;j++){tmp_w(i,j) = 0.0;}
  //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  wm);

  for (int dir = -dir_limit; dir < dir_limit; ++dir) {
    const Coordinate xl1 = coordinate_shifts(xl, dir);
    long src_off = geo_v.offset_from_coordinate(xl1);
    //const WilsonMatrixT<T>& wm1 = prop1.get_elem(xl1);
    const T* wm1 = &psrc[src_off*12*12];
    const T* lp = &gf[(dir+3)*Nvol*9 + index*9 + 0];
    //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wm1.p);
    //Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) &wm1(0,0));

    Eigen::Matrix<T, 3, 3, Eigen::RowMajor>&  lE = *((Eigen::Matrix<T, 3, 3, Eigen::RowMajor>*) lp);
    Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wm1);

    wmE += bw * (lE * pE);
    //tmp_w += bw * (lE * pE);
    //for(int c0=0;c0<3;c0++)
    //for(int c1=0;c1<3;c1++)
    //for(int c2=0;c2<3;c2++)
    //for(int j=0;j<16;j++)
    //{
    //  //wm[(c0*3+c2)*16+j] += (bw * lp[c0*3+c1] * wm1[(c1*3+c2)*16 + j ]);
    //  wm[(c0*3+c2)*16+j] += 0.0;
    //}
  }

  //for(int i=0;i<3;i++)for(int j=0;j<3*16;j++){wm[i*3*16+j] += tmp_w(i,j);}

  }

}

template <class T>
void smear_propagator2(Propagator4dT<T>& prop, const qlat::vector<T >& gf,
                      const double width, const int step)
{
  const double aw   = 3.0*width*width/(2*step);
  const double bw = width*width/(4.0*step - 6.0*width*width);

  const Geometry& geo = prop.geo();
  long Nvol = geo.local_volume();
  const Geometry geo1 = geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  rotate_prop(prop);

  Propagator4dT<T> prop1;
  prop1.init(geo1);

  //long sn = long(std::sqrt(Nvol))+1;
  //dim3 dimGrid(  sn, sn, 1);
  //dim3 dimBlock(16,   1, 1);

  int nt = 32;
  long nb = (Nvol+nt-1)/nt;
  dim3 dimBlock( nt,  1, 1);
  dim3 dimGrid(  nb,  1, 1);


  for (int i = 0; i < step; ++i) {
    print0("step %d \n",i);
    TIMER("Matrix multiply");
    fflush_MPI();
    //prop1 = prop;
    qacc_for(index,  geo.local_volume(),{
      const Coordinate xl = geo.coordinate_from_index(index);
      //WilsonMatrixT<T>& wm = prop.get_elem(xl);
      //const Coordinate xl = geo.coordinate_from_index(index);
      //Vector<M> v = f.get_elems(xl);
      WilsonMatrixT<T>& v1 = prop1.get_elem(xl);
      WilsonMatrixT<T>& v0 = prop.get_elem(xl);
      v1 = v0;
    });

    refresh_expanded_1(prop1);
    //const int dir_limit = 3;
    //qacc_for(index,  geo.local_volume(),{
    //  const Coordinate xl = geo.coordinate_from_index(index);
    //  WilsonMatrixT<T>& wm = prop.get_elem(xl);
    //  Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  wm.p);

    //  for (int dir = -dir_limit; dir < dir_limit; ++dir) {
    //  //WilsonMatrixT<T>& wm = prop.get_elem(index);
    //  //WilsonMatrixT<T>& wm1 = prop1.get_elem(index);

    //  const Coordinate xl1 = coordinate_shifts(xl, dir);
    //  WilsonMatrixT<T>& wm1 = prop1.get_elem(xl1);

    //  const T* lp = &gf[(dir+3)*Nvol*9 + index*9 + 0]; 
    //  Eigen::Matrix<T, 3, 3, Eigen::RowMajor>&  lE = *((Eigen::Matrix<T, 3, 3, Eigen::RowMajor>*) lp);
    //  Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wm1.p);
    //  wmE += bw * (lE * pE);
    //  }   
    //}); 

    //smear_global2<T ><<< dimGrid, dimBlock >>>(prop, prop1, &gf[0], bw, Nvol, geo);
    T* pres = &prop.get_elem(0)(0,0);
    T* psrc = &prop1.get_elem(0)(0,0);
    smear_global2<T ><<< dimGrid, dimBlock >>>(pres, psrc, &gf[0], bw, Nvol, geo, prop1.geo());
    //for(long index=0;index<Nvol;index++){
    //  for(int j=0;j<144;j++)pres[index*12*12+j] = 0;
    //}
    qacc_barrier(dummy);
  }

  qacc_for(index,  geo.local_volume(),{
    ////const Coordinate xl = geo.coordinate_from_index(index);
    ////WilsonMatrixT<T>& wm = prop.get_elem(xl);
    WilsonMatrixT<T>& wm = prop.get_elem(index);
    wm *= std::pow(1-aw, step);
  });

  rotate_prop(prop,1);
}


template <class T >
void smear_propagator1(Propagator4dT<T>& prop, const qlat::vector<T >& gf,
                      const double width, const int step)
{

  const double aw   = 3.0*width*width/(2*step);
  const double bw = width*width/(4.0*step - 6.0*width*width);

  const Geometry& geo = prop.geo();
  long Nvol = geo.local_volume();
  const Geometry geo1 = geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int dir_limit = 3;
  rotate_prop(prop);

  Propagator4dT<T> prop1;
  prop1.init(geo1);
  //prop1 = prop;
  //refresh_expanded_1(prop1);
  for (int i = 0; i < step; ++i) {
    print0("step %d \n",i);
    TIMER("Matrix multiply");
    fflush_MPI();
    prop1 = prop;
    refresh_expanded_1(prop1);

    qacc_for(index,  geo.local_volume(),{
      const Coordinate xl = geo.coordinate_from_index(index);
      WilsonMatrixT<T>& wm = prop.get_elem(xl);
      Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>& wmE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*)  wm.p);

      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      //WilsonMatrixT<T>& wm = prop.get_elem(index);
      //WilsonMatrixT<T>& wm1 = prop1.get_elem(index);

      const Coordinate xl1 = coordinate_shifts(xl, dir);
      WilsonMatrixT<T>& wm1 = prop1.get_elem(xl1);

      const T* lp = &gf[(dir+3)*Nvol*9 + index*9 + 0];
      Eigen::Matrix<T, 3, 3, Eigen::RowMajor>&  lE = *((Eigen::Matrix<T, 3, 3, Eigen::RowMajor>*) lp);
      Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>&  pE = *((Eigen::Matrix<T, 3, 3*16, Eigen::RowMajor>*) wm1.p);
      wmE += bw * (lE * pE);
      }
    });
  }

  qacc_for(index,  geo.local_volume(),{
    const Coordinate xl = geo.coordinate_from_index(index);
    WilsonMatrixT<T>& wm = prop.get_elem(xl);
    wm *= std::pow(1-aw, step);
  });

  rotate_prop(prop,1);

  //prop[x+-1,y+-1,z+-1]
  //pxp1 + lxp1
  //pxm1 + lxm1
  //pyp1 + lxp1
  //pym1 + lym1
  //pzp1 + lzp1
  //pzm1 = lzm1

}

template <class T>
void smear_propagator_gpu(Propagator4dT<T>& prop, const qlat::vector<T >& gf, const double width, const int step, int mode=1)
{
  // gf1 is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size

  TIMER_FLOPS("==smear propagator");
  long long Tfloat = 0;
  double mem       = 0.0;

  {long long Lat = prop.geo().local_volume();
  int nsrc = 12;
  long long vGb = Lat *nsrc*4;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
  timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));


  const Geometry geo1 = geo_resize(prop.geo(), Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  Propagator4dT<T> prop_buf;
  prop_buf.init(geo1);

  if(mode == 1)smear_propagator1(prop, gf, width, step);
  if(mode == 2)smear_propagator2(prop, gf, width, step);

  rotate_prop(prop);
  if(mode == 3){
    TIMER_FLOPS("==compute time");
    smear_propagator3(prop, gf, width, step, prop_buf);
    timer.flops += Tfloat;
  }
  rotate_prop(prop,1);

}



template <class T>
void smear_propagator_gwu_convension_cpu(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                      const double width, const int step,int mode=0)
{
  // gf1 is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf1, gf)
  // prop is of qnormal size
  ////For complex numbers addition and subtraction require two flops, and multiplication and division require six flops
  ///complex multi 6 + plus 2

  TIMER_FLOPS("==smear propagator");
  long long Tfloat = 0;
  double mem       = 0.0;

  {long long Lat = prop.geo().local_volume();
  int nsrc = 12;
  long long vGb = Lat *nsrc*4;
  int Fcount = 3*(3*6 + 2*2); 
  int direction   = 6;
  Tfloat = step*direction*vGb*Fcount;
  mem = (Lat*nsrc*12 + Lat*4*9)*8.0;}
  timer.flops += Tfloat;
  print0("Memory size %.3e GB, %.3e Gflop \n", 
    mem/(1024.0*1024*1024), Tfloat/(1024.0*1024*1024));

  if (0 == step) {
    return;
  }
  const double aw   = 3.0*width*width/(2*step);
  const double coef = aw;
  if(mode == 0)smear_propagator(prop, gf1, coef, step);
  ////if(mode == 1)smear_propagator1(prop, gf1, coef, step);

}


template <class T>
void prop_to_EigenV(Propagator4dT<T>& prop, EigenV &pE, int dir=0)
{
  TIMER("Copy prop");
  const Geometry& geo = prop.geo();
  long long Nvol =  geo.local_volume();
  if(dir == 0){pE.resize(12*12* Nvol);}

  qacc_for(isp, long(Nvol),{
    qlat::WilsonMatrixT<T>& v0 =  prop.get_elem(isp);

    for(int c0 = 0;c0 < 3; c0++)
    for(int d0 = 0;d0 < 4; d0++)
    for(int c1 = 0;c1 < 3; c1++)
    for(int d1 = 0;d1 < 4; d1++)
    {
      LInt off = (d1*3+c1)*12+d0*3+c0;
      if(dir==0)pE[off*Nvol + isp] = v0(d0*3 + c0, d1*3 + c1);
      if(dir==1)v0(d0*3 + c0, d1*3 + c1) = pE[off*Nvol + isp];
    }
  });
}


template <class T>
void links_to_EigenV(EigenV &link, GaugeFieldT<T>& gf, int dir = 0)
{
  const Geometry& geo = gf.geo();
  long long Nvol =  geo.local_volume();
  if(dir == 0){link.resize(4*9* Nvol);}

  qacc_for(index,  geo.local_volume(),{
    for(int dir=0;dir<4;dir++){
      ColorMatrixT<T>  cm = gf.get_elem(index, dir);
      for (int c1 = 0; c1 < NUM_COLOR; ++c1) 
      for (int c2 = 0; c2 < NUM_COLOR; ++c2) 
      {
         long long off = (dir*9 + c1*3 + c2)*Nvol + index;
         if(dir == 0)link[ off] = cm(c1, c2);
         if(dir == 1)cm(c1, c2) = link[ off];
      }
    }
  });

  //  for (int s1 = 0; s1 < 4; ++s1) {
  //    for (int s2 = 0; s2 < 4; ++s2) {
  //      for (int c1 = 0; c1 < NUM_COLOR; ++c1) {
  //        for (int c2 = 0; c2 < NUM_COLOR; ++c2) {
  //          for (int c3 = 0; c3 < NUM_COLOR; ++c3) {
  //            ret(s1 * NUM_COLOR + c1, s2 * NUM_COLOR + c2) +=
  //                cm(c1, c3) * m(s1 * NUM_COLOR + c3, s2 * NUM_COLOR + c2);
  //          }
  //        }
  //      }
  //    }
  //  }

}


//template <class M>
//void refresh_expanded_GPU1(Field<M>& f, const CommPlan& plan, M* send_buffer, M* recv_buffer)
//{
//  const long total_bytes =
//      (plan.total_recv_size + plan.total_send_size) * sizeof(M);
//  if (0 == total_bytes) {
//    return;
//  }
//  TIMER_FLOPS("refresh_expanded");
//  timer.flops += total_bytes / 2;
//
//  {
//  TIMER("shift copy");
//  #pragma omp parallel for
//  for (long i = 0; i < (long)plan.send_pack_infos.size(); ++i) {
//    const CommPackInfo& cpi = plan.send_pack_infos[i];
//    cudaMemcpyAsync(&send_buffer[cpi.buffer_idx], &f.get_elem(cpi.offset),
//           cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);
//    //Memcpy(&send_buffer[cpi.buffer_idx], &f.get_elem(cpi.offset),
//    //       cpi.size * sizeof(M));
//  }
//  }
//  qacc_barrier(dummy);
//
//  //{
//  //  sync_node();
//  //  TIMER_FLOPS("refresh_expanded-comm");
//  //  timer.flops +=
//  //      (plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;
//  //  vector<MPI_Request> send_reqs(plan.send_msg_infos.size());
//  //  vector<MPI_Request> recv_reqs(plan.recv_msg_infos.size());
//  //  {
//  //    TIMER("refresh_expanded-comm-init");
//  //    const int mpi_tag = 10;
//  //    for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
//  //      const CommMsgInfo& cmi = plan.recv_msg_infos[i];
//  //      MPI_Irecv(&recv_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
//  //                cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);
//  //    }
//  //    for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
//  //      const CommMsgInfo& cmi = plan.send_msg_infos[i];
//  //      MPI_Isend(&send_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
//  //                cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);
//  //    }
//  //  }
//  //  MPI_Waitall(recv_reqs.size(), recv_reqs.data(), MPI_STATUS_IGNORE);
//  //  MPI_Waitall(send_reqs.size(), send_reqs.data(), MPI_STATUS_IGNORE);
//  //  sync_node();
//  //}
//
//  int local_MPI = qlat::get_num_node();
//  qlat::vector<int > send;
//  qlat::vector<int > spls;
//  qlat::vector<int > recv;
//  qlat::vector<int > rpls;
//
//  send.resize(local_MPI);
//  recv.resize(local_MPI);
//  spls.resize(local_MPI);
//  rpls.resize(local_MPI);
//  for(int n=0;n<local_MPI;n++){
//    send[n] = 0;
//    recv[n] = 0;
//    spls[n] = 0;
//    rpls[n] = 0;
//  }
//
//  double v0 = 0.0;
//  double v1 = 0.0;
//  for (size_t i = 0; i < plan.recv_msg_infos.size(); ++i) {
//    const CommMsgInfo& cmi = plan.recv_msg_infos[i];
//    recv[cmi.id_node] = cmi.size * sizeof(M) ;
//    rpls[cmi.id_node] = cmi.buffer_idx * sizeof(M);
//    //rpls[cmi.id_node] = cmi.id_node*cmi.size * sizeof(M);
//    //rpls[cmi.id_node] = cmi.size * sizeof(M);
//    v0 += cmi.size * sizeof(M) * 1.0;
//    //MPI_Irecv(&recv_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
//    //          cmi.id_node, mpi_tag, get_comm(), &recv_reqs[i]);
//  }
//  for (size_t i = 0; i < plan.send_msg_infos.size(); ++i) {
//    const CommMsgInfo& cmi = plan.send_msg_infos[i];
//    send[cmi.id_node] = cmi.size * sizeof(M) ;
//    spls[cmi.id_node] = cmi.buffer_idx * sizeof(M);
//    //spls[cmi.id_node] = cmi.size * sizeof(M);
//    v1 += cmi.size * sizeof(M) * 1.0;
//    //MPI_Isend(&send_buffer[cmi.buffer_idx], cmi.size * sizeof(M), MPI_BYTE,
//    //          cmi.id_node, mpi_tag, get_comm(), &send_reqs[i]);
//  }
//  ////v0 = v0/(1024.0*1024.0*1024.0);
//  ////v1 = v1/(1024.0*1024.0*1024.0);
//  //print0("rank %d, Size memory recv %.3e GB, send %.3e GB \n", 
//  //  qlat::get_id_node(), v0/(1024.0*1024.0*1024.0), v1/(1024.0*1024.0*1024.0));
// 
//  {
//    TIMER_FLOPS("refresh_expanded-comm");
//    ////(plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;
//    ///timer.flops += (v0);
//    timer.flops += (plan.total_recv_size + plan.total_send_size) * sizeof(M) / 2;
//
//    //for(int ranklocal=0;ranklocal<local_MPI;ranklocal++)
//    //{
//    //  send[ranklocal] = biva*N0*N1*N2*civ*2;
//    //  spls[ranklocal] = biva*N0*N1*N2*civ*2*ranklocal;
//    //  recv[ranklocal] = biva*N0*N1*N2*civ*2;
//    //  rpls[ranklocal] = biva*N0*N1*N2*civ*2*ranklocal;
//    //}
//
//    //Ftype* src = (Ftype*) &send_buffer[0];
//    //Ftype* res = (Ftype*) &recv_buffer[0];
//    //if(src == NULL){src = (Ftype*) &sendbuf[0];}
//    //if(res == NULL){res = (Ftype*) &recvbuf[0];}
//    // cmi.size * sizeof(M)
//    // cmi.size * sizeof(M)
//    //#pragma acc kernels deviceptr(send_buffer, recv_buffer)
//    //#pragma acc host_data use_device (send_buffer, recv_buffer, send, spls, recv, rpls)
//    MPI_Alltoallv(&send_buffer[0],(int*) &send[0],(int*) &spls[0], MPI_BYTE,
//                  &recv_buffer[0],(int*) &recv[0],(int*) &rpls[0], MPI_BYTE, get_comm());
//  
//  }
//
//  {
//  TIMER("shift copy");
//  #pragma omp parallel for
//  for (long i = 0; i < (long)plan.recv_pack_infos.size(); ++i) {
//    const CommPackInfo& cpi = plan.recv_pack_infos[i];
//    cudaMemcpyAsync(&f.get_elem(cpi.offset), &recv_buffer[cpi.buffer_idx],
//           cpi.size * sizeof(M), cudaMemcpyDeviceToDevice);
//    //memcpy(&f.get_elem(cpi.offset), &recv_buffer[cpi.buffer_idx],
//    //       cpi.size * sizeof(M));
//  }
//  qacc_barrier(dummy);
//  }
//}
//template <class M>
//void refresh_expanded_1_GPU(Field<M>& f)
//{
//  const CommPlan& plan = get_comm_plan(set_marks_field_1, "", f.geo());
//  refresh_expanded_GPU(f, plan);
//}

struct smear_fun{

  const Geometry* geop;
  const Geometry* geop_v;
  CommPlan plan;

  void* send_bufferV;
  void* recv_bufferV;
  qlat::vector<long > mapvq_send;
  qlat::vector<long > mapvq_recv;
  qlat::vector<MPI_Request> send_reqs;
  qlat::vector<MPI_Request> recv_reqs;
  unsigned int bfac;
  unsigned int Tsize;
  unsigned long Nvol;

  qlat::vector<long > map_buf;
  qlat::vector<long > map_bufD;
  qlat::vector<int> Nv,nv,mv;

  int fft_copy;


  smear_fun(){
    send_bufferV = NULL;
    recv_bufferV = NULL;
    bfac  = 0;
    Tsize = 0;
    Nvol  = 0;
    fft_copy = 0;
  }

  void init_mem(const Geometry& geo_a, const Geometry& geov_a, const int multiplicity_a, const int Tsize_a)
  {
    geop   = &geo_a;
    geop_v = &geov_a;
    plan = get_comm_plan(set_marks_field_1, "", geov_a);
    bfac = multiplicity_a;
    Tsize = Tsize_a;

    Nv.resize(4);nv.resize(4);mv.resize(4);
    for(int i=0;i<4;i++){Nv[i]=geop->node_site[i];nv[i] = geop->node_site[i] * geop->geon.size_node[i];}
    for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}

    cudaMalloc((void **)&send_bufferV, plan.total_send_size*bfac*Tsize);
    cudaMalloc((void **)&recv_bufferV, plan.total_recv_size*bfac*Tsize);

    get_mapvq(plan.send_pack_infos, mapvq_send, 0);
    get_mapvq(plan.recv_pack_infos, mapvq_recv, 1);
    send_reqs.resize(plan.send_msg_infos.size());
    recv_reqs.resize(plan.recv_msg_infos.size());

    Nvol = geop->local_volume();

    map_buf.resize(Nvol);
    for(long index=0;index<Nvol;index++){
      const Coordinate xl = geop->coordinate_from_index(index);
      map_buf[index] = geop_v->offset_from_coordinate(xl);
    }

    map_bufD.resize(8*Nvol);
    for(int dir=-4;dir<4; dir++)
    for(long index=0;index<Nvol;index++)
    {
      const Coordinate xl = geop->coordinate_from_index(index);
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      map_bufD[(dir+4)*Nvol + index] = geop_v->offset_from_coordinate(xl1);
    }
    /////WilsonMatrixT<T> * send_buffer;WilsonMatrixT<T> * recv_buffer;
  }

  void init_distribute(const Geometry& geo_a)
  {
    fft_copy = 1;
    ///////T distributed, spatial z,y,x
    geop   = &geo_a;
    Nv.resize(4);nv.resize(4);mv.resize(4);
    for(int i=0;i<4;i++){Nv[i]=geop->node_site[i];nv[i] = geop->node_site[i] * geop->geon.size_node[i];}
    for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}


    Nvol = Nv[3] * nv[2]*nv[1]*nv[0];
    std::vector<int > Nn;Nn.resize(4);
    for(int i=0;i<3;i++){Nn[i] = nv[i];}Nn[3] = Nv[3];

    map_buf.resize(Nvol);
    for(long index=0;index<Nvol;index++){
      map_buf[index] = index;
    }

    std::vector<int > xl;xl.resize(4);

    map_bufD.resize(8*Nvol);
    for(int dir=-4;dir<4; dir++)
    for(long index=0;index<Nvol;index++)
    {

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

      map_bufD[(dir+4)*Nvol + index] = ((xl[3]*Nn[2]+xl[2])*Nn[1] + xl[1])*Nn[0] + xl[0];
    }
    /////WilsonMatrixT<T> * send_buffer;WilsonMatrixT<T> * recv_buffer;
  }


  ~smear_fun(){
    /////if(plan != NULL){delete plan;}
    if(send_bufferV != NULL)cudaFree(send_bufferV);
    if(recv_bufferV != NULL)cudaFree(recv_bufferV);
  }

};

template <class T>
__global__ void shift_copy_global(T* pres, T* psrc, long* map, long long Nvol, unsigned int bfac)
{
  unsigned long  index =  blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int  tid    =  threadIdx.x;
  unsigned int   nt    = blockDim.x;

  __shared__ long m[2];
  if(tid < 2){
  m[tid] = map[index*2+tid];
  m[tid] = map[index*2+tid];
  }
  __syncthreads();

  double* r = (double*) &pres[m[0]*bfac];
  double* s = (double*) &psrc[m[1]*bfac];
  unsigned int Mend = (bfac*sizeof(T))/sizeof(double);
  unsigned int off = tid;
  while(off < Mend){r[off] = s[off];off += nt;}
}


template <class T>
void shift_copy(T* res, T* src,qlat::vector<long >& mapvq, unsigned int bfac)
{

  TIMER("shift copy");
  long Mend = (bfac*sizeof(T));qassert(Mend%sizeof(double) == 0);
  long nt = 24;
  long long Nvol = mapvq.size()/2;

  dim3 dimBlock(  nt,   1, 1);
  dim3 dimGrid( Nvol,  1, 1);
  shift_copy_global<T ><<< dimGrid, dimBlock >>>(res, src, &mapvq[0], Nvol, bfac);
  qacc_barrier(dummy);
}


template <class T>
void refresh_expanded_GPU(T* f, smear_fun& smf)
//(const CommPlan& plan, M* send_buffer, M* recv_buffer, qlat::vector<long > &mapvq_send, qlat::vector<long > &mapvq_recv)
{
  if(smf.fft_copy == 1){return ;}
  const long total_bytes =
      (smf.plan.total_recv_size + smf.plan.total_send_size) * smf.bfac * sizeof(T);
  if (0 == total_bytes) {
    return;
  }
  TIMER_FLOPS("refresh_expanded");
  timer.flops += total_bytes / 2;
  T* send_buffer = (T*) smf.send_bufferV;
  T* recv_buffer = (T*) smf.recv_bufferV;
  shift_copy(send_buffer, f , smf.mapvq_send, smf.bfac);

  ///vector<MPI_Request> send_reqs(plan.send_msg_infos.size());
  ///vector<MPI_Request> recv_reqs(plan.recv_msg_infos.size());

  {
    //sync_node();
    TIMER_FLOPS("refresh_expanded-comm");
    timer.flops +=
        (smf.plan.total_recv_size + smf.plan.total_send_size) * smf.bfac * sizeof(T) / 2;
    {
      TIMER("refresh_expanded-comm-init");
      const int mpi_tag = 10;
      for (size_t i = 0; i < smf.plan.recv_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = smf.plan.recv_msg_infos[i];
        MPI_Irecv(&recv_buffer[cmi.buffer_idx*smf.bfac], cmi.size * smf.bfac * sizeof(T), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), &smf.recv_reqs[i]);
      }
      for (size_t i = 0; i < smf.plan.send_msg_infos.size(); ++i) {
        const CommMsgInfo& cmi = smf.plan.send_msg_infos[i];
        MPI_Isend(&send_buffer[cmi.buffer_idx*smf.bfac], cmi.size * smf.bfac * sizeof(T), MPI_BYTE,
                  cmi.id_node, mpi_tag, get_comm(), &smf.send_reqs[i]);
      }
    }
    MPI_Waitall(smf.recv_reqs.size(), smf.recv_reqs.data(), MPI_STATUS_IGNORE);
    MPI_Waitall(smf.send_reqs.size(), smf.send_reqs.data(), MPI_STATUS_IGNORE);
    //sync_node();
  }
  shift_copy(f, recv_buffer , smf.mapvq_recv, smf.bfac);
}







}
#endif
