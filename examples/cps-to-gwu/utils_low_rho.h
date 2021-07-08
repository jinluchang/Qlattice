#ifndef utils_low_rho_h
#define utils_low_rho_h
#pragma once


#include "io_gwu.h"
#include <qlat/reduce_V.h>
///#include "cach_reduce.h"
//#include "reduce_V_dev.h"
#include "gammas.h"

#define EigenVq Eigen::Matrix< Complexq, 1, Eigen::Dynamic ,Eigen::RowMajor>
#define EigenAq Eigen::Array< Complexq, Eigen::Dynamic , 1>
#define Aoper  16

__device__ __constant__  signed char  Gmap0C[32];
__device__ __constant__  signed char  Gmap1C[32];
__device__ __constant__  signed char  Gmap2C[32];
__device__ __constant__  signed char  Gmap3C[32];
///__device__ __constant__  Complexq  NabC[16*128];


__global__ void multiplyNab_global(const Complexq* Nab, Ftype *Mres,const Ftype *Mvalues,const int nt,const int nmass,const unsigned long bufN0){

  __shared__ signed char G0[32];
  __shared__ signed char G1[32];
  __shared__ signed char G2[32];
  __shared__ signed char G3[32];
  extern __shared__ Complexq NMv_values[];
  ///__shared__ Complexq Nabv_multi[16*3];

  unsigned int tid = threadIdx.x;
  long bi = blockIdx.x;
  long bN = gridDim.x;
  ///unsigned int bi = ji/Aoper;
  ///unsigned int ipr  = ji%Aoper;

  ////Load gammas
  int off = tid;
  while(off < 32){G0[off] = Gmap0C[off];off += blockDim.x;}
  off = tid;
  while(off < 32){G1[off] = Gmap1C[off];off += blockDim.x;}
  off = tid;
  while(off < 32){G2[off] = Gmap2C[off];off += blockDim.x;}
  off = tid;
  while(off < 32){G3[off] = Gmap3C[off];off += blockDim.x;}

  //off = tid;unsigned long offAB = bi*nt*16;
  //while(off < nt*16){NMv_values[off] = Nab[offAB + off];off += blockDim.x;}
  ////off = tid;
  unsigned long offAB = bufN0*16;off = bi*16;
  if(tid < 16){for(int ti=0;ti<nt;ti++)NMv_values[ti*16 + tid] = Nab[ti*offAB + off + tid];}
  ///off = tid;int offT = nt*16;
  ///while(off<nmass*2){NMv_values[offT + off] = Mvalues[off];off += blockDim.x;}
  __syncthreads();

  const Complexq *Av,*Bv;
  Complexq v0[16];
  Complexq v1[16];
  for(int ipr=0;ipr<16;ipr++){v0[ipr] = 0.0;}
  for(int ipr=0;ipr<16;ipr++){v1[ipr] = 0.0;}

  unsigned int toff = tid;
  /////unsigned int t0 = tid;
  for(unsigned int t0=0;t0<nt;t0++)
  {
    unsigned int t1 = (t0+toff)%nt;
    ///Bv = &Avs[0];
    Bv = &NMv_values[t0*16 + 0];
    Av = &NMv_values[t1*16 + 0];

    for(int ipr=0;ipr<16;ipr++){
      v0[ipr] += Av[G0[ipr*2+0]]*qlat::qconj(Bv[G1[ipr*2+0]])
       *Ftype(G0[ipr*2+1])*Ftype(         G1[ipr*2+1]);
      v1[ipr] += Av[G2[ipr*2+0]]*qlat::qconj(Bv[G3[ipr*2+0]])
       *Ftype(G2[ipr*2+1])*Ftype(         G3[ipr*2+1]);
    }
  }
  __syncthreads();

  off = tid;offAB = bi*nmass*2;
  Ftype *buf = (Ftype*) &NMv_values[0];
  Ftype *src = (Ftype*) &Mvalues[offAB];
  while(off<nmass*2){buf[off] = src[off];off += blockDim.x;}
  __syncthreads();

  unsigned long shiftM = bufN0*nt;
  offAB = ((0*nmass+0)*bufN0 + bi)*nt + toff;
  for(int ipr=0;ipr<16;ipr++)
  {
    Ftype v00 = v0[ipr].real();
    Ftype v10 = v1[ipr].real();
    for(int mi=0;mi<nmass;mi++)
    {
      ////long jobN = Aoper*nt;
      Mres[offAB] += (buf[mi*2+0]*v00 + buf[mi*2+1]*v10);
      offAB += shiftM;
      ///mi += blockDim.x;
    }
  }

}

__global__ void prodab_global(const Complexq *a,const Complexq *b, Complexq *fd,const int Nvol,const int Nsum)
{
  Complexq as[12];
  Complexq bs[12];
  extern __shared__ Complexq resab[];
  unsigned int tid = threadIdx.x;
  unsigned long isp = blockIdx.x*blockDim.x + threadIdx.x;


  if(isp < Nvol){
    for(int dc=0;dc<12;dc++){as[dc] = a[dc*Nvol+isp];}
    for(int dc=0;dc<12;dc++){bs[dc] = b[dc*Nvol+isp];}

    {
    for(int bi=0;bi<4;bi++)
    {
      int iv = bi*4 + 0;
      for(int ai=0;ai<4;ai++)
      {

        Eigen::Map<const EigenVq > aM(&as[ai*3+0],3);
        Eigen::Map<const EigenVq > bM(&bs[bi*3+0],3);
        resab[tid*16 + iv] =  bM.dot(aM);
        iv += 1;
      }
    }
    }

  }else{
    for(int iv=0;iv<16;iv++){resab[tid*16 + iv] = 0.0;}
  }
  __syncthreads();

  ///resab --> isp --> 16 --> reduce by a factor of 4/8/16
  ////Assume thread number 32
  if(tid<16){for(int is= 1;is<16;is++){resab[ 0*16 +tid   ] += resab[is*16+tid   ];}}
  else{      for(int is=17;is<32;is++){resab[16*16 +tid-16] += resab[is*16+tid-16];}}
  __syncthreads();

  if(tid < 16){resab[ 0*16 + tid] += resab[16*16 + tid];}
  __syncthreads();

  const int it  = isp/Nsum;
  ////blockIdx.x*blockDim.x
  const unsigned long offv = Nsum/32;
  const unsigned long off0 = it*16*offv + blockIdx.x%offv;
  Complexq *f0 = &fd[off0];
  if(tid==0){for(int iv=0;iv<16;iv++)f0[iv*offv] = resab[iv];}
  ////__syncthreads();

}

namespace qlat{

inline void prodab(qlat::FermionField4dT<Complexq > &a0,qlat::FermionField4dT<Complexq > &b0, Complexq *fd){
  const qlat::Geometry &geo = a0.geo();
  const Coordinate vg = geo.total_site();
  int nt = vg[3];

  unsigned long Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  long Nsum = Nvol/Nt;

  //print0("===nthreads %8d %8d \n",qlat::qacc_num_threads(),omp_get_num_threads());

  {
  ////TIMER("Prod core 1");

  Complexq* a = (Complexq* ) &(a0.get_elem(0));
  Complexq* b = (Complexq* ) &(b0.get_elem(0));
  Complexq *f0 = &fd[0];

  //int nthreads = qlat::qacc_num_threads();
  const int nthreads = 32;
  size_t bN = (Nvol+nthreads-1)/nthreads;
  ///#pragma acc host_data use_device (a,b,f0)
  size_t bSize = nthreads*16*sizeof(Complexq);
  prodab_global<<< bN, nthreads, bSize >>>(a,b,f0, Nvol, Nsum);

  #ifdef QLAT_USE_ACC
  qacc_barrier(dummy);
  #endif
  }

  return ;

}

inline void reducefd(qlat::vector<Complexq > &fd,qlat::vector<Complexq > &NabL,qlat::vector<ga_M > &gL,const Geometry &geo,const int nvec,const Ftype facvol, inputpara& in, unsigned long bufi){
  /////const qlat::Geometry &geo = a0.geo();
  const Coordinate vg = geo.total_site();
  int nt = vg[3];

  unsigned long Nvol = geo.local_volume();
  int Nt = geo.node_site[3];
  long Nsum = Nvol/Nt;
  if(Nsum%32 !=0){print0("Assumed Nsum divice 32 == 0, %8d \n",Nsum%32);qassert(false);}
  Nsum = Nsum/32;

  qlat::vector<Complexq > reduce_sum;reduce_sum.resize((nvec*Nt)*16);
  set_zero(reduce_sum);

  int bSum = 256/nvec;
  int cutN  = 32;
  reduce_gpu2d_6(&fd[0],&reduce_sum[0],Nsum,nvec*Nt*16,  1, bSum, cutN);

  unsigned long bufN = NabL.size()/(nt*16);

  #pragma omp parallel for
  for(int op0=0;op0<nvec*Nt*16;op0++){
    int ivec = op0/(Nt*16);int op = op0%(Nt*16);
    int it = op/16; int gi = op%16;
    ////Nab[ivec*Nt*16 + it*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;

    NabL[ivec*Nt*bufN*16 + it*bufN*16 + bufi*16 + gi] += reduce_gamma(&reduce_sum[ivec*Nt*16 + it*16+0],gL[gi])/facvol;
  }

}

inline void multiplyNab_Global(const qlat::vector<Complexq > &Nab, qlat::vector<Ftype > &Mres,std::vector<int > avL, std::vector<int > bvL,const qlat::vector<Complexq > &values,const int &n_vec,const int &nt,const int nzero,const unsigned long bufN0)
{
  unsigned long bufN = avL.size();
  if(bufN == 0)return;

  int nmass = values.size()/n_vec;
  /////int nt = Nab.size()/(Aoper);

  /////Set up Mvalues
  qlat::vector<Ftype > Mvalues;Mvalues.resize(bufN*nmass*2);
  #pragma omp parallel for
  for(unsigned long bmi=0;bmi<bufN*nmass;bmi++){
    unsigned long bi = bmi/nmass;
    int mi = bmi%nmass;
    int av = avL[bi];
    int bv = bvL[bi];

    int caseab = 2;
    Ftype fac_ab = 2.0;
    if(av==bv)fac_ab = 1.0;
    if(av >= nzero and bv >= nzero ){caseab = 2;}
    if(av >= nzero and bv <  nzero ){caseab = 1;}
    if(av <  nzero and bv <  nzero ){caseab = 0;}

    Complexq la = values[av*nmass + mi];
    Complexq lb = values[bv*nmass + mi];

    unsigned long offM = bi*nmass*2 + mi*2;
    if(caseab == 2){
      Mvalues[offM+0] = 2*fac_ab*(la*lb).real();
      Mvalues[offM+1] = 2*fac_ab*(la*qconj(lb)).real();
    }
    if(caseab == 1)
    {
      Mvalues[offM+0] = 2*fac_ab*(la*lb).real();
      Mvalues[offM+1] = 0.0;
    }

    if(caseab == 0)
    {
      Mvalues[offM+0] = fac_ab*(la*lb).real();
      Mvalues[offM+1] = 0.0;
    }
  }

  int  nthreads = nt;
  long nB = bufN;

  long largeB = nt*16;if(nmass > largeB){largeB = nmass;}
  int sizeB = largeB*sizeof(Complexq);
  if(nt < 16){print0("time too short for production. \n");qassert(false);}

  multiplyNab_global<<< nB, nthreads, sizeB >>>(&Nab[0],&Mres[0],&Mvalues[0],nt,nmass,bufN0);

  #ifdef QLAT_USE_ACC
  qacc_barrier(dummy);
  #endif

}




inline qlat::vector<unsigned long > get_loop_cut(int Nx,int Ny, int Nycut, int Nxcut){
  qlat::vector<unsigned long >jobL;
  jobL.resize(Ny*Nx);
  int Nx_bound = (Nx+Nxcut-1)/Nxcut;
  int Ny_bound = (Ny+Nycut-1)/Nycut;
  long count = 0;
  for(int lyi=0;lyi<Ny_bound;lyi++)
  for(int lxi=0;lxi<Nx_bound;lxi++)
  for(int ayi=0;ayi<Nycut;ayi++)
  for(int axi=0;axi<Nxcut;axi++)
  {
    long yi = lyi*Nycut + ayi;
    long xi = lxi*Nxcut + axi;
    if(xi < Nx and yi < Ny){
      jobL[count] = yi*Nx + xi;count += 1;
    }
  }
  return jobL;
}

inline void get_map_gammaL(qlat::vector<ga_M > &g0,qlat::vector<ga_M > &gL,qlat::vector<signed char > &Gmap){
  Gmap.resize(32);
  for(int i=0;i<16;i++){
    unsigned long r0;unsigned long r1;
    unsigned long a0;unsigned long a1;
    signed char findi =-1;
    signed char sign = 1;
    g0[i].check_sum(r0,r1);
    for(int j=0;j<16;j++){
      gL[j].check_sum(a0,a1);
      if(r0==a0){
        if(findi != -1){print0("WRONG!!!!\n");}
        findi = j;
      }
      if(r0==a1){
        if(findi != -1){print0("WRONG!!!!\n");}
        findi = j;
        sign  = -1;
      }
    }
    if(findi == -1){print0("WRONG!!!! %5d \n",findi);}
    Gmap[i*2+0] = findi;
    Gmap[i*2+1] = sign;
  }
  
}

inline void get_low_rho(std::vector<qlat::FermionField4dT<Complexq > > &eigen,const qlat::vector<Complexq > &values,const int &nzero,qlat::vector<Ftype > &Mres,const qlat::Geometry &geo, inputpara& in)
{
  ////Input must be chiral vectors, eigen_chi, n_vec --> chi --> d/2 --> t,y,z,x --> c --> complex
  ////values --> massi, n_vec
  ////
  const int n_vec = eigen.size();
  const int nmass = values.size()/n_vec;
  const Coordinate vg = geo.total_site();
  const int nt = vg[3];

  Ftype facvol = std::sqrt(vg[0]*vg[1]*vg[2]);

  /////Get map list
  unsigned short Nt = geo.node_site[3];

  ga_matrices_cps   ga_cps;
  qlat::vector<ga_M > gL;gL.resize(Aoper);
  qlat::vector<ga_M > g0;g0.resize(Aoper);
  qlat::vector<ga_M > g05;g05.resize(Aoper);
  qlat::vector<ga_M > g1;g1.resize(Aoper);
  qlat::vector<ga_M > g15;g15.resize(Aoper);
  //////0 , 1, 2, 3, 4, 5, 6
  //////1-2, 1-3, 1-4, 1-5
  //////2-3, 2-4, 2-5
  //////3-4, 3-5
  //////4-5
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}
  
  ////GL

  for(int i=0;i<Aoper;i++){
     g0[i] = gL[i];
    g05[i] = (gL[5]*gL[i])*gL[5];
     g1[i] = gL[i]*gL[5];
    g15[i] = gL[5]*gL[i];
  }

  ///std::vector<std::vector<int > > Gmap;Gmap.resize(4);
  ///for(int gi=0;gi<4;gi++){Gmap[gi].resize(32);}
  qlat::vector<signed char> Gmap0;///g0
  qlat::vector<signed char> Gmap1;///g05
  qlat::vector<signed char> Gmap2;///g1
  qlat::vector<signed char> Gmap3;///g15
  get_map_gammaL(g0 ,gL, Gmap0);
  get_map_gammaL(g05,gL, Gmap1);
  get_map_gammaL(g1 ,gL, Gmap2);
  get_map_gammaL(g15,gL, Gmap3);

  cudaMemcpyToSymbol(Gmap0C, &Gmap0[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(Gmap1C, &Gmap1[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(Gmap2C, &Gmap2[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(Gmap3C, &Gmap3[0],32*sizeof(signed char),0 , cudaMemcpyHostToDevice);

  int Ncut = n_vec;

  int  noden   = qlat::get_num_node();
  long npoints = eigen[0].geo().local_volume()*12;
  double vGb_vec = npoints*noden*2.0/(1024.0*1024*1024);

  int meas = 4;int Fcount = 3 + 1;////((complex multi 6 + plus 2)/2)
  ////int meas = 4;int Fcount = 2 + 1;////Original
  double vGb     = vGb_vec*meas*Fcount;
  print0("==total Eigen %.3e Gb \n",vGb_vec*(sizeof(Complexq)/2.0)*n_vec);

  ////double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  int bufa0 = -1;int bufa1 = -1;
  int bufb0 = -1;int bufb1 = -1;
  int modeCopy = 0;
  #ifdef QLAT_USE_ACC
  size_t freeM = 0;size_t totalM = 0;
  cudaMemGetInfo(&freeM,&totalM);
  double freeD = freeM*pow(0.5,30);double totalD = totalM*pow(0.5,30);
  int Nfull = (freeD*noden/(vGb_vec*sizeof(Complexq)/2.0));
  if(n_vec < Nfull/1.3){
    Ncut = n_vec;
    //Ncut = n_vec/10;
    //Ncut = 5;
  }
  else{Ncut = Nfull/2.5;}///(or 2.5)
  print0("==n_vec %8d, Ncut %5d/%5d , Fac %.3e , free %.3e GB, total %.3e GB \n",
      n_vec,Ncut,Nfull,n_vec*1.0/Ncut,freeD, totalD);
  modeCopy = 0;
  #endif

  std::vector<qlat::FermionField4dT<Complexq > > bufE;
  if(modeCopy == 1){
    bufE.resize(2*Ncut);
    for(unsigned long iv=0;iv<bufE.size();iv++){bufE[iv].init(eigen[0].geo());}
  }

  ///int N_bound = (n_vec+Ncut-1)/Ncut;
  ///Ncut = n_vec;
  qlat::FermionField4dT<Complexq > *a0p;
  qlat::FermionField4dT<Complexq > *b0p;

  qlat::vector<Complexq > fd;fd.resize(noden*geo.local_volume()*16/32);

  ///int bufN = 1;
  long total_vol = noden*geo.local_volume();
  //unsigned long bufN = 2;
  unsigned long bufN = 2*((total_vol*12)/(nmass*16*nt));
  unsigned long bufi = 0;


  //Mres.resize(nmass*16*nt*nt);set_zero(Mres);
  qlat::vector<Ftype > MresL;
  MresL.resize(bufN*nmass*16*nt);set_zero(MresL);
  qlat::vector<Complexq > Nab;Nab.resize(16*nt);
  set_zero(Nab);

  qlat::vector<Complexq > NabL;NabL.resize(bufN*16*nt);
  set_zero(NabL);


  qlat::vector<unsigned long > jobL = get_loop_cut(n_vec,n_vec,Ncut,Ncut);
  int countrun = 0;int totrun =  0;timeval tm0,tm1,tm2;
  gettimeofday(&tm0, NULL);gettimeofday(&tm1, NULL);gettimeofday(&tm2, NULL);
  int eachrun  = 0;
  for(int jobi=0;jobi<jobL.size();jobi++){
    int avi = jobL[jobi];
    int av = avi/n_vec;
    int bv = avi%n_vec;
    if(bv > av){continue;}totrun +=1;
  }

  ////qlat::vector<Complexq > NabS;NabS.resize(16*nt);
  /////Rederive chiral forms
  std::vector<int > avL,bvL;
  avL.resize(0);bvL.resize(0);
  std::vector<int > avL_local,bvL_local;
  avL_local.resize(0);bvL_local.resize(0);



  ////Buffer index for av,bv
  //////#pragma omp parallel for
  for(int jobi=0;jobi<jobL.size();jobi++)
  {
    int avi = jobL[jobi];
    int av = avi/n_vec;
    int bv = avi%n_vec;
    if(bv > av)continue;

    avL.push_back(av);
    bvL.push_back(bv);

    if(modeCopy == 0){
      a0p = &eigen[av];
      b0p = &eigen[bv];
    }

    ////Buffer for bv
    #ifdef QLAT_USE_ACC
    if(modeCopy == 1)
    {
    TIMER("Copy memory to Device");

    if(bv >= bufb0 and bv < bufb1){b0p = &bufE[bv%Ncut];}else{
      if(bv % Ncut == 0){
        for(int iv=0;iv<Ncut;iv++){
          if(bv + iv < n_vec)cudaMemcpy(&bufE[iv].get_elem(0), &eigen[bv+iv].get_elem(0),
            npoints*sizeof(Complexq), cudaMemcpyHostToDevice);
        }
        bufb0 = (bv/Ncut)*Ncut;bufb1 = bufb0 + Ncut;
      }
      b0p = &bufE[bv%Ncut];
    }

    if(av >= bufb0 and av < bufb1){a0p = &bufE[av%Ncut];}else{
    if(av >= bufa0 and av < bufa1){a0p = &bufE[Ncut+av%Ncut];}else{
      if(av % Ncut == 0){
        for(int iv=0;iv<Ncut;iv++){
          if(av + iv < n_vec)cudaMemcpy(&bufE[Ncut + iv].get_elem(0),&eigen[av + iv].get_elem(0),
            npoints*sizeof(Complexq), cudaMemcpyHostToDevice);
        }
        bufa0 = (av/Ncut)*Ncut;bufa1 = bufa0 + Ncut;
      }
      a0p = &bufE[Ncut + av%Ncut];
    }
    }
    }
    #endif

    ////Vector reduce to Nab
    long off_fd = ((countrun%noden)*geo.local_volume()*16)/32;
    {
      TIMER("Prod core a b");
      prodab(*a0p,*b0p,&fd[off_fd]);
      //prodab(*a0p,*b0p,Nab,mapT,fd, gL);
      //{TIMER("Touch 1");touchmem(*a0p,fd);}
      //{TIMER("Touch 2");touchmem(*b0p,fd);}
    }
    countrun += 1;

    if((countrun%noden == 0) or countrun == totrun)
    {
      int nvec = noden;if(countrun%noden != 0)nvec = countrun%noden;
      {TIMER("Reduce fd");reducefd(fd,NabL,gL,geo,nvec,facvol, in, bufi);}
      //bufi += nvec;
      bufi += 1;


      int rank = qlat::get_id_node();
      if(rank < nvec){
        avL_local.push_back(avL[rank]);bvL_local.push_back(bvL[rank]);
      }

      if(avL_local.size() == bufN or countrun == totrun )
      {
        //////nvec --> it --> bufi --> 16
        {TIMER("Reduce Nab");Redistribute_all_Nt(reinterpret_cast<Ftype* > (&NabL[0]),2*NabL.size(),geo);}
        /////rank  --> nt --> bufi --> 16
        TIMER("Sum at zero node.");
        multiplyNab_Global(NabL,MresL,avL_local,bvL_local,values,n_vec,nt,nzero,bufN);
        avL_local.resize(0);bvL_local.resize(0);
        bufi = 0;
      }
      avL.resize(0);bvL.resize(0);
      set_zero(Nab);
    }

    double perc = countrun*1.0/totrun;
    eachrun  += 1;

    if(jobi%(n_vec)==0){
      gettimeofday(&tm1, NULL);
      double time0 = tm1.tv_sec - tm0.tv_sec;
      time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;

      double time1 = tm1.tv_sec - tm2.tv_sec;
      time1 += (tm1.tv_usec - tm2.tv_usec)/1000000.0;

      double flops_pers = vGb*countrun/(1.0*time0);
      double flops_pers_round = vGb*eachrun/(1.0*time1);eachrun=0;gettimeofday(&tm2, NULL);
      print0("==jobi %10d, ai %5d , bi %5d , per %.3f, use %.3e sec, %.3f Gflops, %.3f Gflops/r . \n",
        jobi,av,bv, perc, time0,flops_pers,flops_pers_round);
    }
  }

  {
    gettimeofday(&tm1, NULL);
    double time0 = tm1.tv_sec - tm0.tv_sec;
    time0 += (tm1.tv_usec - tm0.tv_usec)/1000000.0;

    double flops_pers = vGb*totrun/(1.0*time0);
    print0("==Total use %.3e sec, average %.3f Gflops . \n", time0,flops_pers);
  }

  #ifdef QLAT_USE_ACC
  qacc_barrier(dummy);
  #endif

  ///ipr --> mi --> bi -- > t0
  Mres.resize(nmass*16*nt);set_zero(Mres);
  /////long Msum = Mres.size();
  for(int mi=0;mi<nmass;mi++)
  for(int ipr=0;ipr<16;ipr++)
  {
    Ftype *res =  &Mres[(mi*16 + ipr)*nt + 0];
    for(unsigned long bi=0;bi<bufN;bi++){
    Ftype *src =  &MresL[((ipr*nmass+mi)*bufN+bi)*nt + 0];
    #pragma omp parallel for
    for(int ti=0;ti<nt;ti++){
      res[ti] += src[ti];
    }
    }
  }

  {TIMER("Final sum Mres");sum_all_size(reinterpret_cast<Ftype* > (&Mres[0]),Mres.size());}

}

}

#endif

