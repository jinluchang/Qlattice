// utils_FFT_GPU.h
// Gen Wang
// Sep. 2021

#ifndef UTILS_FFT_GPU_H
#define UTILS_FFT_GPU_H
#pragma once

#include "general_funs.h"
#include "utils_Vec_redistribute.h"
#ifdef __QLAT_WITH_FFT_MPI__
#include "fftw3-mpi.h"
#endif
//#include <inc/cufft.h>
#ifdef QLAT_USE_ACC
#include <cufftXt.h>
#endif

/////FFT for complex to complex for 3D on CPU and GPU, 4D only on CPU

namespace qlat
{

struct FFT_Vecs{

  int dim;
  std::vector<int > nv;
  int vol;

  void *fft_dat;

  fftw_plan   plan_cpuD0,plan_cpuD1;
  fftwf_plan  plan_cpuF0,plan_cpuF1;
  #ifdef QLAT_USE_ACC
  cufftHandle plan_gpu;
  #endif

  bool GPU;
  int single_type;

  size_t datasize;
  ////size_t doffsize;
  bool flag_mem_set;
  int bsize;
  int civ;

  int istride_fac_default, idist_default;

  ////MPI fft paras
  ptrdiff_t block0, color_xyz, ranku;
  ptrdiff_t alloc_local,local_n0,local_0_start;
  /////ptrdiff_t i0, j0, k0;
  ptrdiff_t* nrank;
  MPI_Comm fft_comm;
  std::vector<size_t > MPI_para;
  size_t MPI_datasize;
  ////MPI fft paras

  template<typename Ty>
  void set_plan(std::vector<int>& nv_set, int civ_set=2, std::vector<size_t > MPI_para_set = std::vector<size_t>());

  inline void clear_plan();

  template<typename Ty>
  void do_fft(Ty* inputD, bool fftdir = true, bool dummy = true);

  ////only one mode (cpu/gpu) per FFT struct
  FFT_Vecs(bool GPU_set=true){
    (void)GPU_set;
    /////qassert(GPU == 0 or GPU == 1);
    #ifndef QLAT_USE_ACC
    GPU = false;
    #else
    GPU = GPU_set;
    #endif
    fft_dat = NULL;
    //va_dat = NULL;
    //vb_dat = NULL;
    flag_mem_set = false;
    vol = 1;
    single_type = -1;
    bsize = -1;
    ////biva = -1;
    civ = -1;
    dim = 0;nv.resize(0);

    MPI_para.resize(0);
    ////block0 = -1; color_xyz = -1; ranku=-1;
    nrank = new ptrdiff_t[10];
    istride_fac_default = 1;
    idist_default = 1;
    //fftw_mpi_init();
  }

  ~FFT_Vecs()
  {
    clear_plan();
    //fftw_mpi_cleanup();
    delete [] nrank;
  }

};

template<typename Ty>
void FFT_Vecs::set_plan(std::vector<int>& nv_set, int civ_set, std::vector<size_t > MPI_para_set)
{
  TIMERB("FFT_Vecs Set up plan");
  qassert(nv_set.size()!=0 and civ_set > 0);
  ///qassert(nv_set.size()!=0 and biva_set>0 and civ_set > 0);
  int do_clear = 0;
  if(do_clear == 0 and nv_set != nv){do_clear = 1;}
  ///if(do_clear == 0 and biva_set != biva){do_clear = 1;}
  if(do_clear == 0 and civ_set != civ){do_clear = 1;}
  if(do_clear == 0 and sizeof(Ty) != bsize){do_clear = 1;}
  if(do_clear == 0 and MPI_para_set != MPI_para){do_clear = 1;} 

  if(flag_mem_set == true and do_clear == 0){return ;}
  else{
    if(do_clear == 1){clear_plan();}

    nv  = nv_set;
    dim = nv.size();
    ////biva = biva_set;
    civ = civ_set;

    bsize = sizeof(Ty);
    if(bsize == 2*sizeof(double))single_type = 0;
    if(bsize == 2*sizeof(float ))single_type = 1;
    qassert(single_type != -1);

    vol = 1;for(int i=0;i<dim;i++){vol = vol*nv[i];}qassert(vol > 0);
    #ifdef QLAT_USE_ACC
    if(MPI_para.size() !=0 ){abort_r("GPU FFT MPI version not supported. \n");}
    #endif

    MPI_para  = MPI_para_set;
    qassert(MPI_para.size() == 3 or MPI_para.size() == 0);
    if(MPI_para.size() == 3){
      block0 = MPI_para[0];
      ////color_xyz = cxyz;ranku = ranku_set;
      color_xyz = MPI_para[1];
      ranku = MPI_para[2];
    }
  
  }

  int howmany = civ;   //// the number of transforms in fft_dat
  int istride = howmany* istride_fac_default;int idist = idist_default;
  //doffsize  = vol * howmany * bsize;
  //datasize  = doffsize;
  datasize  = vol * howmany * bsize;

  if(GPU){
  #ifdef QLAT_USE_ACC
  /////====GPU parts
  CUFFT_CALL( cufftCreate( &plan_gpu ) );
  if(dim == 4){abort_r("dim 4 on GPU not supported! \n");}

  cufftType cutype = CUFFT_Z2Z;
  if(single_type == 0)cutype = CUFFT_Z2Z;
  if(single_type == 1)cutype = CUFFT_C2C;

  CUFFT_CALL( cufftPlanMany(&plan_gpu, dim, &nv[0],
      &nv[0], istride, idist,
      &nv[0], istride, idist, cutype, howmany) );

  /////====GPU parts
  #endif
  }else{

  if(MPI_para.size() == 3)
  {
    #ifdef __QLAT_WITH_FFT_MPI__
    MPI_Comm_split(get_comm(), color_xyz, ranku, &fft_comm);
    for(int i=0;i<int(nv.size());i++){nrank[i] = nv[i];}

    if(single_type == 0){
    alloc_local = fftw_mpi_local_size_many(dim, nrank, howmany, block0, fft_comm, &local_n0, &local_0_start);
    fft_dat = fftw_malloc(alloc_local*bsize);
    plan_cpuD0 = fftw_mpi_plan_many_dft(dim,nrank, howmany, block0, block0, (fftw_complex*) fft_dat, (fftw_complex*) fft_dat,
              fft_comm, FFTW_FORWARD , FFTW_MEASURE);
    plan_cpuD1 = fftw_mpi_plan_many_dft(dim,nrank, howmany, block0, block0, (fftw_complex*) fft_dat, (fftw_complex*) fft_dat,
              fft_comm, FFTW_BACKWARD, FFTW_MEASURE);
    }

    if(single_type == 1){
    alloc_local = fftwf_mpi_local_size_many(dim, nrank, howmany, block0, fft_comm, &local_n0, &local_0_start);
    fft_dat = fftwf_malloc(alloc_local*bsize);
    plan_cpuF0 = fftwf_mpi_plan_many_dft(dim,nrank, howmany, block0, block0, (fftwf_complex*) fft_dat, (fftwf_complex*) fft_dat,
              fft_comm, FFTW_FORWARD , FFTW_MEASURE);
    plan_cpuF1 = fftwf_mpi_plan_many_dft(dim,nrank, howmany, block0, block0, (fftwf_complex*) fft_dat, (fftwf_complex*) fft_dat,
              fft_comm, FFTW_BACKWARD, FFTW_MEASURE);
    }
  
    /////each node has data vol/("Nv[0]")
    MPI_datasize = datasize/(nrank[0]/block0);
    if(local_0_start != ranku or local_n0 != block0){abort_r("fft_mpi not correct !\n");}
    #else
    abort_r("fft_mpi not set! \n");
    #endif
  }
  else{

    ////void* vb_dat = NULL;
    if(single_type == 0)fft_dat =  fftw_malloc(datasize);
    if(single_type == 1)fft_dat = fftwf_malloc(datasize);

    /////////"n" same layout for input, "istride" -- zyx(j) -- j*istride+k*idist , "idist" -- howmany(k) -- in+k*idist
    if(single_type == 0){
    plan_cpuD0 = fftw_plan_many_dft(dim,&nv[0], howmany,
        (fftw_complex*)  fft_dat,&nv[0],istride,idist,
        (fftw_complex*)  fft_dat,&nv[0],istride,idist, FFTW_FORWARD,  FFTW_MEASURE);
    plan_cpuD1 = fftw_plan_many_dft(dim,&nv[0], howmany,
        (fftw_complex*)  fft_dat,&nv[0],istride,idist,
        (fftw_complex*)  fft_dat,&nv[0],istride,idist,FFTW_BACKWARD, FFTW_MEASURE);}

    if(single_type == 1){
    plan_cpuF0 = fftwf_plan_many_dft(dim,&nv[0], howmany,
        (fftwf_complex*) fft_dat,&nv[0],istride,idist,
        (fftwf_complex*) fft_dat,&nv[0],istride,idist, FFTW_FORWARD,  FFTW_MEASURE);
    plan_cpuF1 = fftwf_plan_many_dft(dim,&nv[0], howmany,
        (fftwf_complex*) fft_dat,&nv[0],istride,idist,
        (fftwf_complex*) fft_dat,&nv[0],istride,idist,FFTW_BACKWARD, FFTW_MEASURE);}

    ////if(rotate_buf == false)
    if(fft_dat!=NULL){
      if(single_type == 0){fftw_free(fft_dat) ;}
      if(single_type == 1){fftwf_free(fft_dat);}fft_dat = NULL;}
  }////MPI verstion
  }

  ////qassert(fft_dat != NULL);
  flag_mem_set = true;

}

inline void FFT_Vecs::clear_plan()
{
  TIMERB("FFT_Vecs clear plan");
  if(flag_mem_set == false){
    return ;
  }

  if(GPU){
  #ifdef QLAT_USE_ACC
  CUFFT_CALL( cufftDestroy(plan_gpu) );
  #endif
  }else{

  if(single_type == 0){
    fftw_destroy_plan(plan_cpuD0);
    fftw_destroy_plan(plan_cpuD1);}

  if(single_type == 1){
    fftwf_destroy_plan(plan_cpuF0);
    fftwf_destroy_plan(plan_cpuF1);}
  }

  if(fft_dat != NULL){
    if(single_type == 0){fftw_free( fft_dat);}
    if(single_type == 1){fftwf_free(fft_dat);}
    fft_dat = NULL;
  }

  flag_mem_set = false;
  single_type = -1;
  vol = 1;
  bsize = -1;
  civ = -1;

  dim=0;nv.resize(0);

  MPI_para.resize(0);
  ////color_xyz = -1;ranku = -1;

}

template<typename Ty>
void FFT_Vecs::do_fft(Ty* inputD, bool fftdir, bool dummy)
{
  TIMERB("FFT excute");
  (void)dummy;
  if(flag_mem_set != true or sizeof(Ty) != bsize){
    print0("%d %d \n", int(sizeof(Ty)), int(bsize));
    abort_r("FFT_Vecs memory not set ! \n");
  }

  if(MPI_para.size() == 0)fft_dat = (void*) &inputD[0];

  if(GPU == true){
  #ifdef QLAT_USE_ACC

  //if(data_on_cpu_only){CUDA_RT_CALL( cudaMemcpy( fft_dat, src, datasize, cudaMemcpyHostToDevice   ) );}
  //else{                CUDA_RT_CALL( cudaMemcpy( fft_dat, src, datasize, cudaMemcpyDeviceToDevice ) );}

  if(fftdir == true )CUFFT_CALL( cufftXtExec( plan_gpu, fft_dat, fft_dat, CUFFT_FORWARD  ) );
  if(fftdir == false)CUFFT_CALL( cufftXtExec( plan_gpu, fft_dat, fft_dat, CUFFT_INVERSE  ) );
  if(dummy)qacc_barrier(dummy);

  //if(data_on_cpu_only){CUDA_RT_CALL( cudaMemcpy( res, fft_dat, datasize, cudaMemcpyDeviceToHost   ) );}
  //else{                CUDA_RT_CALL( cudaMemcpy( res, fft_dat, datasize, cudaMemcpyDeviceToDevice ) );}

  #endif
  }else{

  /////copy data from cpu
  ////print0("civ %d, nt %d, nz %d, ny %d, nx %d, block0 %d \n", civ, nv[0], nv[1], nv[2], nv[3], int(block0));
  ////print0("alloc_local %d, MPI_datasize %d . \n ", int(alloc_local), int(MPI_datasize/bsize));
  ////abort_r("check point. \n");
  if(MPI_para.size() != 0)cpy_data_thread((Ty*) fft_dat, inputD, MPI_datasize/bsize, GPU);

  if(MPI_para.size() != 0){
  if(single_type == 0){
  if(fftdir == true )fftw_execute( plan_cpuD0);
  if(fftdir == false)fftw_execute( plan_cpuD1);}

  if(single_type == 1){
  if(fftdir == true )fftwf_execute(plan_cpuF0);
  if(fftdir == false)fftwf_execute(plan_cpuF1);}}

  else{
  if(single_type == 0){
  if(fftdir == true )fftw_execute_dft( plan_cpuD0, (fftw_complex* ) fft_dat, (fftw_complex* ) fft_dat);
  if(fftdir == false)fftw_execute_dft( plan_cpuD1, (fftw_complex* ) fft_dat, (fftw_complex* ) fft_dat);}

  if(single_type == 1){
  if(fftdir == true )fftwf_execute_dft(plan_cpuF0, (fftwf_complex*) fft_dat, (fftwf_complex*) fft_dat);
  if(fftdir == false)fftwf_execute_dft(plan_cpuF1, (fftwf_complex*) fft_dat, (fftwf_complex*) fft_dat);}}

  if(MPI_para.size() != 0)cpy_data_thread(inputD, (Ty*) fft_dat, MPI_datasize/bsize, GPU);

  }

  if(MPI_para.size() == 0)fft_dat = NULL;

}


struct customLess{
    inline bool operator()(const std::vector<int >& a, const std::vector<int >& b) const {
    bool flag = false;
    if(a[3] < b[3]){flag = true;}
    ////c0 to be large
    if(a[3] == b[3] and a[2] > b[2]){flag = true;}
    if(a[3] == b[3] and a[2] == b[2] and a[0] < b[0]){flag = true;}
    //////need copy or not
    //if(a[3] == b[3] and a[4] < b[4]){flag = true;} 
    //////c0 to be large
    //if(a[3] == b[3] and a[4] == b[4] and a[2] > b[2]){flag = true;} 
    //if(a[3] == b[3] and a[4] == b[4] and a[2] == b[2] and a[0] < b[0]){flag = true;} 
    return flag;}
};


inline std::vector<int > get_factor_jobs(int nvec, int civ, int N0=-1, int N1=-1, int maxN0 = 16, int dataB=1)
{
  if(nvec < 0 or civ < 0 or maxN0 <= 0){abort_r("Vectors wrong. \n");}
  if(N1 != -1){if(N0 < N1 or N0%N1 != 0){abort_r("Vectors wrong. \n");}}
  std::vector<int > Nl;
  if(N1 != -1){Nl.resize(2);Nl[0] = N0;Nl[1] = N1;}
  if(N1 == -1){Nl.resize(1);Nl[0] = N0;}
  int total = nvec * civ;
  std::vector<std::vector<int > > jobL;
  for(LInt i=0;i<Nl.size();i++)
  {
    int Ne = Nl[i];
    int maxN = maxN0;
    //if(N1 != -1){maxN = maxN * (Nl[0]/Nl[i]);}
    //if(N1 != -1){maxN = maxN;}

    int Bmax = (total+Ne-1)/Ne;
    ////print0("N %d, maxN %d, Bmax %d, civ %d \n", Ne, maxN, Bmax, civ);

    for(int bz = 1;bz <= Bmax ; bz++)
    for(int cz = 1;cz <= maxN; cz++)
    {
      if( bz*cz <= maxN and civ%cz == 0){
        std::vector<int > tmp;tmp.resize(5);
        int one_job = bz*Ne*cz;
        int count_job = (total + one_job -1)/one_job;
        int cost = count_job*bz*Ne*cz;
        int NEED_COPY = 1;
        if((dataB*(civ/cz))%(Ne) == 0 ){NEED_COPY = 0;}
        tmp[0] = i;tmp[1] = bz;tmp[2] = cz;tmp[3] = cost;tmp[4] = NEED_COPY;
        jobL.push_back(tmp);
      }
    }
  }
  if(jobL.size() == 0){abort_r("input errors for schedule! \n");}

  std::sort(jobL.begin(), jobL.end(), customLess());
  ///int cost = jobL0[0][2];
  ///for(int i=0;i<jobL0.size();i++){if(jobL0[i][2] == cost)}
  return jobL[0];
  ////b0 = jobL[0][0];c0 = jobL[0][1];

}

struct fft_schedule{
  fft_desc_basic fd;
  int default_MPI;int enable_MPI;
  Rotate_vecs rot;
  FFT_Vecs fft;
  FFT_Vecs fft_gpu;
  bool GPU;
  std::vector<int > dimN;
  std::vector<int > dimN0, dimN1;
  int dim, nvec, civ, N_extra;
  int b0, c0;
  int maxN;
  bool flag_mem_set;
  move_index mv_civ;
  int NEED_COPY;
  int dataB;
  int bsize;

  fft_schedule(fft_desc_basic& fd_set, bool GPU_set=true):fd(),rot(fd_set,GPU_set),fft(GPU_set),fft_gpu(GPU_set)
  {
    #ifndef QLAT_USE_ACC
    GPU = false;
    #else
    GPU = GPU_set;
    #endif
    copy_fft_desc(fd, fd_set);

    enable_MPI = -1;
    dimN.resize(0);dim = -1;nvec=-1;civ=-1;b0=-1;c0=-1;N_extra=-1;
    flag_mem_set = false;
    maxN = -1;
    NEED_COPY = 1;
    dataB = 1;
    bsize = -1;
  }

  ////-1 auto choose MPI and civ
  ////-2, nvec, civ with MPI, -3 nvec, civ without MPI 
  template<typename Ty>
  void set_mem(const int nvec_set, const int civ_set, const std::vector<int >& dimN_set=std::vector<int >(), int default_MPI_set = -1 , int maxN_set=16, int dataB_set=1)
  {
    TIMERA("fft_schedule set mem");
    if(nvec == nvec_set and civ == civ_set and default_MPI == default_MPI_set and maxN == maxN_set and dataB == dataB_set){
      if(dimN_set.size() == 0 and dimN.size() == 0){abort_r("Dimension not set. \n");}
      if(dimN_set.size() == 0  or dimN == dimN_set){return ;}
    }
    bsize = sizeof(Ty);
    maxN = maxN_set;dataB = dataB_set;
    if(maxN < 1 or dataB < 1){abort_r("maxN wrong. \n");}
    if(dimN_set.size() !=0 ){dimN = dimN_set;} dim = dimN.size();
    if(dim < 3 or dim > 4  ){abort_r("Dimension not supported. \n");}

    civ = civ_set;nvec = nvec_set;default_MPI= default_MPI_set;
    if(!(default_MPI >= -3 and default_MPI < 2)){abort_r("Wrong default_MPI ! \n");};
    ////long total  = nvec*civ;
    N_extra = -1;int N0 = -1;int N1 = -1;

    /////default vector options
    if(dim == 3){N0 = fd.mz * fd.my * fd.mx;N1 = 1;}
    if(dim == 4){N0 = fd.mt * fd.mz * fd.my * fd.mx;N1 = fd.mz * fd.my * fd.mx;}

    /////====Set up b0, c0
    std::vector<int > job(5);

    if(default_MPI == -2 or default_MPI == -3)
    {
      enable_MPI = 0;
      if(!GPU){if(default_MPI == -2){enable_MPI = 1;}}
      if(enable_MPI==0)job[1] = (nvec+N0-1)/N0;
      if(enable_MPI==1)job[1] = (nvec+N1-1)/N1;
      job[2] = civ;
    }else{
      std::vector<int > job0 = get_factor_jobs(nvec, civ, N0, -1, maxN, dataB);
      std::vector<int > job1 = get_factor_jobs(nvec, civ, N1, -1, maxN, dataB);
      ///////small rotation for 3D not defined, need define it from start;
        
      if(default_MPI == -1 and dim == 3 and fd.my * fd.mx != 1){enable_MPI = 0;}

      if(enable_MPI == -1){
        if(GPU){enable_MPI = 0;job = job0;}else{
        if(default_MPI == -1){job = get_factor_jobs(nvec, civ, N0, N1, maxN, dataB); enable_MPI = job[0];}
        else{
          enable_MPI = default_MPI;
        }}
      }
      if(default_MPI >= -1){if(enable_MPI == 0){job = job0;}if(enable_MPI == 1){job = job1;}}
    }
    ////print0("job size %d \n", int(job.size()));
    if(enable_MPI == 1 and dim == 3 and fd.my * fd.mx != 1){abort_r("mode not supported ! \n");}
    if(!(enable_MPI == 0 or enable_MPI == 1)){abort_r("enable_MPI not set yes! \n");};

    b0 = job[1];c0 = job[2];///////NEED_COPY = job[4];
    /////====Set up b0, c0

    #ifndef __QLAT_WITH_FFT_MPI__
    enable_MPI = 0;
    #endif

    if(GPU){qassert(enable_MPI == 0);}
    NEED_COPY = 0;
    if(civ%(N_extra*c0) != 0 or GPU){NEED_COPY = 1;}

    int mode_rot = -2;
    if(enable_MPI == 0){N_extra = N0;if(dim == 3){mode_rot =  0;}if(dim == 4){mode_rot = 1;}}
    if(enable_MPI == 1){N_extra = N1;if(dim == 3){mode_rot = -1;}if(dim == 4){mode_rot = 0;}}

    //print_info();
    //print_mem_info("before mem fft");
    ///ckpoint;
    rot.set_mem<Ty>(b0, c0, mode_rot);

    if(enable_MPI == 0){
      if(GPU == false){
        fft.set_plan<Ty >(dimN, c0);
      }
      if(GPU == true ){
        if(dim == 3){fft.set_plan<Ty >(dimN, c0);}
        if(dim == 4){
          dimN0.resize(3);dimN1.resize(1);
          for(int i=0;i<3;i++){dimN0[i] = dimN[i + 1];}
          dimN1[0] = dimN[0];
          fft.set_plan<Ty >(dimN0, c0);
          fft_gpu.set_plan<Ty >(dimN1, fd.nz*fd.ny*fd.nx * c0);
      }}
    }
    if(enable_MPI == 1){
      std::vector<size_t > MPI_para(3);
      if(dim == 3){
      MPI_para[0] = fd.Nz;
      MPI_para[1] = fd.init; ////colorxyz for MPI
      MPI_para[2] = fd.iniz; ////initial point of the MPI FFT
      }
      
      if(dim == 4){
      MPI_para[0] = fd.Nt; 
      MPI_para[1] = (fd.iniz*fd.ny + fd.iniy)*fd.nx + fd.inix; ////colorxyz for MPI
      MPI_para[2] = fd.init; ////initial point of the MPI FFT 
      }

      fft.set_plan<Ty >(dimN, c0, MPI_para);
    }
    //print_mem_info("after mem fft");

    flag_mem_set = true;
  }

  void print_info()
  {
    print0("==Jobs dim %d, ", dim);
    for(int di=0;di<dim;di++){print0("%d ", dimN[di]);}
    print0(", nvec %d, civ %d, bsize %d. \n", nvec, civ, bsize);
    print0("==N_extra %d, enable_MPI %d, GPU %d, b0 %d, c0 %d, mode_rot %d, need copy %d. \n", 
      N_extra , enable_MPI, int(GPU), b0, c0, rot.mode, NEED_COPY );
    fflush_MPI();
  }

  template<typename Ty>
  void call_fft(Ty* src, bool fftdir=true)
  {
    TIMERB("schedule call FFT");
    Ty* tsrc;
    bool flag_GPU_4d = false;
    if(GPU == true and dim == 4){flag_GPU_4d = true;}
    if(!flag_GPU_4d){
    if(dim == 3){for(int ti=0;ti<fd.Nt;ti++){
      tsrc=&src[ti * (rot.vol_buf/fd.Nt)*c0];fft.do_fft(tsrc, fftdir, false);
    }
    qacc_barrier(dummy);
    }
    if(dim == 4){tsrc=src;fft.do_fft(tsrc, fftdir);}
    ////if(dim==4){for(int ti=0;ti<fd.nt;ti++){tsrc=&src[ti * (fd.nz*fd.ny*fd.nx) *c0]; fft.do_fft(tsrc, fftdir, false);}}
    }
    
    if(flag_GPU_4d){
      for(int ti=0;ti<fd.nt;ti++){tsrc=&src[ti * (fd.nz*fd.ny*fd.nx) *c0]; fft.do_fft(tsrc, fftdir, false);}
      qacc_barrier(dummy);
      tsrc=src;fft_gpu.do_fft(tsrc, fftdir);
    }
    tsrc = NULL;
  }

  template<typename Ty>
  void dojob(std::vector<Ty* >& data, bool fftdir=true, int civ_set = -1 ){
    TIMERC("fft dojob");
    int bN = data.size();
    int each = b0*N_extra*c0;

    ////update nvec, civ
    nvec = bN;if(civ_set != -1)civ = civ_set;

    /////if(civ_set != -1){set_mem<Ty>(bN, civ);}
    if(!flag_mem_set or b0 < 0 or c0 < 0){abort_r("Memory not set for fft_schedule! \n");}
    if(rot.flag_mem_set != true or fft.flag_mem_set != true){abort_r("Memory not set for fft_schedule! \n");}
    if(civ%c0 != 0){abort_r("civ cannot by divided by c0 ! \n");}

    int cN = civ/c0;
    /////print0("bN %d, cN %d, civ %d, c0 %d, vol %d \n", bN, cN, civ, c0, int(fd.noden));
    /////rotate civ
    if(civ != c0){
      TIMERB("fft mem reorder");
      //for(int bi=0;bi<bN;bi++)mv_civ.dojob((char*) &data[bi][0],(char*) &data[bi][0], 1, cN, fd.noden, 1, c0*sizeof(Ty));
      for(int bi=0;bi<bN;bi++)mv_civ.dojob((Ty*) &data[bi][0],(Ty*) &data[bi][0], 1, cN, fd.noden, 1, c0, GPU);
    }

    NEED_COPY = 0;
    if(civ%(N_extra*c0) != 0 or GPU){NEED_COPY = 1;}

    //////if(NEED_COPY == 0){if(civ%(N_extra*c0) != 0)abort_r("set ups wrong. !\n");}
    //////NEED_COPY = 1;
    //////if(GPU)NEED_COPY = 1;

    ////print0("==buf size %zu \n", rot.Bsize);
    int jobN = (bN*civ + each-1) /each;
    int bm = 0;
    int cm = 0;
    Ty* src;Ty* dp;
    for(int ji=0;ji<jobN;ji++)
    {
      /////print0("jobN %d, bm %d, cm %d. \n", jobN, bm, cm);
      ////Ty* tsrc = NULL;
      if(bm == bN){break;}
      int bma = bm;int cma = cm;
      ////copy result to buf
      if(NEED_COPY == 1){
        src = (Ty*) rot.src;

        for(int bi=0;bi<b0*N_extra;bi++){
          dp = &(data[bma][cma*fd.noden]);
          cpy_data_thread(&src[bi*fd.noden*c0], dp, fd.noden*c0, GPU , false);
          cma += c0;if(cma == civ){cma=0;bma+=1;if(bma == bN)break;}
        }
        qacc_barrier(dummy);

        /////rotate data for FFT
        rot.reorder<Ty >(true , src);
        ////for(int bi=0;bi<b0;bi++)fft.do_fft(&src[bi*fd.noden*c0], fftdir);
        for(int bi=0;bi<b0;bi++){
          call_fft(&src[bi*rot.vol_buf*c0], fftdir);
        }
        rot.reorder<Ty >(false, src);

        for(int bi=0;bi<b0*N_extra;bi++){
          dp = &(data[bm ][cm *fd.noden]);
          cpy_data_thread(dp, &src[bi*fd.noden*c0], fd.noden*c0, GPU , false);
          cm  += c0;if(cm  == civ){cm =0;bm +=1;if(bm  == bN)break;}
        }
        qacc_barrier(dummy);
      }

      if(NEED_COPY == 0){
        for(int bi=0;bi<b0;bi++){
          src = &(data[bm ][cm * fd.noden]);
          rot.reorder<Ty >(true , src);
          call_fft(src, fftdir);
          rot.reorder<Ty >(false , src);
          cm += N_extra*c0;
          if(cm  == civ){cm =0;bm +=1;if(bm  == bN)break;}
        }
      }
    }
    src = NULL;dp = NULL;

    if(civ != c0){
      TIMERB("fft mem reorder");
      //for(int bi=0;bi<bN;bi++)mv_civ.dojob((char*) &data[bi][0],(char*) &data[bi][0], 1, cN, fd.noden, 0, c0*sizeof(Ty));
      for(int bi=0;bi<bN;bi++)mv_civ.dojob((Ty*) &data[bi][0],(Ty*) &data[bi][0], 1, cN, fd.noden, 0, c0, GPU);
    }
  }

  void clear_mem(){
    if(flag_mem_set == false){return;}
    default_MPI = -1;enable_MPI = -1;
    
    dimN.resize(0);
    dim = -1;nvec=-1;civ=-1;b0=-1;c0=-1;N_extra=-1;
    maxN = -1;

    fft.clear_plan();
    rot.clear_mem();
    flag_mem_set = false;
    NEED_COPY = 1;
  }

  ~fft_schedule(){
    clear_mem();
  }

};

template <class Ty, int civ>
void fft_fieldM_test(std::vector<qlat::FieldM<Ty, civ> > &src, bool fftdir=true, bool GPU=true , bool fft4D = false)
{
  TIMERD("fft fieldM");
  ////need bfac job create
  ////need rotation of data or buf memory

  if(src.size() < 1)return;
  Geometry& geo = src[0].geo();
  fft_desc_basic fd(geo);

  std::vector<int > nv,Nv,mv;
  geo_to_nv(geo, nv,Nv,mv);

  std::vector<int > dimN;int nvec = src.size();
  std::vector<Ty* > data;data.resize(nvec);
  for(int si=0;si<nvec;si++){data[si] = (Ty*) qlat::get_data(src[si]).data();}
  /////3D
  if(!fft4D){dimN.resize(3);dimN[0] = nv[2];dimN[1] = nv[1];dimN[2] = nv[0];}
  if( fft4D){dimN.resize(4);dimN[0] = nv[3];dimN[1] = nv[2];dimN[2] = nv[1];dimN[3] = nv[0];}

  fft_schedule fft(fd, GPU);
  fft.set_mem<Ty >(nvec, civ, dimN, -1 );
  fft.print_info();
  fft.dojob(data, fftdir);

}

template <class Ty, int civ>
void fft_fieldM(fft_schedule& fft, std::vector<qlat::FieldM<Ty, civ> >& src, bool fftdir=true)
{
  #if PRINT_TIMER>4
  TIMER_FLOPS("fft fieldM");
  timer.flops += src.size() * get_data(src[0]).data_size()/(sizeof(Ty)) * 64;
  #else
  TIMER("fft fieldM");
  #endif

  if(src.size() < 1)return;
  int nvec = src.size();
  std::vector<Ty* > data;data.resize(nvec);
  for(int si=0;si<nvec;si++){data[si] = (Ty*) qlat::get_data(src[si]).data();}

  fft.dojob(data, fftdir);
}

///Failed with pair
struct FFTGPUPlanKey {
  Geometry geo;
  bool GPU;
  bool fft4D;
  int nvec ;
  int civ;
  DATA_TYPE prec;
};

struct fft_gpu_copy{
  //void* fftP;
  //fft_gpu_copy(){fftP = NULL;}
  //~fft_gpu_copy(){if(fftP != NULL){delete ((fft_schedule*) fftP); fftP = NULL;}}
  fft_schedule* fftP;
  DATA_TYPE prec;
  bool is_copy;  // do not free memory if is_copy=true

  fft_gpu_copy(){fftP = NULL;is_copy = false;prec = Complex_TYPE;}
  fft_gpu_copy(const fft_gpu_copy& fft) 
  {
    #ifndef QLAT_USE_ACC
    qassert(false);
    #endif
    is_copy = true;
    fftP = fft.fftP;
    prec = fft.prec;
  }
  fft_gpu_copy(fft_gpu_copy&& fft) noexcept
  {
    is_copy = fft.is_copy;
    fftP = fft.fftP;
    prec = fft.prec;
    fft.is_copy = true;
  }

  qacc void swap(fft_gpu_copy& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    fft_schedule* tmp = fftP;
    fftP   = x.fftP;
    x.fftP = tmp;
    DATA_TYPE tmp_prec = prec;
    prec   = x.prec;
    x.prec = tmp_prec;
  }
  //
  const fft_gpu_copy& operator=(const fft_gpu_copy& fft)
  {
    qassert(not is_copy);
    if(fftP != NULL){delete (fftP); fftP = NULL;}

    prec = fft.prec;
    fftP = new fft_schedule(fft.fftP->fd, fft.fftP->GPU);
    int nvec = fft.fftP->nvec;
    int civ  = fft.fftP->civ;
    std::vector<int > dimN = fft.fftP->dimN;

    if(prec == Complex_TYPE ){     fftP->set_mem<Complex  >(nvec, civ, dimN, -1 );}
    else if(prec == ComplexF_TYPE){fftP->set_mem<ComplexF >(nvec, civ, dimN, -1 );}
    else{print0("Only Complex and ComplexF supported for fft on GPU! \n");qassert(false);}
    ///fft.fftP->print_info();
    ///fftP->print_info();

    return *this;
  }

  ~fft_gpu_copy(){if(fftP != NULL and is_copy == false){delete (fftP); fftP = NULL;}}
};

inline bool operator<(const FFTGPUPlanKey& x, const FFTGPUPlanKey& y)
{
  if(x.GPU   < y.GPU  ){return true;}
  if(y.GPU   < x.GPU  ){return false;}

  if(x.fft4D < y.fft4D){return true;}
  if(y.fft4D < x.fft4D){return false;}

  if(x.nvec  < y.nvec ){return true;}
  if(y.nvec  < x.nvec ){return false;}

  if(x.civ   < y.civ  ){return true;}
  if(y.civ   < x.civ  ){return false;}

  if(x.prec  < y.prec ){return true;}
  if(y.prec  < x.prec ){return false;}

  return false;
}


inline fft_gpu_copy make_fft_gpu_plan(const Geometry& geo, int nvec, int civ , bool GPU, bool fft4D , const DATA_TYPE prec)
{
  TIMER_VERBOSE("make_fft_gpu_plan");

  fft_desc_basic fd(geo);

  std::vector<int > nv,Nv,mv;
  geo_to_nv(geo, nv,Nv,mv);

  std::vector<int > dimN;
  /////3D
  if(!fft4D){dimN.resize(3);dimN[0] = nv[2];dimN[1] = nv[1];dimN[2] = nv[0];}
  if( fft4D){dimN.resize(4);dimN[0] = nv[3];dimN[1] = nv[2];dimN[2] = nv[1];dimN[3] = nv[0];}

  fft_gpu_copy ft;
  ft.fftP = new fft_schedule(fd, GPU);
  ft.prec = prec;

  if(prec == Complex_TYPE ){     ft.fftP->set_mem<Complex  >(nvec, civ, dimN, -1 );}
  else if(prec == ComplexF_TYPE){ft.fftP->set_mem<ComplexF >(nvec, civ, dimN, -1 );}
  else{print0("Only Complex and ComplexF supported for fft on GPU! \n");qassert(false);}
  ft.fftP->print_info();

  ///int nvec = src.size();
  //std::vector<Ty* > data;data.resize(nvec);
  //for(int si=0;si<nvec;si++){data[si] = (Ty*) qlat::get_data(src[si]).data();}
  //fft.dojob(data, fftdir);

  return ft;
}

inline fft_gpu_copy make_fft_gpu_plan(const FFTGPUPlanKey& fkey)
{
  return make_fft_gpu_plan(fkey.geo, fkey.nvec, fkey.civ, fkey.GPU, fkey.fft4D, fkey.prec);
}

inline Cache<FFTGPUPlanKey, fft_gpu_copy >& get_fft_gpu_plan_cache()
{
  static Cache<FFTGPUPlanKey, fft_gpu_copy > cache("FFTGPUPlanCache", 16);
  return cache;
}

inline const fft_gpu_copy& get_fft_gpu_plan(const FFTGPUPlanKey& fkey)
{
  if (!get_fft_gpu_plan_cache().has(fkey)) {
    get_fft_gpu_plan_cache()[fkey] = make_fft_gpu_plan(fkey);
  }
  //get_fft_gpu_plan_cache()[fkey].set();
  //get_fft_gpu_plan_cache()[fkey].fftP->print_info();
  return get_fft_gpu_plan_cache()[fkey];
}


template<typename Ty, int civ>
inline FFTGPUPlanKey get_fft_gpu_plan_key(std::vector<qlat::FieldM<Ty, civ> >& src, bool fft4d=false)
{
  qassert(src.size() > 0);
  FFTGPUPlanKey fkey;
  fkey.geo = src[0].geo();
  fkey.GPU = true;
  fkey.nvec = src.size();
  fkey.civ = civ;
  fkey.prec = get_data_type<Ty >();

  fkey.fft4D = fft4d;
  return fkey;
}

bool check_fft_mode(const int nfft, const Geometry& geo, const bool fft4d)
{
  bool use_qlat = false;
  std::vector<int > nv, Nv, mv;
  geo_to_nv(geo, nv, Nv, mv);
  if(fft4d == true ){
    if(mv[0] * mv[1] * mv[2] > nfft)
    {use_qlat = true;}
  }
  if(fft4d == false){
    if(mv[0] * mv[1] != 1 and (mv[0] * mv[1] * mv[2] > nfft))
    {use_qlat = true;}
  }

  #ifdef QLAT_USE_ACC
  use_qlat = false;
  #endif
  return use_qlat;
}

template <class Ty>
void fft_fieldM(std::vector<Ty* >& data, int civ, const Geometry& geo, bool fftdir=true, bool fft4d = false)
{
  if(data.size() < 1)return;

  #if PRINT_TIMER>4
  TIMER_FLOPS("fft fieldM");
  timer.flops += data.size() * geo.local_volume() * civ * 64;
  #else
  TIMER("fft fieldM");
  #endif

  FFTGPUPlanKey fkey;
  fkey.geo  = geo;
  fkey.GPU  = true;
  fkey.nvec = data.size();
  fkey.civ = civ;
  fkey.prec = get_data_type<Ty >();
  fkey.fft4D = fft4d;
  ////std::vector<Ty* > data;data.resize(nvec);
  ////for(int si=0;si<nvec;si++){data[si] = (Ty*) qlat::get_data(src[si]).data();}
  /////int nvec = data.size();
  get_fft_gpu_plan(fkey).fftP->dojob(data, fftdir);
}

template <class Ty>
void fft_fieldM(Ty* src, int nvec, int civ, const Geometry& geo, bool fftdir=true, bool fft4d = false)
{
  std::vector<Ty* > data;data.resize(nvec);
  const long offD = geo.local_volume() * civ;
  for(int iv=0;iv<nvec;iv++){data[iv] = &src[iv*offD];}
  fft_fieldM<Ty >(data, civ, geo, fftdir, fft4d);
}

template <class Ty, int civ>
void fft_fieldM(std::vector<qlat::FieldM<Ty, civ> >& src, bool fftdir=true, bool fft4d = false)
{
  if(src.size() < 1)return;

  int nfft = src.size() * civ;
  const Geometry& geo = src[0].geo();
  bool use_qlat = check_fft_mode(nfft, geo, fft4d);
  if(use_qlat){
    TIMER("fft_complex_field_dir fieldM");
    for(unsigned int i=0;i<src.size();i++)
    {
      qlat::FieldM<Ty, civ> ft;
      int ndir = 3;if(fft4d){ndir = 4;}
      for (int k = 0; k < ndir; k++) {
        ft = src[i];
        fft_complex_field_dir(src[i], ft, k, fftdir);
      }
    }
    return ;
  }

  {
  std::vector<Ty* > data;data.resize(src.size());
  for(unsigned int si=0;si<src.size();si++){data[si] = (Ty*) qlat::get_data(src[si]).data();}
  fft_fieldM<Ty >(data, civ, geo, fftdir, fft4d);
  }
}

template<class M>
void fft_fieldM(std::vector<Handle<qlat::Field<M> > >& src, bool fftdir=true, bool fft4d = false)
{
  if(src.size() < 1)return;
  bool is_double     = get_data_type_is_double<M >();
  DATA_TYPE prec = Complex_TYPE;int civ = 1;
  const Geometry& geo = src[0]().geo();
  if( is_double){prec = Complex_TYPE ; civ = geo.multiplicity * sizeof(M)/sizeof(Complex ); }
  if(!is_double){prec = ComplexF_TYPE; civ = geo.multiplicity * sizeof(M)/sizeof(ComplexF); }

  int nfft = src.size() * civ;
  bool use_qlat = check_fft_mode(nfft, geo, fft4d);
  if(use_qlat){
    TIMER("fft_complex_field_dir fieldM");
    for(unsigned int i=0;i<src.size();i++)
    {
      qlat::Field<M > ft;
      int ndir = 3;if(fft4d){ndir = 4;}
      for (int k = 0; k < ndir; k++) {
        ft = src[i]();
        fft_complex_field_dir(src[i](), ft, k, fftdir);
      }
    }
    return ;
  }

  {
  #if PRINT_TIMER>4
  TIMER_FLOPS("fft fieldM");
  timer.flops += src.size() * get_data(src[0]).data_size()/(sizeof(M)) * 64;
  #else
  TIMER("fft fieldM");
  #endif

  FFTGPUPlanKey fkey;
  fkey.geo = geo;
  fkey.GPU = true;
  fkey.nvec = src.size();
  fkey.civ = civ;
  fkey.prec = prec;
  fkey.fft4D = fft4d;

  int nvec = src.size();
  if( is_double){
    std::vector<Complex* > data;data.resize(nvec);
    for(int si=0;si<nvec;si++){data[si] = (Complex*) qlat::get_data(src[si]()).data();}
    get_fft_gpu_plan(fkey).fftP->dojob(data, fftdir);}

  if(!is_double){
    std::vector<ComplexF* > data;data.resize(nvec);
    for(int si=0;si<nvec;si++){data[si] = (ComplexF*) qlat::get_data(src[si]()).data();}
    get_fft_gpu_plan(fkey).fftP->dojob(data, fftdir);}
  }

}

template<typename Ty>
void FFT_vecs_corr(qlat::vector_gpu<Ty >& src, std::vector<qlat::FieldM<Ty, 1> >& FFT_data, int off=0)
{
  const Geometry& geo = FFT_data[0].geo();
  const long volume = geo.local_volume();
  const long Sdata = src.size();const int nvecs = Sdata/volume;

  fft_fieldM<Ty >(src.data(), nvecs, 1, geo, true);
  for(int iv=0;iv<nvecs;iv++)
  {
    qlat::FieldM<Ty, 1>& buf = FFT_data[off*nvecs + iv];qassert(buf.initialized);
    Ty* s0 = &src[iv*volume];
    Ty* r0 = (Ty*) qlat::get_data(buf).data();
    cpy_data_thread(r0, s0, volume, 1, true);
  }
}


}

#endif
