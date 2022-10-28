#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <qlat/vector_utils/general_funs.h>
#include <qlat/vector_utils/utils_fft_desc.h>
#include <qlat/vector_utils/utils_smear_vecs.h>

using namespace qlat;

void simple_tests()
{
  TIMER_VERBOSE("simple_structures");
  RngState rs(get_global_rng_state(), fname);
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);

  {
    //qlat::vector<double > buf;buf.resize(16);
    //for(int i=0;i<buf.size();i++){buf[i] = 0;}
    //sum_all_size(buf.data(), buf.size());
    //displayln_info(ssprintf("CHECK: sum_all_size: OK") );
    TIMER_VERBOSE("test-fft-sec-basic");
    //VectorGPUKey gkey(100*sizeof(qlat::Complex), std::string("test_buf"), true);
    //vector_gpu<char >& buf = get_vector_gpu_plan<char >(gkey);
    //buf[0] = 0;
    //buf[100*sizeof(qlat::Complex)-1] = 1;

    vector_gpu<char > buf; buf.resize(500);
    buf[0] = 0;

    displayln_info(ssprintf("CHECK: vector gpu: OK") );

    //fft_desc_basic fd(geo);
    //(void) fd;
    //////const fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
    //////size_t offv = fd.index_g_from_local(0 , 0);
    //////(void) offv;
    //displayln_info(ssprintf("CHECK: fft-sec-basic: OK") );
  }

  //{
  //  TIMER_VERBOSE("test-vector-gpu");
  //  VectorGPUKey gkey(100*sizeof(qlat::Complex), std::string("test_buf"), true);
  //  vector_gpu<char >& buf = get_vector_gpu_plan<char >(gkey);
  //  buf[0] = 0;
  //  buf[100*sizeof(qlat::Complex)-1] = 1;
  //  displayln_info(ssprintf("CHECK: vector gpu: OK") );
  //}

  //{
  //  TIMER_VERBOSE("test-move-index");
  //  move_index mv_idx;
  //  qlat::vector<qlat::Complex > buf;buf.resize(800);
  //  mv_idx.dojob(buf.data(), buf.data(), 2, 50, 4, 1,   2, true);
  //  displayln_info(ssprintf("CHECK: move index: OK") );
  //}

  //{
  //  TIMER_VERBOSE("test-shift-vec-cov");
  //  GaugeField gf;gf.init(geo);
  //  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  //  Propagator4d propS;propS.init(geo);
  //  Propagator4d propT;propT.init(geo);
  //  set_g_rand_double(propS, RngState(rs, "prop"));
  //  std::vector<double > norm(4);
  //  for(int di=0;di<4;di++){norm[di] = 0;}

  //  fft_desc_basic fd(geo);
  //  shift_vec svec(fd, true);
  //  qlat::vector_gpu<qlat::Complex > gfE;
  //  extend_links_to_vecs(gfE, gf);
  //  svec.set_gauge(qlat::get_data(gfE).data(), 4, 12);

  //  for(int di=0;di<4;di++){
  //    std::vector<int > iDir(4);for(int i=0;i<4;i++){iDir[i] = 0;}
  //    iDir[di] = 1;

  //    propT  = propS;
  //    shift_fieldM(svec, propT, propT, iDir);

  //    propT -= propS;
  //    norm[di] = qnorm(propT);
  //  }
  //  displayln_info(ssprintf("CHECK: Consistency: orig qnorm: %.10E ; shift qnorm %.10E %.10E %.10E %.10E",
  //                          qnorm(propS), norm[0],  norm[1], norm[2], norm[3]));
  //}

  //{
  //  TIMER_VERBOSE("test-rotate-vec");
  //  GaugeField gf;gf.init(geo);

  //  fft_desc_basic fd(geo);

  //  const long NVmpi = fd.mz*fd.my*fd.mx;
  //  const long Nsize = fd.Nx* fd.Ny* fd.Nz* fd.Nt * 9;

  //  qlat::vector_gpu<qlat::Complex > gauge;gauge.resize(Nsize);
  //  random_Ty(gauge.data(), gauge.size(), 1, int(qlat::u_rand_gen(rs) * 100) );

  //  qlat::vector_gpu<qlat::Complex > gfT;gfT.resize(NVmpi*Nsize);
  //  qlat::vector_gpu<qlat::Complex > gfT_buf;gfT_buf.resize(NVmpi*Nsize);
  //  for(long vi=0;vi<NVmpi;vi++){cpy_data_thread(&(gfT.data()[vi*Nsize]), gauge.data(), Nsize);}

  //  Vec_redistribute vec_rot(fd);
  //  vec_rot.reorder(gfT.data(), gfT_buf.data(), 1, 9 ,   0);

  //  double gnorm = gauge.norm().real();
  //  double rnorm = gfT.norm().real();
  //  qlat::vector_gpu<qlat::Complex > diff;diff.resize(Nsize);
  //  qacc_for(isp, Nsize, {
  //    diff[isp] = qnorm(gfT[isp] - gauge[isp]);
  //  });
  //  displayln_info(ssprintf("CHECK: Consistency: orig qnorm: %.10E ; rotate qnorm %.10E, diff qnorm %.10E",
  //                          gnorm, rnorm/gnorm, diff.norm().real()));
  //}

}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qlat-basic-structure");
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
