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
  geo.init(total_site);

  {
    vector_gpu<char > buf; buf.resize(500, false);
    buf[0] = 0;
    displayln_info(ssprintf("CHECK: vector gpu: OK") );
  }

  {
    VectorGPUKey gkey(100*sizeof(qlat::ComplexD), std::string("test_buf"), false);
    vector_gpu<qlat::ComplexD>& buf = get_vector_gpu_plan<qlat::ComplexD >(gkey);
    buf[0] = 0;
    buf[100-1] = 1;
    displayln_info(ssprintf("CHECK: vector gpu buf: OK") );
  }

  {
    TIMER_VERBOSE("test-move-index");
    move_index mv_idx;
    qlat::vector_acc<qlat::ComplexD > buf;buf.resize(800);
    mv_idx.dojob(buf.data(), buf.data(), 2, 50, 4, 1,   2, true);
    displayln_info(ssprintf("CHECK: move index: OK") );
  }

  {
    fft_desc_basic& fd = get_fft_desc_basic_plan(geo);
    size_t offv = fd.index_g_from_local(0 , 0);
    (void) offv;
    displayln_info(ssprintf("CHECK: fft-sec-basic: OK") );
  }

  {
    qlat::vector<double > buf;buf.resize(16);
    for(int i=0;i<buf.size();i++){buf[i] = 0;}
    sum_all_size(buf.data(), buf.size());
    displayln_info(ssprintf("CHECK: sum_all_size: OK") );
  }

  {
    TIMER_VERBOSE("test-shift-vec-cov");
    GaugeField gf;gf.init(geo);
    set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
    Propagator4d propS;propS.init(geo);
    Propagator4d propT;propT.init(geo);
    set_g_rand_double(propS, RngState(rs, "prop"));
    std::vector<double > norm(4);
    for(int di=0;di<4;di++){norm[di] = 0;}

    fft_desc_basic fd(geo);
    shift_vec svec(fd, true);
    qlat::vector_gpu<qlat::ComplexD > gfE;
    extend_links_to_vecs(gfE, gf);
    svec.set_gauge(qlat::get_data(gfE).data(), 4, 12);

    for(int di=0;di<4;di++){
      std::vector<int > iDir(4);for(int i=0;i<4;i++){iDir[i] = 0;}
      iDir[di] = 1;

      propT  = propS;
      shift_fieldM(svec, propT, propT, iDir);

      propT -= propS;
      norm[di] = qnorm(propT);
    }
    displayln_info(ssprintf("CHECK: Consistency: orig qnorm: %.10E ; shift qnorm %.10E %.10E %.10E %.10E",
                            qnorm(propS), norm[0],  norm[1], norm[2], norm[3]));
  }

  {
    TIMER_VERBOSE("test-rotate-vec");
    GaugeField gf;gf.init(geo);

    fft_desc_basic fd(geo);

    const Long NVmpi = fd.mz*fd.my*fd.mx;
    const Long Nsize = fd.Nx* fd.Ny* fd.Nz* fd.Nt * 9;

    qlat::vector_gpu<qlat::ComplexD > gauge;gauge.resize(Nsize);
    random_Ty(gauge.data(), gauge.size(), 1, int(qlat::u_rand_gen(rs) * 100) );

    qlat::vector_gpu<qlat::ComplexD > gfT;gfT.resize(NVmpi*Nsize);
    qlat::vector_gpu<qlat::ComplexD > gfT_buf;gfT_buf.resize(NVmpi*Nsize);
    for(Long vi=0;vi<NVmpi;vi++){cpy_data_thread(&(gfT.data()[vi*Nsize]), gauge.data(), Nsize, 1);}

    Vec_redistribute vec_rot(fd, true);
    vec_rot.reorder(gfT.data(), gfT_buf.data(), 1, 9 ,   0);

    double gnorm = gauge.norm2().real();
    double rnorm = gfT.norm2().real();
    qlat::vector_gpu<qlat::ComplexD > diff;diff.resize(Nsize);
    ComplexD* p1 = gfT.p;
    ComplexD* p2 = gauge.p;
    ComplexD* p3 = diff.p;
    qacc_for(isp, Nsize, {
      p3[isp] = qnorm(p1[isp] - p2[isp]);
    });
    displayln_info(
        ssprintf("CHECK: Consistency: orig qnorm: %.10E ; rotate qnorm %.10E, "
                 "diff qnorm %.10E",
                 gnorm, rnorm / gnorm, diff.norm2().real()));
  }

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
