#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
#include "check_fun.h"
#include "utils_FFT_GPU.h"
#include "utils_shift_vecs.h"

#define TyD qlat::Complex
#define TyF qlat::Complex

int main(int argc, char* argv[])
{
  using namespace qlat;

  inputpara in;int mode_dis = -1;
  begin_Lat(&argc, &argv, in, mode_dis);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  ///int icfg  = in.icfg;
  ///int ionum = in.ionum;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  std::vector<int > nv,Nv,mv;
  geo_to_nv(geo, nv,Nv,mv);

  fft_desc_basic fd(geo);

  //int Nvec = in.nvec;
  int Nvec = fd.Nmpi;
  const int inner = 6;
  ////int ns   = in.bfac;

  std::vector<qlat::FieldM<TyD, inner> > src ;src.resize(Nvec) ;
  std::vector<qlat::FieldM<TyD, inner> > res ;res.resize(Nvec) ;
  std::vector<qlat::FieldM<TyF, inner> > srcF;srcF.resize(Nvec);
  std::vector<qlat::FieldM<TyF, inner> > resF;resF.resize(Nvec);

  for(LInt iv=0;iv<src.size();iv++){
    src[iv].init(geo);
    srcF[iv].init(geo);
    res[iv].init(geo);
    resF[iv].init(geo);
  }
  TyF* PF0 =NULL ;TyD* PD1 =NULL ;
  TyF* PFa =NULL ;
  TyD* PDa =NULL ;
  //qlat::RngState rs(qlat::get_id_node() + 134 );
  //double ini = qlat::u_rand_gen(rs);
  for(int iv=0;iv<Nvec;iv++){
  PD1 = (TyD*) (qlat::get_data(src[iv]).data());
  for(long isp=0;isp < geo.local_volume(); isp++){
    qlat::Coordinate ts = geo.coordinate_from_index(isp);
    qlat::Coordinate p  = geo.coordinate_g_from_l(ts);
    double offV0 = ((p[3]*900+p[2])*900+p[1])*900+p[0];
    for(int di=0;di<inner;di++){
        ////P1[isp*Nvec + di ] = TyD(std::cos((iv + ini+isp+di)*0.5) , (5.0/(isp+1))*(iv + ini+di)*0.1);
        TyD ref0 = TyD(offV0 + std::cos((di)*0.5) , (5.0/(1.0))*(iv + di)*0.1);
        PD1[isp*inner + di ] = ref0;
    }
  }
  
 
  }

  int biva = src.size();int civ = inner;
  bool GPU = true;if(in.GPU == 0){GPU = false;}
  shift_vec svec(fd, GPU);
  //shift_vec svec(fd, true);
  svec.set_MPI_size<TyF >(biva, civ);

  //shift_vec svec(fd, false);
  //svec.set_MPI_size<TyF >(2, 12*12);
  //ckpoint;

  long volE = geo.local_volume() * civ;
  for(LInt ni=0;ni<src.size();ni++){
    PF0 = (TyF*) qlat::get_data(srcF[ni]).data();
    PD1 = (TyD*) qlat::get_data(src[ ni]).data();
    cpy_data_thread(PF0, PD1, volE, 1, true);

    PDa = (TyD*) qlat::get_data(res[ ni]).data();
    cpy_data_thread(PDa, PD1, volE, 1, true);

  }

  std::vector<int > iDir;iDir.resize(4);
  std::vector<std::string > temL = stringtolist(in.paraI);
  if(temL.size() != 4){abort_r("Wrong temL. \n");};
  for(int i=0;i<4;i++){iDir[i] = stringtonum(temL[i]);}

  ///svec.set_MPI_size(2*civ, biva);
   
  //std::vector<qlat::vector_acc<TyF > > srcE;std::vector<qlat::vector_acc<TyF > > resE;
  //srcE.resize(biva);resE.resize(biva);

  //for(int iv=0;iv<biva;iv++){
  //  srcE[iv].resize(volE);
  //  resE[iv].resize(volE);
  //}

  //////===copy src
  //PF0 = (TyF*) &srcE[0][0];
  //PD1 = (TyD*) (qlat::get_data(src[0]).data());
  //cpy_data_thread(PF0, PD1, volE, 1);
  ////===copy src

  ////=====fft_desc shift
  //for(int vi=0;vi<in.nvec;vi++){
  //  svec.shift_Evec(srcE, resE, iDir, civ);
  //  for(int ni=0;ni<srcE.size();ni++)srcE[ni] = resE[ni];
  //}

  ////svec.shift_Evec(srcE, resE, iDir, civ);

  for(int li=0;li<in.nvec;li++){
    shift_fieldM(svec, srcF, resF, iDir);
    for(LInt ni=0;ni<srcF.size();ni++){
      PF0 = (TyF*) qlat::get_data(srcF[ni]).data();
      PFa = (TyF*) qlat::get_data(resF[ ni]).data();
      cpy_data_thread(PF0, PFa, volE, 1, true);
    }
  }

  svec.print_info();

  //////=====fft_desc shift
  //////ckpoint;

  //////=====qlat shift
  for(int li=0;li<in.nvec;li++){
  for(LInt ni=0;ni<src.size();ni++){
  for(int dir=0;dir<4;dir++){
    field_shift_dir(src[ni], src[ni], dir, iDir[dir]);
  }}}
  //////=====qlat shift

  ////===check diff

  double diff = 0.0;double sum = 0.0;double sum0 = 0.0;double sum1 = 0.0;
  double diff0 = 0.0;double diff1 = 0.0;
  qlat::vector_acc<int > iDir_acc;iDir_acc.resize(4);for(int i=0;i<4;i++){iDir_acc[i] = iDir[i]*in.nvec;}
  qlat::vector_acc<int > nv_acc;nv_acc.resize(4);for(int i=0;i<4;i++){nv_acc[i] = fd.nv[i];}
  for(LInt iv=0;iv<src.size();iv++){
    PF0 = (TyF*) qlat::get_data(resF[iv]).data();
    PD1 = (TyD*) qlat::get_data(src[ iv]).data();
    for(long isp=0;isp<geo.local_volume();isp++){
      qlat::Coordinate ts = geo.coordinate_from_index(isp);
      qlat::Coordinate p  = geo.coordinate_g_from_l(ts);
      double offV0 = ((p[3]*900+p[2])*900+p[1])*900+p[0];
      for(int i=0;i<4;i++){p[i] = (p[i] - iDir_acc[i] + nv_acc[i])%nv_acc[i];}
      //for(int i=0;i<4;i++){p[i] = (p[i] + iDir_acc[i] + nv_acc[i])%nv_acc[i];}
      double offV1 = ((p[3]*900+p[2])*900+p[1])*900+p[0];

      for(int di=0;di<inner;di++){
          long offP = isp*inner + di ;
          ////P1[isp*Nvec + di ] = TyD(std::cos((iv + ini+isp+di)*0.5) , (5.0/(isp+1))*(iv + ini+di)*0.1);
          TyD ref0 = TyD(offV0 + std::cos((di)*0.5) , (5.0/(1.0))*(iv + di)*0.1);
          TyD ref1 = TyD(offV1 + std::cos((di)*0.5) , (5.0/(1.0))*(iv + di)*0.1);

          diff += qnorm(PF0[offP] - PD1[offP]);

          diff0 += qnorm(PF0[offP] - ref1);
          diff1 += qnorm(PD1[offP] - ref1);

          sum  += qnorm(ref0)    ;
          sum0 += qnorm(PF0[offP]);
          sum1 += qnorm(PD1[offP]);
          //sum1 += qnorm(PD1[offP] - ref0);
      }
    }
  }
  sum_all_size(&diff, 1);
  sum_all_size(&diff0, 1);
  sum_all_size(&diff1, 1);
  sum_all_size(&sum, 1);
  sum_all_size(&sum0, 1);
  sum_all_size(&sum1, 1);
  print0("diff shift %.3e, diff check %.3e %.3e, sum %.3e %.3e %.3e \n", diff, diff0, diff1, sum, sum0, sum1);
  ////===check diff


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

