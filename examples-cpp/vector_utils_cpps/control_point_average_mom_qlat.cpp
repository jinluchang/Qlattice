#include <sys/sysinfo.h>
#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_eigen_ov.h"
#include "utils_check_fun.h"
#include "utils_smear_vecs.h"
#include "utils_lms_funs.h"

using namespace qlat;

inline void write_corr_zero(corr_dat<Ftype >& res, qlat::vector_gpu<Complexq >& FFT_data, const long off,
    std::vector<bool >& Bsite, std::vector<long >& Zsite)
{
  ////qmessage("write data %d \n", int(off));
  const int nmass = res.key_T[4];
  const int nt    = res.key_T[5];
  Complexq* Pres = (Complexq*) (&res.dat[off]);
  qassert(int(Bsite.size()) == nt);
  qassert(int(Zsite.size()) == nt);
  const int nvec = 32*nmass;
  qlat::vector_acc<Complexq > buf_gpu;buf_gpu.resize(nvec * nt);
  qlat::set_zero(buf_gpu);
  Complexq* Pbuf = (Complexq*) (qlat::get_data(buf_gpu).data());
  Complexq* PFFT = (Complexq*) (qlat::get_data(FFT_data).data());
  qlat::vector_acc<long> PZsite;PZsite.resize(Zsite.size());
  for(Long t=0;t<PZsite.size();t++){PZsite[t] = Zsite[t];}

  ////sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", nsource, 2, 3, 32, int(massL.size()), in.nt, 2);
  ////long corr_zero_off = (((si + npoint)*2 + sm)*3 + LHF) * 32 * nmass * in.nt * 2;
  for(int t=0;t<nt;t++){if(Bsite[t]){
  qacc_for(iv, nvec,{
    Pbuf[iv*nt + t] += PFFT[ PZsite[t]*nvec + iv ] ;
  });
  }}

  for(int t=0;t<nt;t++){if(Bsite[t]){
  qthread_for(iv, nvec,{
      Pres[iv*nt + t] += Pbuf[iv*nt + t];
  });}}

}

inline std::vector<Coordinate > read_positions(corr_dat<Ftype >& mom_info){
  std::vector<Coordinate > src_pos;
  const int Ninfo = mom_info.INFOA.size();
  bool short_info = true;
  for(int i=0;i<Ninfo;i++)
  {
    if(mom_info.INFOA[i].find("Positions") != std::string::npos)
    {
      int ini = mom_info.INFOA[i].find("Positions") + std::string("Positions").size();
      std::string str5 = mom_info.INFOA[i].substr(ini);     // get from "live" to the end
      //qmessage("test %s \n", str5.c_str());
      std::vector<Coordinate > srcs = string_to_Coordinates(str5);
      for(unsigned int j=0;j<srcs.size();j++)
      {
        src_pos.push_back(srcs[j]);
      }
      short_info = false;
    }
  }
  if(short_info == true){
    src_pos = string_to_Coordinates(mom_info.INFO_LIST);
  }

  return src_pos;
}

int main(int argc, char* argv[])
{

  inputpara in;begin_Lat(&argc, &argv, in);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  int icfg  = in.icfg;

  {

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  print_mem_info();
  fflush_MPI();

  int mass_group = 1;
  in.find_para(std::string("mass_group"), mass_group);

  std::string flag_add = std::string("");
  in.find_para(std::string("flag_add"), flag_add);
  std::string output_vecs = std::string("NONE");
  in.find_para(std::string("output_vecs"), output_vecs);
  if(output_vecs == std::string("NONE")){
    output_vecs = in.output_vec;
  }

  std::string mom_shift_ = std::string("0 0 0 0");
  in.find_para(std::string("mom_shift"), mom_shift_);
  std::vector<Coordinate > mom_shiftL;
  {
    std::vector<std::string > Li = stringtolist(mom_shift_);
    Qassert(Li.size() % 4 == 0);
    const int Nmom = Li.size() / 4;
    mom_shiftL.resize(Nmom);
    for(int momi=0;momi<Nmom;momi++)
    {
      for(int i=0;i<4;i++){mom_shiftL[momi][i] = stringtonum(Li[momi*4 + i]);}
    }
  }

  int mom_shifti = 0;
  in.find_para(std::string("mom_shifti"), mom_shifti);
  Qassert(mom_shifti < int(mom_shiftL.size()));
  const Coordinate mom_shiftp = mom_shiftL[mom_shifti];

  for(int mg = 0; mg < mass_group; mg ++)
  {
    momentum_dat mdat(geo, in.mom_cut, mom_shiftL);
    char name[450],namecr[500],namesm[500],nameQ[550],name_mom[550], namei[500];

    std::string outputG      = ssprintf("%s%s", in.output.c_str(), flag_add.c_str());
    std::string output_vecG  = ssprintf("%s%s", in.output_vec.c_str(), flag_add.c_str());
    std::string output_vecGs = ssprintf("%s%s",   output_vecs.c_str(), flag_add.c_str());

    if(mass_group > 1){
      if(outputG != std::string("NONE")){
        outputG = ssprintf("%s.mgroup%02d", outputG.c_str(), mg);
      }

      if(output_vecG != std::string("NONE")){
        output_vecG  = ssprintf("%s.mgroup%02d", output_vecG.c_str(), mg);
        output_vecGs = ssprintf("%s.mgroup%02d", output_vecGs.c_str(), mg);
      }

    }

    std::vector<std::string > pt_infoS = stringtolist(in.paraIA);
    std::vector<std::string > gr_infoS = stringtolist(in.paraIB);
    std::vector<int > pt_info, gr_info;
    for(unsigned int i=0;i<pt_infoS.size();i++){pt_info.push_back(stringtonum(pt_infoS[i]));}
    for(unsigned int i=0;i<gr_infoS.size();i++){gr_info.push_back(stringtonum(gr_infoS[i]));}

    qassert(pt_info.size() == 2);qassert(gr_info.size() == 3);
    int npoint = pt_info[1] - pt_info[0];
    int ngrid  = gr_info[1] - gr_info[0];
    int nsource = npoint + ngrid ;

    ////==check setup and scales
    qassert(gr_info[1] > gr_info[0]);
    qassert(gr_info[2] > 1);
    ////===setup basic info and factors
    const int oper = 32;
    int nmass = 0;
    int nt    = 32;
    Coordinate pL;for(int i=0;i<3;i++){pL[i] = 9;}pL[3] = nt;
    double factor_pt = 0.0;int gr_start = 1;
    if(pt_info[1] == pt_info[0]){factor_pt = 0.0;gr_start = 0;}
    else{
      factor_pt = 1.0/( (pt_info[1]-pt_info[0]));
      gr_start  = 1;
    }
    double factor_gr = 1.0/( (gr_info[1]-gr_info[0])*(gr_info[2]-1) );
    ////===setup basic info and factors

    std::vector<Coordinate > momT;momT.resize(in.nt);
    for(int t=0;t<in.nt;t++){
      momT[t]    = mom_shiftp;
      momT[t][3] = t;
    }
    std::vector<bool > Bsite;
    std::vector<long > Zsite;
    mdat.get_mom_pos(momT, Bsite, Zsite);

    fft_desc_basic fd(geo);
    //for(unsigned int bi=0;bi<Bsite.size();bi++)
    //{
    //  if(Bsite[bi]){printf("node %d %3d \n", fd.rank, bi);}
    //}

    ////==check setup and scales

    //  sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", nsource, 2, 3, 32, int(massL.size()), in.nt, 2);
    //  sprintf(dimN , "src  sm  LHF operator masses nt complex");

    corr_dat<Ftype > average(std::string(""));
    corr_dat<Ftype > corr_zero(std::string(""));

    std::vector<bool > data_sm(2);
    ////===check data
    for(int sm = 0; sm < 2;sm ++)
    {
      data_sm[sm] = true;
      for(int si = gr_info[0]; si < gr_info[1]; si++)
      {
        sprintf(name, output_vecG.c_str(), icfg, si);
        if(sm == 0){sprintf(namesm, "%s.pt", name);}
        if(sm == 1){sprintf(namesm, "%s.sm", name);}
        sprintf(nameQ   , "%s.Gsrc" , namesm );
        sprintf(name_mom, "%s.GInfo", namesm );
        if(get_file_size_MPI(name_mom) == 0)
        {
          qmessage("File %s not found! \n",name_mom);
          data_sm[sm] = false;
        }

        if(average.dim == 0){
          corr_dat<Ftype > mom_info(name_mom);
          std::vector<Coordinate > src_pos = read_positions(mom_info);
          //std::vector<Coordinate > src_pos = string_to_Coordinates(mom_info.INFO_LIST);

          std::string ktem = std::string("32 ")       + mom_info.get_key_T();
          std::string dtem = std::string("operator ") + mom_info.get_dim_name();
          average.create_dat(ktem, dtem);

          qassert(mom_info.dim == 6);
          qassert(average.key_T[0] == oper);
          ////nmass = average.key_T[1];
          nt = average.key_T[2];
          for(int i=0;i<3;i++){pL[i] = average.key_T[i+3];}pL[3] = nt;

          char tem[500];sprintf(tem, "point %s, grid %s ", in.paraIA.c_str(), in.paraIB.c_str());
          average.INFO_LIST = std::string("Average data, ") + std::string(tem);
          average.INFOA = mom_info.INFOA;
          average.INFOA.push_back(std::string("mom shift ") + Coordinate_to_string(mom_shiftp));
        }

        ////omit high mode for grid source
        if(data_sm[sm])
        for(int gi = 0; gi < gr_info[2]; gi++)
        {
          fflush_MPI();sprintf(namei, "%09d", gi);
          if(!mdat.check_fn_momcut(std::string(nameQ), std::string(namei))){data_sm[sm] = false;}
        }
      }

      for(int si = pt_info[0]; si < pt_info[1]; si++)
      {
        fflush_MPI();
        sprintf(name, output_vecG.c_str(), icfg, si);

        if(sm == 0){sprintf(namesm, "%s.pt", name);}
        if(sm == 1){sprintf(namesm, "%s.sm", name);}
        sprintf(nameQ, "%s.Gsrc", namesm );
        sprintf(name_mom, "%s.GInfo", namesm );
        if(get_file_size_MPI(name_mom) == 0)
        {
          qmessage("File %s not found! \n",name_mom);
          data_sm[sm] = false;
        }

        if(data_sm[sm])
        for(int gi=0;gi<2;gi++){
          fflush_MPI();sprintf(namei, "%09d", gi);
          if(!mdat.check_fn_momcut(std::string(nameQ), std::string(namei))){data_sm[sm] = false;}
        }
      }
    }

    /////===zero momentum corr
    int Nsm = 2; 
    {
      char key_T[1000], dimN[1000];
      std::vector<double> massL =  string_to_mass(average.INFOA[0]);
      nmass = massL.size();
      if(!data_sm[1]){Nsm = 1;}
      sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", nsource, Nsm, 3, 32, int(massL.size()), in.nt, 2);
      sprintf(dimN , "src  sm  LHF operator masses nt complex");
      std::string ktem = std::string(key_T);
      std::string dtem = std::string(dimN);
      corr_zero.create_dat(ktem, dtem);
      corr_zero.INFOA = average.INFOA;
      corr_zero.set_zero();
    }
    ////===

    qassert(output_vecG!= std::string("NONE"));

    qlat::vector_gpu<Complexq > FFT_data;
    qlat::vector_gpu<Complexq > FFT_data_average;
    ////qlat::vector_gpu<Complexq > FFT_data_average;
    qlat::vector_gpu<Complexq > FFT_data_global;
    //qlat::vector_gpu<Complexq > FFT_data_global_cpu;

    print_mem_info();
    for(int sm = 0; sm < 2;sm ++)
    {
      if(data_sm[sm] == false){qmessage("Continue sm %2d \n", sm);continue ;}

      FFT_data_average.set_zero();

      ///std::vector<Complexq > phases;

      //qassert(gr_info[1] > gr_info[0]);
      //qassert(gr_info[2] > 1);

      int nsi = 0;
      for(int si = gr_info[0]; si < gr_info[1]; si++)
      {
        print_mem_info();
        sprintf(name, output_vecG.c_str(), icfg, si);
        if(sm == 0){sprintf(namesm, "%s.pt", name);}
        if(sm == 1){sprintf(namesm, "%s.sm", name);}
        sprintf(nameQ, "%s.Gsrc", namesm );
        sprintf(name_mom, "%s.GInfo", namesm );
        print_time();
        qmessage("  %s \n", nameQ);

        corr_dat<Ftype > mom_info(name_mom);
        std::vector<Coordinate > src_pos = read_positions(mom_info);
        //std::vector<Coordinate > src_pos = string_to_Coordinates(mom_info.INFO_LIST);
        mom_info.print_info();

        for(int gi = 0; gi < gr_info[2]; gi++)
        {
          fflush_MPI();sprintf(namei, "%09d", gi);
          if(mdat.read_momcut(FFT_data, std::string(nameQ), std::string(namei)) == 0){qassert(false);}
          if(FFT_data_average.size() ==  0){FFT_data_average.resize(FFT_data.size());FFT_data_average.set_zero();}

          mdat.apply_src_phases(FFT_data, src_pos[gi] );
          mdat.shift_t(FFT_data, FFT_data, src_pos[gi][3]);

          if(gi != 0){
          long  corr_zero_off = (((long(nsi) + npoint)*Nsm + sm)*3 +  2) * 32 * nmass * in.nt * 2;
                write_corr_zero(corr_zero, FFT_data, corr_zero_off, Bsite, Zsite);}
          if(gi == 0){
          long  corr_zero_off = (((long(nsi) + npoint)*Nsm + sm)*3 +  1) * 32 * nmass * in.nt * 2;
                write_corr_zero(corr_zero, FFT_data, corr_zero_off, Bsite, Zsite);}

          ////omit high mode for grid source
          if(gi >= gr_start){
            cpy_data_thread(FFT_data_average.data(), FFT_data.data(), FFT_data.size(), 1, QTRUE, factor_gr);
          }
        }
        nsi += 1;
      }

      nsi = 0;
      for(int si = pt_info[0]; si < pt_info[1]; si++)
      {
        print_mem_info();
        fflush_MPI();
        sprintf(name, output_vecG.c_str(), icfg, si);

        if(sm == 0){sprintf(namesm, "%s.pt", name);}
        if(sm == 1){sprintf(namesm, "%s.sm", name);}
        sprintf(nameQ, "%s.Gsrc", namesm );
        sprintf(name_mom, "%s.GInfo", namesm );

        corr_dat<Ftype > mom_info(name_mom);
        std::vector<Coordinate > src_pos = read_positions(mom_info);
        //std::vector<Coordinate > src_pos = string_to_Coordinates(mom_info.INFO_LIST);
        mom_info.print_info();

        qassert(mom_info.dim == 6);

        for(int gi=0;gi<2;gi++){
          fflush_MPI();sprintf(namei, "%09d", gi);
          if(mdat.read_momcut(FFT_data, std::string(nameQ), std::string(namei)) == 0){qassert(false);}

          if(FFT_data_average.size() ==  0){FFT_data_average.resize(FFT_data.size());FFT_data_average.set_zero();}
          mdat.apply_src_phases(FFT_data, src_pos[gi] );
          mdat.shift_t(FFT_data, FFT_data, src_pos[gi][3]);

          if(gi != 0){
          long  corr_zero_off = (((long(nsi) +   0)*Nsm + sm)*3 +  2) * 32 * nmass * in.nt * 2;
                write_corr_zero(corr_zero, FFT_data, corr_zero_off, Bsite, Zsite);}
          if(gi == 0){
          long  corr_zero_off = (((long(nsi) +   0)*Nsm + sm)*3 +  1) * 32 * nmass * in.nt * 2;
                write_corr_zero(corr_zero, FFT_data, corr_zero_off, Bsite, Zsite);}

          ////FFT_data[];
          if(gi == 0){
            cpy_data_thread(FFT_data_average.data(), FFT_data.data(), FFT_data.size(), 1, QTRUE, factor_pt);
          }
        }
        nsi += 1;
      }

      {
        sprintf(name, output_vecGs.c_str(), icfg, 0);
        if(sm == 0)sprintf(namecr, "%s.pt.average", name );
        if(sm == 1)sprintf(namecr, "%s.sm.average", name );

        fft_local_to_global(FFT_data_global, FFT_data_average, mdat, mom_shiftp);

        //average.set_write_lines(namecr);
        average.set_zero();
        average.write_corr(FFT_data_global.data(), FFT_data_global.size(), 3);

        average.print_info();
        average.write_dat(namecr);
      }
    }

    sum_all_size(&corr_zero.dat[0], corr_zero.dat.size(), 0);
    sprintf(namei, outputG.c_str(),icfg);
    corr_zero.print_info();
    corr_zero.write_dat(namei);
  }

  }


  return end_Lat();
}

