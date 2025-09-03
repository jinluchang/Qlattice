#include <sys/sysinfo.h>
#include <unistd.h>

#include "general_funs.h"
#include "utils_io_vec.h"
#include "utils_construction.h"
#include "utils_stagger_contractions.h"
#include "utils_gaugefield.h"
#include "utils_shift_vecs.h"
#include "utils_lms_funs.h"

namespace qlat{

inline void write_corr_zero(corr_dat<Ftype >& res, qlat::vector_gpu<Complexq >& FFT_data,
    std::vector<bool >& Bsite, std::vector<long >& Zsite, const int nt, const int nread)
{
  ////qmessage("write data %d \n", int(off));
  //const int nmass = res.key_T[4];
  //const int nt    = res.key_T[5];
  //Complexq* Pres = (Complexq*) (&res.dat[off]);
  Qassert(int(Bsite.size()) == nt);
  Qassert(int(Zsite.size()) == nt);
  //const int nvec = 1;
  //const int nvec = 32*nmass;
  qlat::vector<Complexq > buf_gpu;buf_gpu.resize(nread * nt);
  qlat::set_zero(buf_gpu);
  Complexq* Pbuf = (Complexq*) (qlat::get_data(buf_gpu).data());
  Complexq* PFFT = (Complexq*) (qlat::get_data(FFT_data).data());
  qlat::vector<long> PZsite;PZsite.resize(Zsite.size());
  for(Long t=0;t<PZsite.size();t++){PZsite[t] = Zsite[t];}

  ////sprintf(key_T, "%d   %d    %d   %d  %d  %d %d", nsource, 2, 3, 32, int(massL.size()), in.nt, 2);
  ////long corr_zero_off = (((si + npoint)*2 + sm)*3 + LHF) * 32 * nmass * in.nt * 2;
  for(int t=0;t<nt;t++){if(Bsite[t]){
  qacc_for(iv, nread,{
    Pbuf[iv*nt + t] += PFFT[ PZsite[t]*nread + iv ] ;
  });
  }}

  //for(int t=0;t<nt;t++){if(Bsite[t]){
  //qthread_for(iv, nread,{
  //    Pres[iv*nt + t] += Pbuf[iv*nt + t];
  //});}}
  res.write_corr(buf_gpu.data(), buf_gpu.size());

}
}

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in;
  begin_Lat(&argc, &argv, in);

  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;

  {
  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 
  fflush_MPI();

  std::string flag_add = std::string("FFT");
  in.find_para(std::string("flag_add"), flag_add);

  int do_corr = 1;
  int do_fft  = 1;
  int corr_with_info = 1;
  in.find_para(std::string("do_corr"), do_corr);
  in.find_para(std::string("do_fft" ), do_fft );
  in.find_para(std::string("corr_with_info" ), corr_with_info );

  int extra_info_flag = 0;
  in.find_para(std::string("extra_info_flag"), extra_info_flag);

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

  momentum_dat mdat(geo, in.mom_cut, mom_shiftL);

  // setup zero corr writting with momentum shift if needed
  std::vector<Coordinate > momT;momT.resize(in.nt);
  for(int t=0;t<in.nt;t++){
    momT[t]    = mom_shiftp;
    momT[t][3] = t;
  }   
  std::vector<bool > Bsite;
  std::vector<long > Zsite;
  mdat.get_mom_pos(momT, Bsite, Zsite);

  fft_desc_basic fd(geo);
  for(unsigned int bi=0;bi<Bsite.size();bi++)
  {   
    if(Bsite[bi]){printf("node %d %3d \n", fd.rank, bi);}
  } 
    
  qlat::vector_gpu<Complexq > FFT_data;
  qlat::vector_gpu<Complexq > FFT_data_global;
  qlat::vector<Complexq > tmp_corr;

  std::vector<std::string > info_tmp;
  
  if(in.paraIA != "NONE"){
    info_tmp = stringtolist(in.paraIA);
    Qassert(long(info_tmp.size()) == in.nsource);
  }

  char key_T[1000], dimN[1000];

  const int    mc = mdat.mom_cut*2 + 1;
  const long mvol = long(mdat.nv[3])*mc*mc*mc;;
  for(int isource = 0; isource < in.nsource; isource++)
  if(in.propN[isource] != std::string("NONE"))
  {
    int with_info = 0;
    corr_dat<Ftype > average(std::string(""));
    corr_dat<Ftype > corr_zero(std::string(""));
    corr_dat<Ftype > corr_info(std::string(""));
    //  long  corr_zero_off = 0;

    std::string nameQ = in.propN[isource];

    std::string name_fft ;
    std::string name_corr;
    std::string name_info;

    // check file done
    {
      int fft_done  = 0;
      int corr_done = 0;
      name_fft = qlat::ssprintf("%s.%s", nameQ.c_str(), flag_add.c_str());
      // need to check writeout sizes 
      if(get_file_size_MPI(name_fft.c_str()) >= size_t(32 * mvol * 2 * 8) and check_qlat_head_exist(name_fft.c_str())){
        fft_done = 1;
      }

      name_corr = qlat::ssprintf("%s.zero.corr", nameQ.c_str());
      if(get_file_size_MPI(name_corr.c_str()) >= size_t(32 * long(mdat.nv[3]) * 2 * 8)){corr_done = 1;}

      bool DONE = true;
      if(do_fft  == 1 and fft_done  != 1){DONE = false;}
      if(do_corr == 1 and corr_done != 1){DONE = false;}

      if(DONE == true)
      {
        std::string namew = qlat::ssprintf("%s", nameQ.c_str());
        qmessage("Pass %s \n", namew.c_str());
        continue ;
      }
    }

    name_info = qlat::ssprintf("%s", nameQ.c_str());
    name_info = name_info + std::string(" ");
    name_info.replace(name_info.end()-5 - extra_info_flag, name_info.end(), qlat::ssprintf("GInfo"));
    //qmessage("%s \n", name_info.c_str());
    if(get_file_size_MPI(name_info.c_str()) > 0){
      with_info = 1;
      corr_info.read_dat(name_info.c_str());
      Qassert(corr_info.dim == 6);
    }

    // continue if info not found
    if(corr_with_info == 1 and with_info == 0)
    {
      qmessage("Pass %s \n", in.propN[isource].c_str());
      continue;
    }

    mdat.get_fn_list(nameQ);
    //// int nvec = -1;
    qmessage("isource %d, current %d \n", isource, int(mdat.fn_list.size()));

    if(info_tmp.size() != 0){
      int check = stringtonum(info_tmp[isource]);
      if(long(mdat.fn_list.size()) != long(check)){
        abort_r("Wrong fn list");
        qmessage("current %d , expect %d \n", int(mdat.fn_list.size()), check);
      }
    }

    for(unsigned int fi=0;fi<mdat.fn_list.size();fi++){
      if(do_fft == 1){
        if(average.dim == 0)
        {
          const int nvec = mdat.read(FFT_data, nameQ, mdat.fn_list[fi]);
          check_nan_GPU(FFT_data);
          //fft_local_to_global(FFT_data_global, FFT_data, mdat, mom_shiftp);
          //nvec = FFT_data_global.size()/mvol;
          ///sprintf(key_T, "%d  %d %d  %d %d %d %d", int(mdat.fn_list.size()), nvec, in.nt, mc, mc, mc, 2);
          sprintf(key_T, "%d  %d %d  %d %d %d %d", int(1), nvec, in.nt, mc, mc, mc, 2);
          sprintf(dimN , "sources ch nt pz py px complex");
          average.INFO_LIST = in.paraI;
          average.create_dat(std::string(key_T), std::string(dimN));
          average.INFOA.resize(0);
          if(with_info == 1){
            average.INFO_LIST = corr_info.INFO_LIST;
            average.INFOA = corr_info.INFOA;
          }
        }
        average.INFOA.push_back(mdat.fn_list[fi]);
        average.INFOA.push_back(std::string("mom shift ") + Coordinate_to_string(mom_shiftp));
      }

      if(do_corr == 1){
        if(corr_zero.dim == 0)
        {
          const int nvec = mdat.read(FFT_data, nameQ, mdat.fn_list[fi]);
          check_nan_GPU(FFT_data);
          //fft_local_to_global(FFT_data_global, FFT_data, mdat, mom_shiftp);
          //nvec = FFT_data_global.size()/mvol;
          ///sprintf(key_T, "%d  %d %d  %d %d %d %d", int(mdat.fn_list.size()), nvec, in.nt, mc, mc, mc, 2);
          sprintf(key_T, "%d  %d %d %d", int(1), nvec, in.nt, 2);
          sprintf(dimN , "sources ch nt complex");
          corr_zero.INFO_LIST = in.paraI;
          corr_zero.create_dat(std::string(key_T), std::string(dimN));
          corr_zero.INFOA.resize(0);
          if(with_info == 1){
            corr_zero.INFO_LIST = corr_info.INFO_LIST;
            corr_zero.INFOA = corr_info.INFOA;
          }
        }
        corr_zero.INFOA.push_back(mdat.fn_list[fi]);
        corr_zero.INFOA.push_back(std::string("mom shift ") + Coordinate_to_string(mom_shiftp));
      }
    }


    //char namew[500];sprintf(namew, "%s.%s", nameQ.c_str(), flag_add.c_str());
    if(do_fft == 1){
      average.set_write_lines(name_fft.c_str());
    }

    for(unsigned int fi=0;fi<mdat.fn_list.size();fi++){
      // should be propotional ch = 32 * nmass
      const int nread = mdat.read(FFT_data, nameQ, mdat.fn_list[fi]);
      //Qassert(nread % 32 == 0);
      //const int nmass = nread / 32;
      check_nan_GPU(FFT_data);
      if(do_fft == 1){
        fft_local_to_global(FFT_data_global, FFT_data, mdat, mom_shiftp);
        Qassert(long(FFT_data_global.size()) == long(nread) * mvol);
        average.write_corr(FFT_data_global.data(), FFT_data_global.size(), 3);
      }
      if(do_corr == 1){
        // initialize memeories
        //tmp_corr.resize(nread * in.nt);
        //corr_zero.write_corr(tmp_corr.data(), tmp_corr.size());
        //for(int ni = 0; ni< nread; ni++)
        {
          // writout offset 
          //const long  corr_zero_off = ni * in.nt * 2;
          write_corr_zero(corr_zero, FFT_data, Bsite, Zsite, in.nt, nread);
          //corr_zero_off += nread * in.nt * 2;
        }
      }
      if(fi%20 == 0){qmessage("  Tag %s END\n", mdat.fn_list[fi].c_str());}
    }

    if(do_fft == 1){
      average.print_info();
      average.write_dat(name_fft.c_str());
    }
    if(do_corr == 1){
      // t maybe on different nodes
      sum_all_size(&corr_zero.dat[0], corr_zero.dat.size(), 0);
      corr_zero.print_info();
      corr_zero.write_dat(name_corr.c_str());
    }

    fflush_MPI();

  }

  }

  //fflush_MPI();
  //qlat::Timer::display();
  //qlat::end();
  //return 0;
  return end_Lat();
}

