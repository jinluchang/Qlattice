#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "utils_smear_vecs.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in0;
  begin_Lat(&argc, &argv, in0);

  {
  int nx,ny,nz,nt;
  nx = in0.nx;ny = in0.ny;nz = in0.nz;nt = in0.nt;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site); 

  int icfg  = in0.icfg;
  /////int ionum = in0.ionum;

  int nini  = in0.nini;
  int nvec  = in0.nvec;
  if(nvec <= nini){abort_r("No vectors for reading!!! \n");}
  int nread = nvec - nini;
  ////int n0 = nini; int n1 = nvec;

  char ename[450],enamev[500];
  char Sname[450],Snamev[500], sDescription[500];
  ////char Ename[500];

  sprintf(ename ,in0.Ename.c_str(), icfg);
  sprintf(enamev,"%s.eigvals"    ,ename);

  sprintf(Sname ,in0.Sname.c_str(), icfg);
  sprintf(Snamev,"%s.eigvals"    ,Sname);
  std::vector<double > values, errors;
  std::vector<double > v0, e0;

  if(in0.debuga == 0)sprintf(sDescription, "Overlap operator: half spectrum");
  if(in0.debuga == 1)sprintf(sDescription, "H_wilson");
  load_gwu_eigenvalues(values,errors, enamev);
  v0.resize(nread*2); e0.resize(nread);
  for(int iv=0;iv<nread*2;iv++){v0[iv] = values[nini*2 + iv];}
  for(int iv=0;iv<nread;iv++){  e0[iv] = errors[nini + iv];}
    
  save_gwu_eigenvalues(v0,e0, Snamev, sDescription);

  ///bool save_single = true;
  ///save_single = get_save_type(in0.save_type);

  print_mem_info();

  double length = (geo.local_volume()*pow(0.5,30))*12*sizeof(Complexq);
  print0("Eign system vector size %.3e GB, total %.3e GB; \n", length, nvec*length);

  int Nv = omp_get_max_threads();
  io_vec io_read( geo, IO_DEFAULT, true, Nv);
  io_vec io_write(geo, IO_DEFAULT, true, Nv);

  int each = io_read.ionum;
  inputpara in_tem;
  int file_type = check_qlat_eigen_file_type(ename, io_read, nvec, in_tem);

  std::vector<qlat::FieldM<qlat::ComplexD , 12> > eigenD;
  std::vector<qlat::FieldM<qlat::ComplexF, 12> > eigenF;
  eigenD.resize(each);
  eigenF.resize(each);

  //std::vector<qlat::FieldM<qlat::ComplexD , 12> > noi_bufD;
  //std::vector<qlat::FieldM<qlat::ComplexF, 12> > noi_bufF;
  //noi_bufD.resize(each);
  //noi_bufF.resize(each);

  bool F_single = false;
  if(file_type == 0 or file_type == 2){F_single = false;}
  if(file_type == 1 or file_type == 3){F_single = true ;}

  for(int iv=0;iv<each;iv++){
    if(!F_single){eigenD[iv].init(geo);}
    if( F_single){eigenF[iv].init(geo);}

    //if(!F_single){eigenD[iv].init(geo);noi_bufD[iv].init(geo);}
    //if( F_single){eigenF[iv].init(geo);noi_bufF[iv].init(geo);}
  }

  inputpara in_read_eigen;
  inputpara in_write_eigen;

  ////bool read = false;
  ////int save_type = 2;
  ////if(save_single == true){save_type = 3;}

  /////===smearing of eigen
  double width = 0.0; int step = 0;
  get_smear_para(in0.src_smear_para, width, step);
  GaugeField gf;
  //GaugeField gfD;

  if(step != 0){
    gf.init(geo);
    char rbc_conf[500];
    sprintf(rbc_conf,in0.Link_name.c_str(), icfg);
    load_gwu_link(rbc_conf, gf);
    ////random_link(gf, 0); 
    //set_left_expanded_gauge_field(gfD, gf);
    /////GaugeFieldT<qlat::ComplexD  > gfD;gfD = gfa;
    //qlat::ComplexD*  tem0 = (qlat::ComplexD*  ) qlat::get_data(gfD).data();
    //qlat::ComplexF* tem1 = (qlat::ComplexF* ) qlat::get_data(gfF).data();
    //LInt cpy_size = (LInt)(qlat::get_data(gfD).data_size()/(sizeof(qlat::ComplexD)));
    //cpy_data_thread(tem1, tem0, cpy_size, 0);
  }

  FILE* file_read  = open_eigensystem_file(ename, nini, nvec, true , io_read , in_read_eigen , 2);
  FILE* file_write = open_eigensystem_file(Sname, nini, nvec, false, io_write, in_write_eigen, stringtonum(in0.save_type));

  std::vector<Long > job =  job_create(nread, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    /////int n0 = nini + job[ji*2 + 0]; int n1 = n0 + job[ji*2 + 1];
    if(!F_single){load_eigensystem_vecs(file_read ,   eigenD, io_read , in_read_eigen , 0, job[ji*2 + 1]);}
    if( F_single){load_eigensystem_vecs(file_read ,   eigenF, io_read , in_read_eigen , 0, job[ji*2 + 1]);}

    ////for(Long iv=0;iv < job[ji*2 + 1];iv++){noi_bufD[iv] = eigenD[iv];}
    ////for(Long iv=0;iv < job[ji*2 + 1];iv++){noi_bufF[iv] = eigenF[iv];}

    /////===smearing of eigen
    if(step != 0){
      for(int iv=0;iv<each;iv++){
      if(!F_single){smear_propagator_gwu_convension(eigenD[iv], gf, width, step);}
      if( F_single){smear_propagator_gwu_convension(eigenF[iv], gf, width, step);}}
    }
    /////===smearing of eigen

    if(!F_single){load_eigensystem_vecs(file_write, eigenD, io_write, in_write_eigen, 0, job[ji*2 + 1]);}
    if( F_single){load_eigensystem_vecs(file_write, eigenF, io_write, in_write_eigen, 0, job[ji*2 + 1]);}

  }

  close_eigensystem_file(file_read , io_read , in_read_eigen );
  close_eigensystem_file(file_write, io_write, in_write_eigen);

  }

  return end_Lat();

}

