#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "utils_smear_vecs.h"

int main(int argc, char* argv[])
{
  using namespace qlat;
  inputpara in0;
  begin_Lat(&argc, &argv, in0);

  int nx,ny,nz,nt;
  nx = in0.nx;ny = in0.ny;nz = in0.nz;nt = in0.nt;

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 

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

  std::vector<qlat::FieldM<qlat::Complex , 12> > eigenD;
  std::vector<qlat::FieldM<qlat::ComplexF, 12> > eigenF;
  eigenD.resize(each);
  eigenF.resize(each);

  //std::vector<qlat::FieldM<qlat::Complex , 12> > noi_bufD;
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
    /////GaugeFieldT<qlat::Complex  > gfD;gfD = gfa;
    //qlat::Complex*  tem0 = (qlat::Complex*  ) qlat::get_data(gfD).data();
    //qlat::ComplexF* tem1 = (qlat::ComplexF* ) qlat::get_data(gfF).data();
    //LInt cpy_size = (LInt)(qlat::get_data(gfD).data_size()/(sizeof(qlat::Complex)));
    //cpy_data_thread(tem1, tem0, cpy_size, 0);
  }

  FILE* file_read  = open_eigensystem_file(ename, nini, nvec, true , io_read , in_read_eigen , 2);
  FILE* file_write = open_eigensystem_file(Sname, nini, nvec, false, io_write, in_write_eigen, stringtonum(in0.save_type));

  std::vector<long > job =  job_create(nread, each);
  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  {
    /////int n0 = nini + job[ji*2 + 0]; int n1 = n0 + job[ji*2 + 1];
    if(!F_single){load_eigensystem_vecs(file_read ,   eigenD, io_read , in_read_eigen , 0, job[ji*2 + 1]);}
    if( F_single){load_eigensystem_vecs(file_read ,   eigenF, io_read , in_read_eigen , 0, job[ji*2 + 1]);}

    ////for(long iv=0;iv < job[ji*2 + 1];iv++){noi_bufD[iv] = eigenD[iv];}
    ////for(long iv=0;iv < job[ji*2 + 1];iv++){noi_bufF[iv] = eigenF[iv];}

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


  //////===Open file for write
  //inputpara in_writeqlat;in_writeqlat.read_geo(geo);
  //int bfac_eigen = 12;
  //bool read = false;
  //std::string VECS_TYPE = std::string("Eigen_system_nvec.12.tzyx.R/I");
  //std::string INFO_LIST = std::string("NONE");
  //bool rotate_bfac = true;
  //open_file_qlat_noisesT(Sname, bfac_eigen, in_writeqlat, read, save_single, nread, VECS_TYPE, INFO_LIST, rotate_bfac);
  //////===Open file for write

  
  //{
  //std::vector<qlat::FermionField4dT<Complex > > eigen;
  //load_gwu_eigen(ename,eigen,io_use, n0 , n1,true, true);
  //double* tmp = (double*)(qlat::get_data(eigen[1]).data());
  //print0("value pre %.3e %.3e \n", tmp[0], tmp[1]);
  //}

  //if(file_type == 2 or file_type == 3){
  //  if(file_type == 2){load_qlat_eigen(ename, eigenD, n0, n1);}
  //  if(file_type == 3){load_qlat_eigen(ename, eigenF, n0, n1);}
  //}else{
  //  bool read = true;bool check = true;

  //  std::vector<double* > respD;respD.resize(nread);
  //  std::vector<float*  > respF;respF.resize(nread);
  //  eigenD.resize(nread);
  //  eigenF.resize(nread);
  //  for(int iv=0;iv<nread;iv++){
  //    if(file_type == 0){eigenD[iv].init(io_use.geop);respD[iv]=(double*)(qlat::get_data(eigenD[iv]).data());}
  //    if(file_type == 1){eigenF[iv].init(io_use.geop);respF[iv]=(float* )(qlat::get_data(eigenF[iv]).data());}
  //  }

  //  if(file_type == 0){load_gwu_eigen(ename, respD, io_use, n0, n1, check, read);}
  //  if(file_type == 1){load_gwu_eigen(ename, respF, io_use, n0, n1, check, read);}
  //}

  //bool single_file = 0;
  //if(file_type == 0){single_file = false;}
  //if(file_type == 1){single_file = true ;}
  //if(file_type == 2){single_file = false;}
  //if(file_type == 3){single_file = true ;}


  //double* tmp = (double*)(qlat::get_data(eigenD[1]).data());
  //print0("values %.3e %.3e \n", tmp[1], tmp[2]);

  //if(!single_file){
  //  io_vec io_use(geo, IO_DEFAULT);
  //  bool read = false;
  //  int save_type = 1;
  //  ////if(save_single == true){save_type = 3;}

  //  int each = io_use.ionum;
  //  std::vector<qlat::FieldM<qlat::Complex , 12> > noi_buf;noi_buf.resize(each);
  //  for(int iv=0;iv<each;iv++){noi_buf[iv].init(geo);}

  //  inputpara in;
  //  FILE* file = open_eigensystem_file(Sname, eigenD.size(), read, io_use, in, save_type);

  //  std::vector<long > job =  job_create(eigenD.size(), each);
  //  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  //  {
  //    for(long iv=0;iv < job[ji*2 + 1];iv++){noi_buf[iv] = eigenD[job[ji*2 + 0] + iv];}
  //    load_eigensystem_vecs(file, noi_buf, io_use, in, job[ji*2 + 0], job[ji*2 + 0] + job[ji*2 + 1]);
  //    /////load_qlat_noisesT(file, noi_buf, io_use, in, job[ji*2 + 0], job[ji*2 + 0] + job[ji*2 + 1]);
  //  }

  //  close_eigensystem_file(file, io_use, in);

  //}
  //if( single_file){save_qlat_eigen(Sname, eigenF, save_single, in_eigen.INFO_LIST);}


  //if(!single_file){
  //  std::vector<qlat::FermionField4dT<qlat::Complex  > > eigen;eigen.resize(eigenD.size());
  //  LInt cpy_size = (LInt)(qlat::get_data(eigenD[0]).data_size()/(sizeof(qlat::Complex)));
  //  for(LInt iv=0;iv<eigenD.size();iv++){
  //    eigen[iv].init(geo);
  //    qlat::Complex* p0 = (qlat::Complex*)(qlat::get_data(eigenD[iv]).data());
  //    qlat::Complex* p1 = (qlat::Complex*)(qlat::get_data(eigen[iv]).data());
  //    cpy_data_thread(p1, p0, cpy_size, 1);
  //  }
  //  io_vec io_use(geo, IO_DEFAULT);
  //  save_gwu_eigen(Sname, eigen, io_use, 0, eigen.size(), true);
  //  /////save_qlat_eigen(Sname, eigenD, true, in_eigen.INFO_LIST);
  //}


  //if(!single_file){
  //  io_vec io_use(geo, IO_DEFAULT);
  //  bool read = true;
  //  int save_type = 2;
  //  //if(save_single == true){save_type = 3;}

  //  int each = io_use.ionum;
  //  std::vector<qlat::FieldM<qlat::Complex , 12> > noi_buf;noi_buf.resize(each);
  //  for(int iv=0;iv<each;iv++){noi_buf[iv].init(geo);}

  //  inputpara in;
  //  FILE* file = open_eigensystem_file(Sname, eigenD.size(), read, io_use, in, save_type);

  //  std::vector<long > job =  job_create(eigenD.size(), each);
  //  for(LInt ji = 0; ji < job.size()/2 ; ji++)
  //  {
  //    load_eigensystem_vecs(file, noi_buf, io_use, in, job[ji*2 + 0], job[ji*2 + 0] + job[ji*2 + 1]);
  //    for(long iv=0;iv < job[ji*2 + 1];iv++){eigenD[job[ji*2 + 0] + iv] = noi_buf[iv];}
  //  }

  //  close_eigensystem_file(file, io_use, in);

  //}
  //if( single_file){load_qlat_eigen(Sname, eigenF);}

  //if(!single_file){save_qlat_eigen(Sname, eigenD, save_single, in_eigen.INFO_LIST);}
  //if( single_file){save_qlat_eigen(Sname, eigenF, save_single, in_eigen.INFO_LIST);}


  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}

