// Gen Wang
// Sep. 2021


#ifndef UTILS_READ_TXT_H
#define UTILS_READ_TXT_H
#pragma once

////#include "general_funs.h"
#include "utils_copy_data.h"

namespace qlat
{

inline bool is_big_endian_gwu(void){
  int num = 1;
  if(*(char *)&num == 1)
  {
    return false;
    //printf("\nLittle-Endian\n");
  }
  else
  {
    return true;
    //printf("Big-Endian\n");
  }
}

inline void swapbytes(void *_object, size_t _size)
{
   unsigned char *start, *end;
   for ( start = (unsigned char *)_object, end = start + _size - 1; start < end; ++start, --end )
   {
      unsigned char swap = *start;
      *start = *end;
      *end = swap;
   }
}

inline void switchendian(void *buffer,size_t length,int dsize)
{
  #pragma omp parallel for
  for(size_t i=0;i<length;i++)
  {
    char* pos = ((char*)buffer) + i*dsize;
    swapbytes(pos,dsize);
  }

}

inline std::vector<std::string > stringtolist(const std::string &tem_string)
{
  std::istringstream iss(tem_string);
  std::vector<std::string> results((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
  return results;
}

inline std::string listtostring(const std::vector<int > src)
{
  char tmp[1000];
  std::string buf;
  for(unsigned int i=0;i<src.size();i++){
    sprintf(tmp,"%d ",src[i]);
    buf += std::string(tmp);
  }
  return buf;
}


inline std::string listtostring(const std::vector<std::string > src)
{
  char tmp[1000];
  std::string buf;
  for(unsigned int i=0;i<src.size();i++){
    sprintf(tmp,"%s ",src[i].c_str());
    buf += std::string(tmp);
  }
  return buf;
}


inline double stringtodouble(std::string &tem_string)
{
  //double use = atof(tem_string.c_str());
  double use = 0.0;
  if(tem_string!="_NONE_")
  {
    use = atof(tem_string.c_str());
  }
  return use;
}

inline int stringtonum(std::string &tem_string)
{
  //int t_Total = 0;
  //if(tem_string!="_NONE_")
  //{
  //  int tem_length = strlen(tem_string.c_str());
  //  for(int i=0;i<tem_length;i++){t_Total = t_Total+(tem_string.c_str()[i]-'0')*std::pow(10,tem_length-i-1);};
  //}
  //return t_Total;
  double use = 0.0;
  if(tem_string!="_NONE_")
  {
    use = atof(tem_string.c_str());
  }
  return int(use);

}

inline Coordinate string_to_Coordinate(const std::string& paraI = std::string("NONE"))
{
  Coordinate sp;for(int i=0;i<4;i++){sp[i] = 0;}
  if(paraI != "NONE"){
    std::vector<std::string > Li = stringtolist(paraI);
    qassert(Li.size() == 4);
    for(int i=0;i<4;i++){sp[i] = stringtonum(Li[i]);}
  }
  return sp;
}

inline std::string mass_to_string(std::vector<double>& massL)
{
  std::string mL;char mnum[500];
  mL += std::string("masses ");
  for(unsigned int mi=0;mi<massL.size();mi++)
  {
    sprintf(mnum, "   %.9f", massL[mi]);
    mL += std::string(mnum);
  }
  return mL;
}

inline std::vector<double> string_to_mass(const std::string& INFO_MASS)
{
  std::vector<std::string > resv = stringtolist(INFO_MASS);
  std::vector<double> massL;
  for(unsigned int i=1;i < resv.size(); i++)
  {
    massL.push_back(stringtodouble(resv[i]));
  }
  return massL;
}

inline unsigned long get_file_size_o(const char *filename)
{
  std::ifstream File(filename);
  if(!File.is_open()){if(qlat::get_id_node() == 0){printf("file is not exist\n");}return 0;}
  unsigned long Begin = File.tellg();
  File.seekg(0, std::ios_base::end);
  unsigned long End = File.tellg();
  File.close();
  return End-Begin;
}

inline size_t get_file_size_MPI(const char *filename, bool silence = false)
{
  size_t sizen = 0;
  if(qlat::get_id_node()==0){
    std::ifstream File(filename);
    if(!File.is_open()){
      if(!silence)if(qlat::get_id_node() == 0){printf("%s file is not exist\n",  filename);}
      sizen = 0;
    }
    else{
      unsigned long Begin = File.tellg();
      File.seekg(0, std::ios_base::end);
      unsigned long End = File.tellg();
      File.close();
      sizen = End-Begin;
    }
  }
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, get_comm());
  return sizen;
}


template<typename Ty>
void write_data(Ty* dat, FILE* file, size_t size, bool read=false, bool single_file = false){
  /////Check whether node zero is necessary
  TIMER("Single node write");
  if(qlat::get_id_node()==0){
    size_t sem = 0;
    qassert(sizeof(Ty) == sizeof(float) or sizeof(Ty) == sizeof(double));
    int bsize = sizeof(double);
    if(single_file == true){bsize = sizeof(float);}

    ///////data for analysis is small endian
    bool Rendian = true;

    //int size = 0;
    //if(read==false)size = dat.size();
    //if(read==true){size_t sizeF = get_file_size_o(filename);size = sizeF/bsize;dat.resize(size);}

    ////Set buf with size
    //char* buf=NULL;
    ////buf = new char[size*sizeof(double)];
    //buf = (char *)aligned_alloc_no_acc(size* bsize);
    std::vector<char > buf; buf.resize(size * bsize);

    ////Open file
    ////FILE* file = NULL;
    //if(read==false)file = fopen(filename, "wb");
    //if(read==true )file = fopen(filename, "rb");

    /////Switch endian of the file write
    if(read==false){
      if(single_file == false)cpy_data_thread((double*)(&buf[0]), &dat[0], size, 0);
      if(single_file == true )cpy_data_thread((float* )(&buf[0]), &dat[0], size, 0);
      /////memcpy(&buf[0],&dat[0],size*sizeof(double));
      if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
      if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
    }

    if(read==false){sem = fwrite(&buf[0], size*bsize, 1, file);}
    if(read==true ){sem =  fread(&buf[0], size*bsize, 1, file);}
    if(sem != 1){printf("Reading/Writing error %zu %zu \n", sem, size_t(1) );}

    /////Switch endian of the file write
    if(read==true ){
      if(Rendian == false)if(!is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
      if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
      ////memcpy(&dat[0],&buf[0],size*sizeof(double));
      if(single_file == false)cpy_data_thread(&dat[0], (double*)&buf[0], size, 0);
      if(single_file == true )cpy_data_thread(&dat[0], (float* )&buf[0], size, 0);
    }

    ////delete []buf;
    ////free(buf);
  }

}

template<typename Ty>
void write_data(std::vector<Ty > dat,const char *filename, bool read=false, bool single_file = false){
  if(qlat::get_id_node()==0)
  {
    int bsize = sizeof(double);
    if(single_file == true){bsize = sizeof(float);}

    int size = 0;
    if(read==false)size = dat.size();
    if(read==true){size_t sizeF = get_file_size_o(filename);size = sizeF/bsize;dat.resize(size);}
    FILE* file = NULL;
    if(read==false)file = fopen(filename, "wb");
    if(read==true )file = fopen(filename, "rb");

    write_data(&dat[0], file, size, read, single_file);
    fclose(file);file = NULL;
  }

}

template<typename Ty>
void read_data(std::vector<Ty > dat,const char *filename, bool single_file = false){
  write_data(dat, filename, true, single_file);
}



inline size_t read_input(const char *filename,std::vector<std::vector<std::string > > &read_f)
{
  read_f.resize(0);
  if(get_file_size_o(filename) == 0){print0("input file size zero %s !\n", filename);return 0;}
  FILE* filer = fopen(filename, "r");
  //////Can only be LINE_LIMIT length string
  char sTemp[LINE_LIMIT+1],tem[LINE_LIMIT+1];
  ///std::string s0(sTemp);
  if (filer == NULL){printf("Error opening file");return 0;}

  int count_line = 0;
  bool binary = 0;
  ////while(!feof(filer))
  /////maximum input line 5000
  for(int i=0;i<5000;i++)
  {
    tem[LINE_LIMIT] = 0;
    if(fgets(tem, LINE_LIMIT + 1, filer) == NULL){binary = true;break;};
    ///for(int j=0;j<1000;j++){if(tem[i] < 0){binary = true;break;}}if(binary){break};
    /////If the file is binary
    /////print0("==%s \n",tem);
    if(tem[LINE_LIMIT] != 0){binary = true;break;}
    if(tem[0]    <  0){binary = true;break;}////need to check whether it works or not
    if(std::string(tem).size() >= LINE_LIMIT){binary = true;break;}
    /////If the file is not binary
    ///printf("line %d %s", count_line, tem);
    if(std::string(tem).size() >= 2){
      sprintf(sTemp,"%s",tem);
      std::string s0(sTemp);
      if(s0 == std::string("END_OF_HEAD\n")){break;}
      std::vector<std::string > resv = stringtolist(s0);
      read_f.push_back(resv);
    }
    count_line += 1;
    if(feof(filer)){break;}
  }
  if(count_line == 5000){binary = true;}
  size_t off_file = ftell(filer);
  if(binary){printf("Binary file or file line too long! \n");read_f.resize(0);off_file = 0;}
  fclose(filer);

  return off_file;
}

////Bcast conf_l from zero rank
inline void bcast_vstring(std::vector<std::string> &conf_l, const int Host_rank = 0){

  int rank = get_node_rank_funs0();
  ////Bcast strings
  size_t sizen = 0;if(rank == Host_rank)sizen = conf_l.size();
  MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, Host_rank, get_comm());
  if(rank != Host_rank)conf_l.resize(sizen);
  for(unsigned int is=0;is<conf_l.size();is++){
    if(rank == Host_rank)sizen = conf_l[is].size();
    MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, Host_rank, get_comm());

    if(rank != Host_rank)conf_l[is].resize(sizen);
    MPI_Bcast(&conf_l[is][0], sizen, MPI_CHAR, Host_rank, get_comm());
  }
  ////Bcast strings

}

struct inputpara{
  std::vector<std::vector<std::string > > read_f;
  int bSize;
  int bSum;
  int cutN;
  std::string lat;
  int icfg;
  int icfg_end;
  int icfg_jump;
  int save_prop;

  int bini, ncut0, ncut1;

  int nx;
  int ny;
  int nz;
  int nt;
  int nini;
  int nvec;
  int mode_dis;
  std::string layout;
  int ndouble;
  int nsave;
  int bfac;
  int ionum;
  int seed;
  int hyp;
  int sparsefactor;
  int gridtem;
  int split_save;

  int lms;
  int combineT;
  int mom_cut;
  int mode_FFT_MPI;

  double Eerr;
  int SRC_PROP_WITH_LOW;

  std::string Link_name;
  std::string Ename;
  std::string Ename_Sm;
  std::string Sname;
  std::string Pname;

  std::string paraI;
  std::string paraIA;
  std::string paraIB;
  std::string paraIC;
  std::string paraI_2pt;
  std::string paraI_3pt;
  std::string sink_list;
  std::string para_stag;
  std::string src_smear_para;
  std::string sink_smear_para;
  int write_mode;

  //////qlat data tag
  std::string job_tag;

  int nprop;
  std::string Propname;
  std::string Srcname;
  std::string output;
  std::string output_vec;
  int do_all_low;
  int nmass;
  std::vector<double > masses;

  ////single inversion mass
  double fermion_mass;
  int niter;
  double cg_err;
  int solver_type;
  int inv_deflate;
  int fermion_type;
  int prec_type;

  int niter_sloppy;

  int eig_poly_deg;
  double eig_amin;
  double eig_err;

  int nsource;
  int read_noi;
  std::vector<std::string > propN;
  std::vector<std::string > srcN;
  std::vector<std::string > smearN;

  int debuga;
  bool printlog;
  int GPU;

  ////Double, Single .....
  std::string OBJECT;

  std::string save_type;
  std::string total_size;

  std::string key_T;
  std::string dim_name;
  std::string corr_name;
  std::string INFO_LIST;
  std::string FILE_ENDIAN;
  std::string VECS_TYPE;
  std::vector<std::string > INFOA;

  size_t off_file;
  crc32_t checksum;

  ////===clover paras
  double kappa;
  double clover_csw;
  ////===clover paras

  ////===private usage, not loaded from file head
  int    bsize;
  bool   read;
  bool   single_file;
  size_t Vsize;
  /////size_t file_type;
  /////flag 0,1,2,3 for eigen save
  int file_type;
  size_t end_of_file;
  int N_noi, ncur;
  int bfac_write;
  bool rotate_bfac;
  bool do_checksum;
  std::string filename;
  ////===private usage, not loaded from file head


  //inputpara(bool printlog_set = false){
  //  printlog = printlog_set;
  //}

  ~inputpara(){
    //for(unsigned int is=0;is<read_f.size();is++){
    // for(unsigned int ic=0;ic<read_f[is].size();ic++){read_f[is][ic].resize(0);}
    // read_f[is].resize(0);
    //}
    read_f.resize(0);
  }

  //int find_para(const std::string &str2, crc32_t &res){
  //  for(unsigned int is=0;is<read_f.size();is++){
  //    ////std::string str2("bSize");
  //    std::size_t found = read_f[is][0].find(str2);
  //    if(found != std::string::npos and read_f[is].size() >= 2){
  //      std::sscanf(read_f[is][1].c_str(), "%X", &res);
  //      if(printlog)if(get_node_rank_funs0() == 0)
  //        printf("  %20s %X \n", str2.c_str(), res);
  //      return 1;
  //    }
  //  }
  //  return 0;
  //}

  int find_string(const std::string &str0, const std::string &str2)
  {
    int res = 0;
    //int found = find_string(str0, str2);
    //if(found != std::string::npos and found==0){res = 1;}
    if(str0 == str2 ){res = 1;}
    return res;
  }

  //////TODO Check found == 0 correct for all cases
  int find_para(const std::string &str2, int &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      int found = find_string(read_f[is][0], str2);
      ////std::size_t found = read_f[is][0].find(str2);
      if(found == 1 and read_f[is].size() >= 2){
        res = stringtonum(read_f[is][1]);
        if(printlog)if(get_node_rank_funs0() == 0)
          printf("  %20s %10d \n",str2.c_str(), res);
        return 1;
      }
    }
    return 0;
  }

  int find_para(const std::string &str2, double &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      int found = find_string(read_f[is][0], str2);
      //std::size_t found = read_f[is][0].find(str2);
      if(found == 1 and read_f[is].size() >= 2){
        res = stringtodouble(read_f[is][1]);
        if(printlog)if(get_node_rank_funs0() == 0){
          if(res >= 1e-6){printf("  %20s %.6f \n", str2.c_str(), res);}
          if(res <  1e-6){printf("  %20s %.3e \n", str2.c_str(), res);}
        }
        return 1;
      }
    }
    return 0;
  }


  void read_geo(Geometry& geo)
  {
    std::vector<int > nv(4);
    for(int i=0;i<4;i++){nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
  }

  int find_para(const std::string &str2, std::string &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      int found = find_string(read_f[is][0], str2);
      //std::size_t found = read_f[is][0].find(str2);
      if(found == 1 and read_f[is].size() >= 2){
        if(read_f[is].size()==2){res = read_f[is][1];}else{
          res = "";
          for(unsigned int temi=1;temi<read_f[is].size();temi++){
            res += read_f[is][temi] + " ";}
        }
        ////res = read_f[is][1];
        if(printlog)if(get_node_rank_funs0() == 0)
          printf("  %20s %s \n",str2.c_str(), res.c_str());
        return 1;
      }
    }
    return 0;
  }

  
  template<typename Ty>
  int find_para(const char* str2, Ty &res){
    return find_para(std::string(str2), res);
  }

  void load_para(const char *filename, bool printlog_set = true){
    printlog = printlog_set;
    int rank = get_node_rank_funs0();
    if(rank == 0)off_file = read_input(filename, read_f);
    MPI_Bcast(&off_file, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
    ////===Bcast read_f;
    size_t sizen = 0;if(rank == 0)sizen = read_f.size();
    MPI_Bcast(&sizen   , sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
    if(rank != 0)read_f.resize(sizen);

    for(unsigned int is=0;is<read_f.size();is++)
    {
      if(rank == 0)sizen = read_f[is].size();
      MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
      if(rank != 0)read_f[is].resize(sizen);
      for(unsigned int ic=0;ic<read_f[is].size();ic++)
      {
        if(rank == 0)sizen = read_f[is][ic].size();
        MPI_Bcast(&sizen, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
        if(rank != 0)read_f[is][ic].resize(sizen);
        MPI_Bcast(&read_f[is][ic][0], sizen, MPI_CHAR, 0, MPI_COMM_WORLD);
      }
    }
    ////===Bcast read_f;

    std::string tem;
    if(printlog)if(get_node_rank_funs0() == 0)printf("========Start input \n");
    if(find_para(std::string("nx"),nx)==0)nx  = 0;
    if(find_para(std::string("ny"),ny)==0)ny  = 0;
    if(find_para(std::string("nz"),nz)==0)nz  = 0;
    if(find_para(std::string("nt"),nt)==0)nt  = 0;

    if(find_para(std::string("bSize"),bSize)==0)bSize = 32;
    if(find_para(std::string("bSum"),bSum)==0)bSum  = 512;
    if(find_para(std::string("cutN"),cutN)==0)cutN  = 8;
    if(find_para(std::string("lat"),lat)==0)lat  = std::string("24D");
    if(find_para(std::string("icfg"),icfg)==0)icfg  = 0;
    if(find_para(std::string("icfg_end"),icfg_end)==0)icfg_end  = 999999;
    if(find_para(std::string("icfg_jump"),icfg_jump)==0)icfg_jump  = 10;
    if(find_para(std::string("nprop"),nprop)==0)nprop  = 0;
    if(find_para(std::string("GPU"),GPU)==0)GPU  = 1;

    if(find_para(std::string("bini"),bini)==0)bini = -1;
    if(find_para(std::string("ncut0"),ncut0)==0)ncut0 = 30;
    if(find_para(std::string("ncut1"),ncut1)==0)ncut1 = 30;

    if(find_para(std::string("nini"),nini)==0)nini  = 0;
    if(find_para(std::string("nvec"),nvec)==0)nvec  = 0;
    if(find_para(std::string("write_mode"),write_mode)==0)write_mode  = 0;
    if(find_para(std::string("mode_dis"),mode_dis)==0)mode_dis  = 0;
    if(find_para(std::string("split_save"),split_save)==0)split_save  = 0;
    if(find_para(std::string("ndouble"),ndouble)==0)ndouble  = 200;
    if(find_para(std::string("fermion_type"),fermion_type)==0)fermion_type  = 0;
    if(find_para(std::string("gridtem"),gridtem)==0)gridtem  = 1;
    if(find_para(std::string("sparsefactor"),sparsefactor)==0)sparsefactor  = 16;
    if(find_para(std::string("lms"),lms)==0)lms  = 0;
    if(find_para(std::string("combineT"),combineT)==0)combineT  = 1;
    if(find_para(std::string("mom_cut"),mom_cut)==0)mom_cut     = 4;
    if(find_para(std::string("mode_FFT_MPI"),mode_FFT_MPI)==0)mode_FFT_MPI  = -1;
    if(find_para(std::string("Eerr"),Eerr)==0)Eerr  = 1e-11;
    if(find_para(std::string("nsave"),nsave)==0)nsave  = 0;
    if(find_para(std::string("bfac"),bfac)==0)bfac  = 0;
    if(find_para(std::string("seed"),seed)==0)seed  = 0;
    if(find_para(std::string("hyp"),hyp)==0)hyp  = 0;
    if(find_para(std::string("ionum"),ionum)==0)ionum  = 0;
    if(find_para(std::string("Link_name"),Link_name)==0)Link_name  = std::string("NONE");
    if(find_para(std::string("Ename"),Ename)==0)Ename  = std::string("NONE");
    if(find_para(std::string("Ename_Sm"),Ename_Sm)==0)Ename_Sm  = std::string("NONE");
    if(find_para(std::string("Sname"),Sname)==0)Sname  = std::string("NONE");
    if(find_para(std::string("Pname"),Pname)==0)Pname  = std::string("NONE");
    if(find_para(std::string("do_all_low"),do_all_low)==0)do_all_low  = 1;
    if(find_para(std::string("output"),output)==0)output  = std::string("NONE");
    if(find_para(std::string("output_vec"),output_vec)==0)output_vec  = std::string("NONE");

    if(find_para(std::string("eig_err"),eig_err)==0)eig_err  = 0;
    if(find_para(std::string("eig_poly_deg"),eig_poly_deg)==0)eig_poly_deg  = 0;
    if(find_para(std::string("eig_amin"),eig_amin)==0)eig_amin  = 0;

    if(find_para(std::string("Propname"),Propname)==0)Propname  = std::string("NONE");
    if(find_para(std::string("Srcname"),Srcname)==0)Srcname  = std::string("NONE");

    if(find_para(std::string("layout"),layout)==0)layout = std::string("NONE");
    if(find_para(std::string("paraI"),paraI)==0)paraI  = std::string("NONE");
    if(find_para(std::string("paraIA"),paraIA)==0)paraIA = std::string("NONE");
    if(find_para(std::string("paraIB"),paraIB)==0)paraIB = std::string("NONE");
    if(find_para(std::string("paraIC"),paraIC)==0)paraIC = std::string("NONE");

    if(find_para(std::string("paraI_2pt"),paraI_2pt)==0)paraI_2pt = std::string("NONE");
    if(find_para(std::string("paraI_3pt"),paraI_3pt)==0)paraI_3pt = std::string("NONE");
    if(find_para(std::string("sink_list"),sink_list)==0)sink_list = std::string("NONE");

    if(find_para(std::string("para_stag"),para_stag)==0)para_stag = std::string("NONE");
    if(find_para(std::string("fermion_mass"),fermion_mass)==0)fermion_mass  = 0.11;
    if(find_para(std::string("niter"),niter )==0)niter  = 100000;
    if(find_para(std::string("niter_sloppy"),niter_sloppy )==0)niter_sloppy  = 100000;
    if(find_para(std::string("cg_err"),cg_err)==0)cg_err  = 1e-8;
    if(find_para(std::string("solver_type"),solver_type)==0)solver_type  = 0;
    if(find_para(std::string("inv_deflate"),inv_deflate)==0)inv_deflate  = 0;
    if(find_para(std::string("prec_type"),prec_type)==0)prec_type  = 0;

    if(find_para(std::string("job_tag"),job_tag)==0)job_tag  = std::string("NONE");
    if(find_para(std::string("src_smear_para"),src_smear_para)==0)src_smear_para  = std::string("NONE");
    if(find_para(std::string("sink_smear_para"),sink_smear_para)==0)sink_smear_para  = std::string("NONE");
    if(find_para(std::string("save_type"),save_type)==0)save_type  = std::string("NONE");
    if(find_para(std::string("total_size"),total_size)==0)total_size  = std::string("NONE");

    if(find_para(std::string("kappa"),kappa)==0)kappa  = 0.15;
    if(find_para(std::string("clover_csw"),clover_csw)==0)clover_csw  = 1.1;

    ////temp variables for prop settings
    if(find_para(std::string("SRC_PROP_WITH_LOW"),SRC_PROP_WITH_LOW)==0)SRC_PROP_WITH_LOW  = 0;

    if(find_para(std::string("key_T"),key_T)==0)key_T  = std::string("NONE");
    if(find_para(std::string("dim_name"),dim_name)==0)dim_name  = std::string("NONE");
    if(find_para(std::string("corr_name"),corr_name)==0)corr_name  = std::string("NONE");
    if(find_para(std::string("INFO_LIST"),INFO_LIST)==0)INFO_LIST  = std::string("NONE");
    if(find_para(std::string("FILE_ENDIAN"),FILE_ENDIAN)==0)FILE_ENDIAN  = std::string("NONE");
    if(find_para(std::string("VECS_TYPE"),VECS_TYPE)==0)VECS_TYPE  = std::string("NONE");

    for(int li=0;li<1000;li++){
      std::string tem = std::string("NONE");
      char mname[500];sprintf(mname, "INFOA%02d", li);
      if(find_para(std::string(mname), tem)!=0){INFOA.push_back(tem);}
    }

    if(find_para(std::string("nmass"),nmass)==0)nmass  = 0;
    for(int mi=0;mi<nmass;mi++){
      std::string tem = std::string("NONE");
      char mname[500];sprintf(mname, "mass%02d", mi);
      if(find_para(std::string(mname), tem)==0){masses.push_back(0.2);}
      else{masses.push_back(stringtodouble(tem));}
    }

    if(find_para(std::string("read_noi"),read_noi)==0)read_noi  = 1;
    if(find_para(std::string("nsource"),nsource)==0)nsource  = 0;
    for(int si=0;si<nsource;si++){
      std::string tem = std::string("NONE");
      char sname[500];sprintf(sname, "nois%05d" , si);
      char pname[500];sprintf(pname, "prop%05d", si);
      char Gname[500];sprintf(Gname, "smea%05d", si);
      if(find_para(std::string(sname), tem)==0){srcN.push_back(std::string("NONE"));}
      else{srcN.push_back(tem);}
      if(find_para(std::string(pname), tem)==0){propN.push_back(std::string("NONE"));}
      else{propN.push_back(tem);}
      if(find_para(std::string(Gname), tem)==0){smearN.push_back(std::string("NONE"));}
      else{smearN.push_back(tem);}
    }


    if(find_para(std::string("OBJECT"),tem)==0){OBJECT = std::string("NONE");}
    else{std::vector<std::string > temL = stringtolist(tem);OBJECT = temL[0];}

    if(find_para(std::string("debuga"),debuga)==0)debuga  = 0;
    if(find_para(std::string("save_prop"),save_prop)==0)save_prop  = 0;
    if(find_para(std::string("checksum"),tem)==0){checksum  = 0;}else{std::sscanf(tem.c_str(), "%X", &checksum);}
    if(printlog)if(get_node_rank_funs0() == 0)printf("========End   input \n");

    //if(get_node_rank_funs0() == 0)printf("========sum print %12X %zu \n", checksum, off_file);

  }

  void load_para(int argc, char* argv[]){
    /////load_para("input.txt");
    std::string file = std::string("NONE");
    for (int i = 1; i < argc-1; ++i) {
        if(std::string(argv[i]) == std::string("--input")){file = std::string(argv[i+1]);}
    }
    ////printf("===load %s \n",file.c_str());
    if(file != std::string("NONE")){
      ////printf("%s \n",file.c_str());
      load_para(file.c_str());
    }
    else{load_para("input.txt");}

    for (int i = 1; i < argc-1; ++i) {
      if(std::string(argv[i]) == std::string("--icfg")){
        std::string tem = std::string(argv[i+1]);icfg = stringtonum(tem);
        if(get_node_rank_funs0() == 0)printf("==Current  %20s %10d \n","icfg", icfg);
      }
      if(std::string(argv[i]) == std::string("--lat")){
        std::string tem = std::string(argv[i+1]);lat = tem;
        if(get_node_rank_funs0() == 0)printf("==Current  %20s %20s \n","lat", lat.c_str());
      }
      if(std::string(argv[i]) == std::string("--GPU")){
        std::string tem = std::string(argv[i+1]);GPU = stringtonum(tem);
        if(get_node_rank_funs0() == 0)printf("==Current  %20s %10d \n", "GPU", GPU);
      }

      if(get_node_rank_funs0() == 0)printf("========End Current input \n");
    }

    //for(int i=1;i<argc;i++){
    //  std::string str=argv[i];
    //  std::size_t found = str.find(std::string("txt"));
    //  if(found != std::string::npos)
    //  {
    //    load_para(str.c_str());
    //    return;
    //  }
    //}

    //////if(get_file_size_MPI(std::string("input.txt").c_str()) == 0)return;
    ///load_para("input.txt");
    return;
  }

};

inline void print_time()
{
  char buf[26];
  struct tm* tm_info;
  time_t timer = time(NULL);
  tm_info = localtime(&timer);
  strftime(buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
  print0("%s ", buf);
}

inline size_t vec_head_write(inputpara &in, const char* filename, int type=-1, bool clear=true){
  ////nx,ny,nz,nt, checksum, time
  ////type -1 --> prop, type 0 --> vectors, type 1 --> eigen system
  size_t off_file = 0;
  if(qlat::get_id_node() == 0){

    char buf[26];
    struct tm* tm_info;
    time_t timer = time(NULL);
    tm_info = localtime(&timer);
    strftime(buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    FILE* filew = NULL;

    if(clear == true){filew = fopen(filename, "w+");}
    if(clear == false){
    if(get_file_size_o(filename) == 0){filew = fopen(filename, "w+");}
    else{filew = fopen(filename, "r+");}}

    fseek(filew , 0 , SEEK_SET );
    int version = 0;
    ///if(type == -1)fprintf(filew, "OBJECT BEGIN_Prop_HEAD VER %d, sink 12, src 12, zyxt, R/I \n" , version);
    ///if(type ==  0)fprintf(filew, "OBJECT BEGIN_Noise_HEAD VER %d \n", version);
    ///if(type ==  1)fprintf(filew, "OBJECT BEGIN_Eigen_HEAD VER %d \n", version);

    if(type ==  0)fprintf(filew, "OBJECT BEGIN_Vecs_HEAD VER %d \n" , version);
    if(type ==  1)fprintf(filew, "OBJECT BEGIN_Corr_HEAD VER %d \n" , version);
    fprintf(filew, "VECS_TYPE %s \n", in.VECS_TYPE.c_str());

    if(type == 0){
      fprintf(filew, "nx  %d \nny  %d \nnz  %d \nnt  %d \n", in.nx, in.ny, in.nz, in.nt);}
    if(type == 1){
      fprintf(filew, "key_T     %s \n", in.key_T.c_str());
      fprintf(filew, "dim_name  %s \n", in.dim_name.c_str());
    }

    //////if(type == -1)fprintf(filew, "nprop %d \n", in.nprop);
    if(type == 0){
      fprintf(filew, "nvec %d \n" , in.nvec);
      fprintf(filew, "bfac %d \n" , in.bfac);
    }
    fprintf(filew, "save_type    %s \n", in.save_type.c_str());
    fprintf(filew, "total_size   %s \n", in.total_size.c_str());
    fprintf(filew, "checksum %12X \n"  , in.checksum);

    fprintf(filew, "Save_Date %s \n", buf);
    qassert(in.INFO_LIST.size() < LINE_LIMIT );
    fprintf(filew, "INFO_LIST %s \n", in.INFO_LIST.c_str());

    qassert(in.INFOA.size() < 1000);
    for(unsigned int li=0;li<in.INFOA.size();li++){
      qassert(in.INFOA[li].size() < LINE_LIMIT );
      fprintf(filew, "INFOA%02d %s \n", li, in.INFOA[li].c_str());
    }

    fprintf(filew, "FILE_ENDIAN %s \n", in.FILE_ENDIAN.c_str());
    fprintf(filew, "END_OF_HEAD\n");

    off_file = ftell(filew);
    fclose(filew);
  }

  //if(get_node_rank_funs0() == 0)printf("========sum print %12X %zu \n", checksum, off_file);
  MPI_Bcast(&off_file, sizeof(size_t), MPI_CHAR, 0, get_comm());
  in.off_file = off_file;
  /////if(get_node_rank_funs0() == 0)printf("END of file %zu \n", off_file);
  return off_file;
}


///inline size_t prop_head_write(inputpara &in, const char* filename, bool clear=true){
///  return vec_head_write(in, filename, -1, clear);
///}
///inline size_t noise_head_write(inputpara &in, const char* filename, bool clear=true){
///  return vec_head_write(in, filename,  0, clear);
///}
inline size_t vecs_head_write(inputpara &in, const char* filename, bool clear=true){
  return vec_head_write(in, filename,  0, clear);
}
inline size_t corr_head_write(inputpara &in, const char* filename, bool clear=true){
  return vec_head_write(in, filename,  1, clear);
}

inline int get_save_type(const std::string save_type){
  if(save_type.c_str() == std::string("double")){return 0;}
  if(save_type.c_str() == std::string("Double")){return 0;}
  if(save_type.c_str() == std::string("float") ){return 1;}
  if(save_type.c_str() == std::string("Float") ){return 1;}
  if(save_type.c_str() == std::string("single")){return 1;}
  if(save_type.c_str() == std::string("Single")){return 1;}

  print0("Cannot find type. \n");
  qassert(false);

  return -1;
  ////  if(find_para(std::string("save_type"),save_type)==0)save_type  = std::string("NONE");
}

inline size_t string_to_size(std::string &tem_string)
{
  size_t size = 0.0;
  if(tem_string!="_NONE_")
  {
    std::sscanf(tem_string.c_str(), "%zu", &size);
  }
  return size;
}

inline std::string print_size(size_t size){
  char tem_size[500];
  sprintf(tem_size, "%zu", size_t(size));
  return std::string(tem_size);
  ////qassert(std::string(tem_size) == in.total_size);
}

template <typename Ty >
struct corr_dat
{
  std::vector<int > key_T;
  std::vector<int > c_a_t;
  std::vector<std::string> dim_name;

  int  dim;
  long total;
  ////long off;
  std::vector<Ty > dat;
  std::string corr_name;
  std::string INFO_LIST;
  std::vector<std::string > INFOA;

  inline const Ty& operator[](const long i) const {qassert(i < total); return dat[i]; }
  inline Ty& operator[](const long i) {qassert(i < total); return dat[i]; }

  //corr_dat<Ty >(bool null){
  //  std::string dimN = "NONE";
  //  std::string key("");
  //  create_dat(key, dimN);
  //}

  corr_dat<Ty >(const std::string& key, const std::string& dimN = "NONE", const std::string& corr="NONE"){
    create_dat(key, dimN, corr);
  }

  corr_dat<Ty >(const char* filename, const int node_control = 0){
    read_dat(filename, node_control);
  }

  ~corr_dat(){
    dat.resize(0);
    key_T.resize(0);
    dim_name.resize(0);
  }

  void create_dat(const std::string& key, const std::string& dimN, const std::string& corr="NONE"){
    if(sizeof(Ty) != sizeof(float) and sizeof(Ty) != sizeof(double)){qassert(false);};
    std::vector<std::string > tem = stringtolist(key);
    dim = tem.size();
    key_T.resize(dim);c_a_t.resize(dim);total = 1;
    for(LInt i=0;i<tem.size();i++){
      key_T[i] = stringtonum(tem[i]);qassert(key_T[i] != 0);
      c_a_t[i] = 0;
      total = total * key_T[i];
    }
    qassert(total >= 0 and dim >= 0);
    if(dimN != std::string("NONE")){dim_name = stringtolist(dimN);}else{
      dim_name.resize(dim);
      for(int d=0;d<dim;d++){dim_name[d] = std::string(" ");}
    }
    qassert(int(dim_name.size()) == dim);

    //////memory only on node 0
    //if(qlat::get_id_node() == 0)dat.resize(total);
    dat.resize(total);zero_Ty((Ty*) dat.data(), total, 0);
    corr_name = corr;
    INFO_LIST = std::string("NONE");
    INFOA.resize(0);
  }

  long get_off(){
    long i_num = c_a_t[0];
    for(LInt i=1;i<key_T.size();i++){
      i_num = (i_num)*key_T[i] + c_a_t[i];
    }
    qassert(i_num <= total);
    return i_num;
  }

  long get_off(std::string &site){
    std::vector<std::string > tem = stringtolist(site);
    qassert(int(tem.size()) == dim);
    for(LInt i=0;i<tem.size();i++){
      c_a_t[i] = stringtonum(tem[i]);
      qassert(c_a_t[i] < key_T[i]);
    }
    return get_off();
  }

  std::vector<int > get_site(long n){
    std::vector<int > site;site.resize(dim);
    for(int iv=0;iv<dim;iv++){site[iv] = 0;}
    long tem_i = n;
    for(int Ni=0; Ni < dim; Ni++)
    {
      long N_T = 1;
      for(int numi = Ni+1;numi < dim;numi++)
      {
        N_T = N_T*key_T[numi];
      }
      site[Ni] = tem_i/N_T;
      tem_i = tem_i%N_T;
    }
    return site;
  }

  void shift_zero(){
    c_a_t = get_site(0);
  }

  void shift_off(long off){
    long cur = get_off();
    c_a_t = get_site(cur + off);
  }

  long shift_off(std::vector<int > c_a_t_off){
    if(c_a_t_off.size() == 0){return get_off();}
    qassert(c_a_t_off.size() == (LInt) dim);
    c_a_t = c_a_t_off;
    return get_off();
  }

  void read_dat(const char* filename, const int node_control = 0){
    inputpara in;
    in.load_para(filename, false);
    qassert(in.OBJECT == std::string("BEGIN_Corr_HEAD"));
    if(in.OBJECT != std::string("BEGIN_Corr_HEAD")){print0("File %s head wrong!\n", filename);
      MPI_Barrier(get_comm());
      fflush(stdout);
      qassert(false);
    }

    ////printf("==OBJECT %s \n", OBJECT);
    size_t off_file = in.off_file;
    create_dat(in.key_T, in.dim_name);
    corr_name = in.corr_name;

    int type = get_save_type(in.save_type);
    int bsize = sizeof(double);if(type == 1){bsize=sizeof(float);}
    //char tem_size[500];
    //printf(tem_size, "%30zu", size_t(total * bsize));
    //qassert(std::string(tem_size) == in.total_size);
    qassert(size_t(total * bsize) == string_to_size(in.total_size));

    INFO_LIST = in.INFO_LIST;
    INFOA     = in.INFOA;

    crc32_t crc32_tem = 0;
    if(qlat::get_id_node()==node_control){
      FILE* file = NULL;
      file = fopen(filename, "rb");
      fseek(file , off_file, SEEK_SET );

      void* buf;
      buf = (void *)aligned_alloc_no_acc(total* bsize);

      if(type==0)write_data((double*) buf, file, total, true, false);
      if(type==1)write_data((float* ) buf, file, total, true, true );

      crc32_tem = crc32_par(buf, total * bsize);

      if(type == 0)cpy_data_thread(&dat[0], (double*)buf, total, 0);
      if(type == 1)cpy_data_thread(&dat[0], (float* )buf, total, 0);

      //for(int i=0;i<total;i++){
      //  if(type==0){dat[i] = tmpD[i];}
      //  if(type==1){dat[i] = tmpF[i];}
      //}

      free(buf);
      fclose(file);file = NULL;
    }

    MPI_Bcast(&crc32_tem, sizeof(crc32_tem), MPI_CHAR, 0, get_comm());
    qassert(crc32_tem == in.checksum);
  }

  std::string get_key_T(){
    //std::string key = std::string("");
    //for(int d=0;d<dim;d++){key += (std::string("  ") + std::to_string(key_T[d]));}
    return qlat::listtostring(key_T);
    //return key;
  }

  std::string get_dim_name(){
    //std::string dim_N = std::string("");
    //for(int d=0;d<dim;d++)(dim_N += (std::string("  ") + dim_name[d]));
    return qlat::listtostring(dim_name);
    //return dim_N;
  }

  void write_dat(const char* filename, const int node_control = 0){
    TIMER("Write corr");

    ////in.key_T, in.dim_name, in.total_size, in.checksum
    ////int type = get_save_type(save_type);
    ////in.save_type = save_type;
    inputpara in;int type = 0;
    if(sizeof(Ty) == sizeof(double)){in.save_type = std::string("Double");type = 0;}
    if(sizeof(Ty) == sizeof(float) ){in.save_type = std::string("Single");type = 1;}
    int bsize = sizeof(double);if(type == 1){bsize=sizeof(float);}

    //in.key_T = std::string("");
    //for(int d=0;d<dim;d++){in.key_T += (std::string("  ") + std::to_string(key_T[d]));}

    //in.dim_name = std::string("");
    //for(int d=0;d<dim;d++)(in.dim_name += (std::string("  ") + dim_name[d]));
    in.key_T = get_key_T();
    in.dim_name = get_dim_name();

    in.total_size = print_size(size_t(total * bsize));

    in.checksum = 0;
    ///////data for analysis is small endian
    in.FILE_ENDIAN = std::string("LITTLEENDIAN");
    in.INFO_LIST = INFO_LIST;
    in.INFOA     = INFOA;
    in.corr_name = corr_name;

    size_t off_file = corr_head_write(in, filename, true);
    //print0("dat off %30zu \n", off_file);

    crc32_t crc32_tem = 0;
    if(qlat::get_id_node()==node_control){
      FILE* file = NULL;
      file = fopen(filename, "wb");
      fseek(file , off_file, SEEK_SET );

      std::vector<char > buf;buf.resize(total * bsize);

      if(type == 0)cpy_data_thread((double*)&buf[0], &dat[0], total, 0);
      if(type == 1)cpy_data_thread((float* )&buf[0], &dat[0], total, 0);

      if(type==0)write_data((double*)&buf[0], file, total, false, false);
      if(type==1)write_data((float* )&buf[0], file, total, false, true );

      crc32_tem = crc32_par(&buf[0], total * bsize);

      fclose(file);file = NULL;
    }

    in.checksum = crc32_tem;

    size_t off_tem = corr_head_write(in, filename, false);
    qassert(off_file == off_tem);

  }

  void add_size(const int n){
    if(key_T.size() < 1){
      print0("key_T size wrong!\n");MPI_Barrier(get_comm());
      fflush(stdout);qassert(false);}

    std::vector<Ty > buf;buf.resize(dat.size());
    cpy_data_thread((Ty*) buf.data(), (Ty*) dat.data(), dat.size(), 0);

    key_T[0] += n;
    total = 1; 
    for(LInt i=0;i<key_T.size();i++){total = total * key_T[i];}

    dat.resize(total);
    cpy_data_thread((Ty*) dat.data(), (Ty*) buf.data(), buf.size(), 0);
    zero_Ty((Ty*) &dat[buf.size()], total - buf.size(), 0);
  }

  template<typename Ta>
  void write_corr(Ta* src, const long size, int mode_copy = 0 ){
    TIMER("write_corr");
    ///if(size > total){abort_r("Write size too larg. \n");}
    const long cur = get_off();
    long long double_size = size;

    const int is_double = get_data_type_is_double<Ta >();
    if( is_double){double_size = size * sizeof(Ta)/sizeof(double);}
    if(!is_double){double_size = size * sizeof(Ta)/sizeof(float );}

    if(double_size + cur >  total){ 
      if(key_T.size() < 1){
        print0("key_T size wrong!\n");MPI_Barrier(get_comm());
        fflush(stdout);qassert(false);}
      long each = total/key_T[0];long base = key_T[0];
      int n = (double_size + cur + each - 1) / (each) - base;
      add_size(n);
    }

    qassert(mode_copy == 0 or mode_copy == 3);
    if( is_double){cpy_data_thread((Ty*) &dat[cur], (double*) src, double_size, mode_copy);}
    if(!is_double){cpy_data_thread((Ty*) &dat[cur], (float* ) src, double_size, mode_copy);}
    shift_off(double_size);
  }

  void set_zero(){
    zero_Ty(&dat[0], dat.size());
    shift_zero();
  }

  void print_info(){
    if(qlat::get_id_node()==0){
      printf("===Corr %s, dim %d, mem size %.3e MB \n", 
            corr_name.c_str(), dim, total * sizeof(double)*1.0/(1024.0*1024.0));
      for(int d=0;d<dim;d++){
        printf("dim %30s   %d \n", dim_name[d].c_str(), key_T[d]);
      }
      if(INFO_LIST != "NONE"){printf("INFO_LIST %s \n", INFO_LIST.c_str());}
      for(unsigned int si=0;si<INFOA.size();si++){
        if(INFOA[si] != "NONE"){printf("INFO_LIST %s \n", INFOA[si].c_str());}
      }
      printf("===End Corr\n");
    }
  }

};



}

#endif

