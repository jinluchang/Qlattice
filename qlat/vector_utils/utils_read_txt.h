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

inline std::vector<std::string > stringtolist(std::string &tem_string)
{
  std::istringstream iss(tem_string);
  std::vector<std::string> results((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
  return results;
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


inline unsigned int get_node_rank_funs0()
{
  int rank;
  //MPI_Comm_rank(get_comm(), &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

#define print0 if(qlat::get_id_node() == 0) printf

inline unsigned long get_file_size_o(const char *filename)
{
  std::ifstream File(filename);
  if(!File.is_open()){if(qlat::get_id_node() == 0) printf("file is not exist\n");return 0;}
  unsigned long Begin = File.tellg();
  File.seekg(0, std::ios_base::end);
  unsigned long End = File.tellg();
  File.close();
  return End-Begin;
}

inline size_t get_file_size_MPI(const char *filename)
{
  size_t sizen = 0;
  if(qlat::get_id_node()==0){
    std::ifstream File(filename);
    if(!File.is_open()){if(qlat::get_id_node() == 0) printf("file is not exist\n");sizen = 0;}
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
  if(qlat::get_id_node()==0){
    qassert(sizeof(Ty) == sizeof(float) or sizeof(Ty) == sizeof(double));
    int bsize = sizeof(double);
    if(single_file == true){bsize = sizeof(float);}

    bool Rendian = true;

    //int size = 0;
    //if(read==false)size = dat.size();
    //if(read==true){size_t sizeF = get_file_size_o(filename);size = sizeF/bsize;dat.resize(size);}

    ////Set buf with size
    char* buf=NULL;
    ////buf = new char[size*sizeof(double)];
    buf = (char *)malloc(size* bsize);

    ////Open file
    ////FILE* file = NULL;
    //if(read==false)file = fopen(filename, "wb");
    //if(read==true )file = fopen(filename, "rb");

    /////Switch endian of the file write
    if(read==false){
      if(single_file == false)cpy_data_thread((double*)buf, &dat[0], size, 1);
      if(single_file == true )cpy_data_thread((float* )buf, &dat[0], size, 1);
      /////memcpy(&buf[0],&dat[0],size*sizeof(double));
      if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
    }

    if(read==false)fwrite(&buf[0], 1, size*bsize, file);
    if(read==true ) fread(&buf[0], 1, size*bsize, file);

    /////Switch endian of the file write
    if(read==true ){
      if(Rendian == true )if( is_big_endian_gwu())switchendian((char*)&buf[0], size, bsize);
      ////memcpy(&dat[0],&buf[0],size*sizeof(double));
      if(single_file == false)cpy_data_thread(&dat[0], (double*)buf, size, 1);
      if(single_file == true )cpy_data_thread(&dat[0], (float* )buf, size, 1);
    }

    ////delete []buf;
    free(buf);
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
  if(get_file_size_o(filename) == 0){return 0;}
  FILE* filer = fopen(filename, "r");
  //////Can only be 1000 length string
  char sTemp[1001],tem[1001];
  ///std::string s0(sTemp);
  if (filer == NULL){printf("Error opening file");return 0;}

  int count_line = 0;
  ////while(!feof(filer))
  /////maximum input line 5000
  for(int i=0;i<5000;i++)
  {
    if(fgets(tem, 1001, filer) == NULL){printf("Binary file or file line too long! \n");break;};
    /////If the file is binary
    if(std::string(tem).size() >= 1000){printf("Binary file or file line too long! \n");break;}
    /////If the file is not binary
    ///printf("line %d %s", count_line, tem);
    if(std::string(tem).size() >= 2){
      sprintf(sTemp,"%s",tem);
      std::string s0(sTemp);
      if(s0 == std::string("END_OF_HEAD\n")){break;}
      std::vector<std::string > resv = stringtolist(s0);
      read_f.push_back(resv);
      count_line += 1;
    }
    if(feof(filer)){break;}
  }
  if(count_line == 5000){printf("Binary file or file line too long! \n");}
  size_t off_file = ftell(filer);
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
  int save_prop;

  int bini, ncut0, ncut1;

  int nx;
  int ny;
  int nz;
  int nt;
  int nvec;
  int bfac;
  int ionum;
  std::string Ename;
  std::string Pname;

  std::string paraI;

  int nprop;
  std::string Propname;
  std::string Srcname;
  std::string output;
  int nmass;
  std::vector<double > masses;

  int debuga;
  bool printlog;

  ////Double, Single .....
  std::string OBJECT;

  std::string save_type;
  std::string total_size;

  std::string key_T;
  std::string dim_name;
  std::string corr_name;
  std::string INFO_LIST;
  std::string VECS_TYPE;

  size_t off_file;
  crc32_t checksum;


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

  int find_para(const std::string &str2, int &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
        res = stringtonum(read_f[is][1]);
        if(printlog)if(get_node_rank_funs0() == 0)
          printf("  %20s %10d \n",str2.c_str(), res);
        return 1;
      }
    }
    return 0;
  }

  int find_para(const std::string &str2, std::string &res){
    for(unsigned int is=0;is<read_f.size();is++){
      ////std::string str2("bSize");
      std::size_t found = read_f[is][0].find(str2);
      if(found != std::string::npos and read_f[is].size() >= 2){
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
    if(find_para(std::string("nprop"),nprop)==0)nprop  = 0;

    if(find_para(std::string("bini"),bini)==0)bini = -1;
    if(find_para(std::string("ncut0"),ncut0)==0)ncut0 = 30;
    if(find_para(std::string("ncut1"),ncut1)==0)ncut1 = 30;

    if(find_para(std::string("nvec"),nvec)==0)nvec  = 0;
    if(find_para(std::string("bfac"),bfac)==0)bfac  = 0;
    if(find_para(std::string("ionum"),ionum)==0)ionum  = 0;
    if(find_para(std::string("Ename"),Ename)==0)Ename  = std::string("NONE");
    if(find_para(std::string("Pname"),Pname)==0)Pname  = std::string("NONE");
    if(find_para(std::string("output"),output)==0)output  = std::string("NONE");

    if(find_para(std::string("Propname"),Propname)==0)Propname  = std::string("NONE");
    if(find_para(std::string("Srcname"),Srcname)==0)Srcname  = std::string("NONE");

    if(find_para(std::string("paraI"),paraI)==0)paraI  = std::string("NONE");
    if(find_para(std::string("save_type"),save_type)==0)save_type  = std::string("NONE");
    if(find_para(std::string("total_size"),total_size)==0)total_size  = std::string("NONE");

    if(find_para(std::string("key_T"),key_T)==0)key_T  = std::string("NONE");
    if(find_para(std::string("dim_name"),dim_name)==0)dim_name  = std::string("NONE");
    if(find_para(std::string("corr_name"),corr_name)==0)corr_name  = std::string("NONE");
    if(find_para(std::string("INFO_LIST"),INFO_LIST)==0)INFO_LIST  = std::string("NONE");
    if(find_para(std::string("VECS_TYPE"),VECS_TYPE)==0)VECS_TYPE  = std::string("NONE");


    if(find_para(std::string("nmass"),nmass)==0)nmass  = 0;
    for(int mi=0;mi<nmass;mi++){
      std::string tem = std::string("NONE");
      char mname[500];sprintf(mname, "mass%02d", mi);
      if(find_para(std::string(mname), tem)==0){masses.push_back(0.2);}
      else{masses.push_back(stringtodouble(tem));}
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
    fprintf(filew, "INFO_LIST %s \n", in.INFO_LIST.c_str());
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

int get_save_type(const std::string save_type){
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

struct corr_dat{
  std::vector<int > key_T;
  std::vector<int > c_a_t;
  std::vector<std::string> dim_name;

  int  dim;
  long off;
  long total;
  std::vector<double > dat;
  std::string corr_name;

  void create_dat(std::string& key, std::string& dimN, std::string corr="NONE"){
    std::vector<std::string > tem = stringtolist(key);
    dim = tem.size();
    key_T.resize(dim);c_a_t.resize(dim);total = 1;
    for(int i=0;i<tem.size();i++){
      key_T[i] = stringtonum(tem[i]);
      c_a_t[i] = 0;
      total = total * key_T[i];
    }
    qassert(total > 0 and dim > 0);
    if(dimN != std::string("NONE")){dim_name = stringtolist(dimN);}else{
      dim_name.resize(dim);
      for(int d=0;d<dim;d++){dim_name[d] = std::string(" ");}
    }

    //////memory only on node 0
    //if(qlat::get_id_node() == 0)dat.resize(total);
    dat.resize(total);
    corr_name = corr;
  }

  long get_off(){
    long i_num = 1;
    for(int i=0;i<key_T.size();i++){
      i_num = (i_num)*key_T[i] + c_a_t[i];
    }
    qassert(i_num < total);
    return i_num;
  }

  long get_off(std::string &site){
    std::vector<std::string > tem = stringtolist(site);
    qassert(tem.size() == dim);
    for(int i=0;i<tem.size();i++){
      c_a_t[i] = stringtonum(tem[i]);
      qassert(c_a_t[i] < key_T[i]);
    }
    return get_off();
  }

  void read_dat(const char* filename){
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

    crc32_t crc32_tem = 0;
    if(qlat::get_id_node()==0){
      FILE* file = NULL;
      file = fopen(filename, "rb");
      fseek(file , off_file, SEEK_SET );

      void* buf;
      buf = (void *)malloc(total* bsize);
      ////double* tmpD = (double*) buf;
      ////float*  tmpF = (float*) buf;

      if(type==0)write_data((double*) buf, file, total, true, false);
      if(type==1)write_data((float* ) buf, file, total, true, true );

      crc32_tem = crc32_par(buf, total * bsize);

      if(type == 0)cpy_data_thread(&dat[0], (double*)buf, total, 1);
      if(type == 1)cpy_data_thread(&dat[0], (float* )buf, total, 1);

      //for(int i=0;i<total;i++){
      //  if(type==0){dat[i] = tmpD[i];}
      //  if(type==1){dat[i] = tmpF[i];}
      //}

      free(buf);
      fclose(file);file = NULL;
    }

    MPI_Bcast(&crc32_tem, sizeof(crc32_tem), MPI_CHAR, 0, get_comm());
    qassert(crc32_tem == in.checksum);

    ///if(read==false)file = fopen(filename, "wb");
    ///if(read==true )

  }

  void write_dat(const char* filename, std::string save_type = std::string("Double")){
    inputpara in;

    //in.key_T, in.dim_name, in.total_size, in.checksum

    in.save_type = save_type;
    int type = get_save_type(save_type);
    int bsize = sizeof(double);if(type == 1){bsize=sizeof(float);}

    char temT[500];
    sprintf(temT,"");
    for(int d=0;d<dim;d++){sprintf(temT, "%s  %d", temT, key_T[d]);}
    in.key_T    = std::string(temT);
    /////print0("key_T %s \n", in.key_T.c_str());

    sprintf(temT,"");
    for(int d=0;d<dim;d++){sprintf(temT, "%s  %s", temT, dim_name[d].c_str());}
    in.dim_name = std::string(temT);
    /////print0("dim %s \n", in.dim_name.c_str());

    in.total_size = print_size(size_t(total * bsize));
    /////print0("total_size %s \n", in.total_size.c_str());

    in.checksum = 0;
    in.corr_name = corr_name;

    size_t off_file = corr_head_write(in, filename, true);
    //print0("dat off %30zu \n", off_file);

    crc32_t crc32_tem = 0;
    if(qlat::get_id_node()==0){
      FILE* file = NULL;
      file = fopen(filename, "wb");
      fseek(file , off_file, SEEK_SET );

      void* buf;
      buf = (void *)malloc(total* bsize);
      //double* tmpD = (double*) buf;
      //float*  tmpF = (float* ) buf;

      if(type == 0)cpy_data_thread((double*)buf, &dat[0], total, 1);
      if(type == 1)cpy_data_thread((float* )buf, &dat[0], total, 1);

      ////for(int i=0;i<total;i++){
      ////  if(type==0){tmpD[i] = dat[i];}
      ////  if(type==1){tmpF[i] = dat[i];}
      ////}

      if(type==0)write_data((double*)buf, file, total, false, false);
      if(type==1)write_data((float* )buf, file, total, false, true );

      crc32_tem = crc32_par(buf, total * bsize);

      free(buf);
      fclose(file);file = NULL;
    }

    in.checksum = crc32_tem;

    ////qassert(crc32_tem == in.checksum);

    size_t off_tem = corr_head_write(in, filename, false);
    qassert(off_file == off_tem);

  }

  corr_dat(std::string key, std::string dimN = "NONE", std::string corr="NONE"){
    create_dat(key, dimN, corr);
  }

  corr_dat(const char* filename){
    read_dat(filename);
  }

  void print_dim(){
    if(qlat::get_id_node()==0){
      printf("===Corr %s, dim %d, mem size %.3e MB \n", 
            corr_name.c_str(), dim, total * sizeof(double)*1.0/(1024.0*1024.0));
      for(int d=0;d<dim;d++){
        printf("dim %30s   %d \n", dim_name[d].c_str(), key_T[d]);
      }
      printf("===End Corr\n");
    }
  }

  ~corr_dat(){
    dat.resize(0);
    key_T.resize(0);
    dim_name.resize(0);
  }

};



}

#endif
