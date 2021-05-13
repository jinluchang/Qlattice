#pragma once
#ifndef fft_desc_h
#define fft_desc_h

#include <stdlib.h> 
#include <time.h>
#include <stdio.h>
#include <algorithm>

#include <complex>
#include <cmath>

#include "general_funs.h"

#define gwuCtype double_complex

namespace qlat{
//void abort_rl(std::string stmp);
//void milc2kentucky_propagator_rotation(vector** cpu_single_mass_prop,bool rotate_vec);
//void kentucky2milc_propagator_rotation(vector** cpu_single_mass_prop,bool rotate_vec);
////void global_sum_all(double *value,int size);
////void global_sum_all(float *value,int size);
//void sum_all(double *value,int size);
//void sum_all(float *value,int size);
//void sum_all(int *value,int size);
//void load_eigenvalues(double_complex pcEigenvalues[],double Residues[],int iNoEig,const char* filename);
//void mem_info(const char *tag);

class fft_desc_basic
{
public:
  int noden;int rank;int Nmpi;
  int nx,ny,nz,nt,vol,Nvol;
  int inix,iniy,iniz,init;
  int Nx,Ny,Nz,Nt;
  int mx,my,mz,mt;

  //std::vector< std::vector<int> > ranged;
  qlat::vector<int> Nv,nv,mv;
  std::vector<qlat::vector<int> > Pos0;
  std::vector<qlat::vector<int>  > mi_list;

  int order_ch;
  int f0,f1,f2,f3;
  qlat::vector<int > orderN;
  std::vector<qlat::vector<LInt > > maporder_Nitoi;
  std::vector<qlat::vector<LInt > > maporder_itoNi;
  std::vector<qlat::vector<LInt > > maporder_Nitoipos;

  const qlat::Geometry* geop;

  fft_desc_basic(const qlat::Geometry& geo,int order_ch_or=0, int check_mem=1)
  {
  
    inix = 1;iniy = 1;iniz = 1;init = 1;
    Nx = 1;  Ny = 1;  Nz = 1;  Nt = 1;
    mx = 1;  my = 1;  mz = 1;  mt = 1;
    ////Order change of memory within a node
    order_ch = order_ch_or;
  
    geop = &geo;
    ///noden = desc->sites_on_node;
    ///rank = get_node_rank();
    //vol = nx*ny*nz;Nmpi = vol*nt/noden;
    Nmpi  = qlat::get_num_node();
    rank  = qlat::get_id_node();

    Nv.resize(4);nv.resize(4);mv.resize(4);
    for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
    nx = nv[0];ny = nv[1];nz = nv[2];nt = nv[3];
    Nx = Nv[0];Ny = Nv[1];Nz = Nv[2];Nt = Nv[3];
    ///Nv[0] = Nx;Nv[1] = Ny;Nv[2] = Nz;Nv[3] = Nt;
    ///nv[0] = nx;nv[1] = ny;nv[2] = nz;nv[3] = nt;
    ///nx=desc->nx;ny=desc->ny;nz=desc->nz;nt=desc->nt;

    mx = nx/Nx;my = ny/Ny;mz = nz/Nz;mt = nt/Nt;
    mv[0] = mx;mv[1] = my;mv[2] = nz;mv[3] = mt;
    vol   =  nx*ny*nz;////Only the spatial volume
    Nvol  =  Nx*Ny*Nz*Nt;
    noden =  Nx*Ny*Nz*Nt;
    print0("=Need=====Size of dim,   x %5d  y %5d  z %5d  t %5d nvol %8d \n",nx,ny,nz,nt,vol*nt);
    print0("=Need=====Size of dim,  Nx %5d Ny %5d Nz %5d Nt %5d Nvol %8d \n",Nx,Ny,Nz,Nt,Nvol);

    qlat::Coordinate ts = geo.coordinate_from_index(0);
    qlat::Coordinate gs = geo.coordinate_g_from_l(ts);
    inix = gs[0];iniy = gs[1];iniz = gs[2];init = gs[3];

    ////Do not assume 0 is the initial positions
    std::vector<int> Pos0_tem;Pos0_tem.resize(Nmpi*4);
    for(int i=0;i<Pos0_tem.size();i++){Pos0_tem[i] = 0;}
    for(int i=0;i<4;i++)Pos0_tem[rank*4 + i] = gs[i];
    sum_all_size(&Pos0_tem[0], Nmpi*4);
    Pos0.resize(Nmpi);
    for(int ri=0;ri<Nmpi;ri++)
    {
      Pos0[ri].resize(4);
      for(int i=0;i<4;i++){Pos0[ri][i] = Pos0_tem[ri*4+i];}
    }
  
    /////Default memory order, nt, nz, ny, nx
    f0 = nt;f1 = nz;f2 = ny;f3 = nx;
    orderN.resize(3);orderN[0] = 2;orderN[1] = 1;orderN[2] = 0;

    /////Define fatest order to inner loop
    if(order_ch == 1){
      int maxD = 0;int maxN = Nv[maxD];
      for(int i=0;i<3;i++){if(Nv[i]>maxN){maxD = i;maxN = Nv[i];}}
      orderN[2] = maxD;
      int ci = 0;
      for(int i=0;i<3;i++){if(i != maxD){orderN[ci] = i;ci++;}}
      f1 = nv[orderN[0]];f2 = nv[orderN[1]]; f3 = nv[orderN[2]];
    }
    print0("=Order====Size of dim,  M0 %5d M1 %5d M2 %5d, d0 %2d d1 %2d d2 %2d. \n",
            Nv[orderN[0]],Nv[orderN[1]],Nv[orderN[2]],orderN[0],orderN[1],orderN[2]);


    /////Check assumption, Initial point is the smallest point
    if(check_mem == 1){
      std::vector< std::vector<int> > ranged;
      ranged.resize(4);
      for(size_t isp= 0; isp < size_t(noden); isp++)
      {
        //////position p = desc->get_position(i,rank);
        qlat::Coordinate t0 = geo.coordinate_from_index(isp);
        qlat::Coordinate g0 = geo.coordinate_g_from_l(t0);
        for(int i=0;i<4;i++){ranged[i].push_back(g0[i]);}
      }
  
      for(int i=0;i<4;i++){
        std::sort(ranged[i].begin(), ranged[i].end());
        ranged[i].erase( unique( ranged[i].begin(), ranged[i].end() ), ranged[i].end() );
      }

      int flagC = 0;
      for(int i=0;i<4;i++){if(Nv[i] != ranged[i].size()){flagC = 1;}}
      for(int i=0;i<4;i++){if(gs[i] != ranged[i][0]){flagC = 1;}}
      sum_all_size(&flagC,1);
      if(flagC>0){
        for(int i=0;i<4;i++)print0("%5d %5d, ", Nv[i], ranged[i].size());
        print0("\n");
        for(int i=0;i<4;i++)print0("%5d %5d, ", gs[i], ranged[i][0]);
        print0("\n");
        abort_r("Layout not continuous in x! \n");
      }

      /////Check continues between nodes
      std::vector<std::vector<LInt > > map_Nitoi;
      map_Nitoi.resize(Nmpi);
      for(int ri=0;ri<Nmpi;ri++){
        map_Nitoi[ri].resize(noden);
      }

      for(size_t isp=0;isp<size_t(noden);isp++){
        qlat::Coordinate t0 = geo.coordinate_from_index(isp);
        qlat::Coordinate g0 = geo.coordinate_g_from_l(t0);
        LInt offv = ((g0[3]*nz+g0[2])*ny+g0[1])*nx+g0[0];
        map_Nitoi[rank][isp] = offv;
      }

      for(int ri=0;ri<Nmpi;ri++)MPI_Bcast(&map_Nitoi[ri][0], noden*sizeof(LInt), MPI_CHAR, ri, get_comm());

      int flag = 0;
      for(int i=0;i<Nmpi;i++)
      for(size_t isp=0;isp<size_t(noden/Nx);isp++){
        LInt inix = map_Nitoi[i][isp*Nx+ 0];
        for(int xi=0;xi<Nx;xi++)if(map_Nitoi[i][isp*Nx+xi] != inix + xi){flag=1;}
      }
      sum_all_size(&flag,1);
      if(flag>0){abort_r("Layout not continuous in x! \n");}
    }

    ////////position p;
  
    maporder_itoNi.resize(Nmpi);
    maporder_Nitoi.resize(Nmpi);
    maporder_Nitoipos.resize(Nmpi);
    std::vector<int > cur0(4),curi(4);
    for(int ri=0;ri<Nmpi;ri++)
    {
      maporder_itoNi[ri].resize(noden);
      maporder_Nitoi[ri].resize(noden);
      maporder_Nitoipos[ri].resize(noden);
      for(size_t i=0;i < size_t(noden);i++)
      {
        maporder_itoNi[ri][i] = 0;
        maporder_Nitoi[ri][i] = 0;
        maporder_Nitoipos[ri][i] = 0;
      }
    }

    cur0[0] = Pos0[rank][0];
    cur0[1] = Pos0[rank][1];
    cur0[2] = Pos0[rank][2];
    cur0[3] = Pos0[rank][3];
  
    for(size_t i= 0; i < size_t(noden); i++)
    {
      qlat::Coordinate t0 = geo.coordinate_from_index(i);
      qlat::Coordinate p = geo.coordinate_g_from_l(t0);
      ///////Take care of this when porting code
      curi[0] = p[0] - cur0[0];
      curi[1] = p[1] - cur0[1];
      curi[2] = p[2] - cur0[2];

      curi[3] = p[3] - cur0[3];
      LInt newi  =  ((curi[3]*Nv[orderN[0]]+curi[orderN[0]])*Nv[orderN[1]]+curi[orderN[1]])*Nv[orderN[2]]+curi[orderN[2]];
      maporder_itoNi[rank][i] = newi;
      maporder_Nitoi[rank][newi] = i;
  
      LInt Newvi = ((curi[3]*nv[orderN[0]]+p[orderN[0]])*nv[orderN[1]]+p[orderN[1]])*nv[orderN[2]]+p[orderN[2]];
      maporder_Nitoipos[rank][newi] = Newvi;
    }

    for(int ri=0;ri<Nmpi;ri++)
    {
      sum_all_size(&maporder_itoNi[ri][0], noden);
      sum_all_size(&maporder_Nitoi[ri][0], noden);
      sum_all_size(&maporder_Nitoipos[ri][0], noden);
    }
  
    std::vector<int> mi_list_tem;
    mi_list_tem.resize(mt*mz*my*mx);
    for(int i=0;i<mi_list_tem.size();i++){mi_list_tem[i] = 0;}
    qlat::Coordinate p_tem;
    for(int tmi=0;tmi<mt;tmi++)
    {
      for(int zi=0;zi<mz;zi++)
      for(int yi=0;yi<my;yi++)
      for(int xi=0;xi<mx;xi++)
      {
        int off = (zi*my+yi)*mx+xi;
        p_tem[0] = xi*Nx;
        p_tem[1] = yi*Ny;
        p_tem[2] = zi*Nz;
        p_tem[3] = tmi*Nt;
        int flagE = 1;
        for(int i=0;i<4;i++){if(p_tem[i] != gs[i]){flagE = 0;}}
        if(flagE == 1){mi_list_tem[tmi*mz*my*mx + off] = rank;}
      }
    }
    sum_all_size(&mi_list_tem[0], mi_list_tem.size());
    std::vector<int> cp_l;
    for(int i=0;i<cp_l.size();i++){cp_l[i] = mi_list_tem[i];}
    std::sort(cp_l.begin(), cp_l.end());
    for(int i=0;i<cp_l.size();i++){
      int v = cp_l[i];
      //if(v < 0 or v > Nmpi){abort_r("Node map Duplicated! \n");}
      if(v != i ){abort_r("Node map Duplicated! \n");}
    }

    mi_list.resize(mt);
    for(int tmi=0;tmi<mt;tmi++)
    {
      mi_list[tmi].resize(mz*my*mx);
      for(int off=0;off<mi_list[tmi].size();off++){
        mi_list[tmi][off] = mi_list_tem[tmi*mz*my*mx + off];
      }
    }
  
  }


  ~fft_desc_basic()
  {
    Pos0.resize(0);
    mi_list.resize(0);
    orderN.resize(0);
    maporder_Nitoi.resize(0);
    maporder_itoNi.resize(0);
    maporder_Nitoipos.resize(0);
  }
};

//class eigen_rotate
//{
//public:
//  //layout_minsurface_eo desc;
//  layout_minsurface_eo *desc;
//  bool typeF;
//  //std::vector< std::vector<int> > ranged;
//  //std::vector<int> b_numb_l;
//
//  /////Use this, 4 will be together
//  //Matrix<int, 3, 4, RowMajor> Arowmajor = Acolmajor;
//  //Ev0.col(i).dot(Ev1.col(i))
//
//  EigenM Mvec;  //  2*bfac*n_keep,b_size*6
//  //std::vector<Eigen::Matrix< double, Eigen::Dynamic, 100 ,Eigen::RowMajor> > MvecC;  //  2*bfac*n_keep,b_size*6
//  EigenM EProp;
//
//  //std::vector<int > n_zeros = {0,0};
//  //EigenM temj; //1,b_size*6
//  EigenM temP;
//  EigenV temj;
//  EigenV alpha;
//
//  EigenM MPItem;
//  EigenM Esource;
//  //Eigen::VectorXcd temk;
//
//  //Normal Size
//  EigenMD eval_list;
//  EigenMD alpha_list;
//
//  //3pt function
//  EigenV alpha_list3pt;
//  EigenV alpha_list3ptC;
//  EigenV alpha3pt;
//
//  EigenV alpha3pt_ker_low;
//
//  int Edouble;
//  int smear_Eigen;
//  std::vector<EigenM > MvecC;  //  2*bfac*n_keep,b_size*6
//  EigenM temjC;
//  std::vector<EigenV >alpha_listC;
//
//  std::vector<double_complex > eval_self;
//  double rho;
//  double Eeigenerror;
//
//  timer Eigen_self;
//  timer Eigen_real;
//  timer Lowmode;
//  timer Vcopy;
//  vector *src_tmp0[12];
//
//  int nx;int ny;int nz;int nt;
//  int Nx;int Ny;int Nz;int Nt;
//
//  int n_vec;
//  int b_size;
//  int bfac;
//  int num_zero;
//  int one_minus_halfD;
//  std::vector< std::vector<int> > ranged;
//  std::vector<int> bi_l;
//  std::vector<int> bp_l;
//  std::vector<int> Cur_Xn;
//
//
//
//  //Projection between gwu and Eigen
//
//  eigen_rotate(layout_minsurface_eo &desc_or,int bsize_or,int n_vec_or);
//  //void gwu2eigen(vector **src,Eigen::Matrix< std::complex< double >, Eigen::Dynamic,  Eigen::Dynamic ,Eigen::RowMajor> &Eres,int ivec);
//  //void eigen2gwu(Eigen::Matrix< std::complex< double >, Eigen::Dynamic,  Eigen::Dynamic ,Eigen::RowMajor> &Esrc,vector **res,int ivec);
//
//  int reorder_flag,reorder_flag_E;
//  void reorder_offM();
//  void reorder();
//  std::vector<int > offM_L;
//  //std::vector<int > chi_L;
//  //std::vector<int > bi_L;
//  //std::vector<int > kn_L;
//
//
//  void copygwu(vector **src,EigenM &Eres);
//
//  void gwu2eigenKY(vector **src,EigenM &Eres);
//
//  void gwu2eigen(vector **src,EigenM &Eres,int ivec);
//  void gwu2eigen(vector **src,EigenM &Eres){
//    //v.rows() << "x" << v.cols()
//    if(Eres.rows()/(2*bfac) !=12){abort_rl("Cannot transform Eigen Prop with wrong size.");}
//    for(int i=0;i<12;i++){gwu2eigen(&src[i],Eres,i);}}
//
//
//  //void eigen_rotate::eigen_reorder(EigenM &Eres);
//  //void eigen_rotate::eigen_reorder_b(EigenM &Eres);
//
//
//
//  void eigen2gwuKY(EigenM &Esrc,vector **res);
//  void eigen2gwu(EigenM &Esrc,vector **res,int ivec,int Mvecflag=0);
//  void eigen2gwu(EigenM &Esrc,vector &res,int ivec,int Mvecflag=0);
//
//  void eigen2gwu(EigenM &Esrc,vector **res,int Mvecflag=0){
//    //v.rows() << "x" << v.cols()
//    if(Esrc.rows()/(2*bfac)!=12){abort_rl("Cannot transform Eigen Prop with wrong size.");}
//    for(int i=0;i<12;i++){eigen2gwu(Esrc,&res[i],i,Mvecflag);}}
//
//
//  //void load_eigen(FILE *pfile, const parallel_io &io,int Edouble_or=0);
//  void load_eigen(char *name_eigen, const parallel_io &io,int Edouble_or=0);
//  void load_eigen(vec_eigen_pair<vector> &eigen_use, const parallel_io &io,int Edouble_or=0);
//
//  void load_eivals(char *name_eigen_evals,double rho_or,double Eerr=1e-9);
//  void load_eivals(std::vector<double_complex > &eval_self_or,double rho_or,double Eerr=1e-9);
//
//  void initiallize(std::vector<double> &mass_or,int one_minus_halfD_or);
//  void initiallize();
//  void load_eigen(int icfg,std::string &ov_evecname,double kappa,std::vector<double> &mass,int one_minus_halfD,double eigenerror,const parallel_io &io,int Edouble_or=0);
//  //void eigen_rotate::load_eigen(vec_eigen_pair<vector> &eigen_use,std::vector<double_complex > &eval_self_or,double kappa,std::vector<double> &mass,int one_minus_halfD,double eigenerror,int Edouble_or=0);
//  void load_eigen(vec_eigen_pair<vector> &eigen_use,std::vector<double_complex > &eval_self_or,double kappa,std::vector<double> &mass,int one_minus_halfD,double eigenerror,int Edouble_or=0);
//
//  void prop_L(vector **src,EigenM &Prop_vec,std::vector<double> &mass_or,int one_minus_halfD_or);
//  void prop_L(vector **src,vector **prop,std::vector<double> &mass_or,int one_minus_halfD_or);
//  void deflation_L(vector **src,vector **res);
//
//  void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,bool bcast,double_complex* ker_low);
//  void seq_L(vector **prop3,double mass_d,double rho,int_vector &map,double_complex* ker_low);
//  void get_MPItem(EigenM &MPItem,std::vector<int > &seq_l,std::vector<std::vector<std::vector<int > > > &nodemap,std::vector<std::vector<int > > &local_update);
//
//  void setup_L(vector **prop2);
//
//
//  //void prop_L(vector **src,EigenM &Prop_vec);
//  void get_Evec(int ni,vector **src);
//  void get_Evec(int ni,vector &src);
//  ~eigen_rotate();
//};


//double sum_Eigen(EigenM &Esrc,int ivec,layout_minsurface_eo *desc);
//bool test_load(lattice_desc& desc,char *name_eigen,const parallel_io& io);
//
//void gwu2eigenKY_E(vector **src, std::vector<Evector> &prop,fft_desc_basic &fd,int chiral=0);
//void eigen2gwuKY_E(std::vector<Evector> &srcE,vector **res, fft_desc_basic &fd,int chiral=0);
//void check_prop_size(std::vector<Evector > &prop);
//
//
//void gwu2eigenKY_E(noise_vectors &src, std::vector<Evector> &res,fft_desc_basic &fd);
//void eigen2gwuKY_E(std::vector<Evector> &src,noise_vectors &res, fft_desc_basic &fd);
//
//void reorder_civ(std::complex<Ftype>* src,std::complex<Ftype>* res,int civ,int biva,int sizeF,int flag=0,int size_inner=1);
//void reorder_civ(std::vector<Evector> &src,std::vector<Evector> &res_civ,fft_desc_basic &fd,int civ,int flag=0);
//
//void memcpy_omp(Ftype* res,Ftype* src,LInt size);

}

#endif
