#include <sys/sysinfo.h>
#include "utils_io_vec.h"
#include "general_funs.h"
#include "utils_smear_vecs.h"
#include "check_fun.h"

int main(int argc, char* argv[])
{
  using namespace qlat;

  int n_node = init_mpi(&argc, &argv);

  //begin(&argc, &argv);
  inputpara in;
  in.load_para("input.txt");
  int nx,ny,nz,nt;
  nx = in.nx;
  ny = in.ny;
  nz = in.nz;
  nt = in.nt;
  Coordinate Lat(in.nx, in.ny, in.nz, in.nt);
  //int n_node =  qlat::get_num_node();
  Coordinate spreadT = guess_nodeL(n_node, Lat);

  ///MPI_Finalize();

  /////begin_comm(get_comm(), spreadT);

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(spreadT);

  //size_node_list.push_back(Coordinate(1, 1, 1,  1));
  //size_node_list.push_back(Coordinate(1, 1, 1,  2));
  //size_node_list.push_back(Coordinate(1, 1, 3,  1));
  //size_node_list.push_back(Coordinate(1, 1, 1,  4));
  ////size_node_list.push_back(Coordinate(1, 1, 3,  2));
  //size_node_list.push_back(Coordinate(1, 2, 3,  1));
  ////size_node_list.push_back(Coordinate(1, 1, 1,  8));
  //size_node_list.push_back(Coordinate(1, 2, 4,  1));
  ////size_node_list.push_back(Coordinate(1, 1, 2,  4));
  //size_node_list.push_back(Coordinate(1, 1, 1, 12));
  //size_node_list.push_back(Coordinate(1, 1, 1, 16));
  ////size_node_list.push_back(Coordinate(1, 1, 1, 24));
  //size_node_list.push_back(Coordinate(1, 1, 6,  4));
  //size_node_list.push_back(Coordinate(1, 1, 1, 32));
  //size_node_list.push_back(Coordinate(1, 1, 4, 16));
  ////size_node_list.push_back(Coordinate(1, 1, 1, 64));
  ////size_node_list.push_back(Coordinate(1, 1, 1, 48));
  ////size_node_list.push_back(Coordinate(1, 1, 1, 96));
  ////size_node_list.push_back(Coordinate(1, 1, 1,128));
  //size_node_list.push_back(Coordinate(4, 4, 8, 16));
  //size_node_list.push_back(Coordinate(4, 8, 8, 16));
  ////size_node_list.push_back(Coordinate(1, 2, 2, 16));
  ////size_node_list.push_back(Coordinate(1, 1, 2, 16));
  ////size_node_list.push_back(Coordinate(1, 2, 2, 16));
  ////size_node_list.push_back(Coordinate(2, 2, 2, 16));
  ////size_node_list.push_back(Coordinate(2, 2, 4, 16));
  ////size_node_list.push_back(Coordinate(2, 4, 4, 16));
  ////size_node_list.push_back(Coordinate(4, 4, 4, 16));

  ////begin_thread(&argc, &argv, size_node_list);
  ////begin(&argc, &argv, size_node_list);
  //qlat::end();
  //begin(&argc, &argv, size_node_list);
  //begin_comm(get_comm(), spreadT);
  begin_comm(MPI_COMM_WORLD , spreadT);
  set_GPU();
  ///set_GPU_threads();

  //fft_desc_basic fd();

  //Coordinate size_node = Coordinate(fd.mx, fd.my, fd.mz, fd.mt);
  //begin(fd.rank, size_node);
  //begin(MPI_COMM_WORLD, size_node);


  int icfg  = in.icfg;
  int ionum = in.ionum;

  ////int vini  = 0;
  int n_vec = in.nvec;

  omp_set_num_threads(omp_get_max_threads());
  print0("===nthreads %8d %8d, max %8d \n",qlat::qacc_num_threads(),omp_get_num_threads(),omp_get_max_threads());

  Coordinate total_site = Coordinate(nx, ny, nz, nt);
  Geometry geo;
  geo.init(total_site, 1); 
  fflush_MPI();

  char rbc_conf[500],prop_name[500],namep[500];
  GaugeField gf;
  //GaugeFieldT<Complexq> gf;


  gf.init(geo);

  //Propagator4d propS;
  //Propagator4d prop_s0;
  //Propagator4d prop_s1;
  Propagator4dT<Complexq > propS;
  Propagator4dT<Complexq > prop_s0;
  Propagator4dT<Complexq > prop_s1;

  propS.init(geo);
  prop_s0.init(geo);prop_s1.init(geo);

  sprintf(rbc_conf,in.Ename.c_str(), icfg);
  /////load_gauge_field(gf_gwu,rbc_conf,true);
  load_gwu_link(rbc_conf, gf);

  if(nx == 24)twist_boundary_at_boundary(gf, -0.5, 3 );
  //if(nx == 24)twist_boundary_at_boundary(gf, EIGEN_PI/1.0, 3 );
  //qacc_for(isp, long(geo1.local_volume_expanded()), {
  //  auto* res = gfF.get_elem(isp).p;
  //  auto* src =  gf.get_elem(isp).p;
  //  for(int m=0;m<3*3*4;m++)res[isp*3*3*4 + m] = src[isp*3*3*4 + m];
  //});
  //gfF = gf1;


  sprintf(prop_name,in.Pname.c_str(), icfg);
  load_gwu_prop(prop_name, propS);

  print0("%s \n",prop_name);

  if(in.paraI != "None"){
    std::vector<std::string > Li = stringtolist(in.paraI);
    print0("Li %s, size %d \n", in.paraI.c_str(),int(Li.size()) );
    fflush_MPI();
    //qassert(Li.size()%2 == 0);

    for(int si=0;si<1;si++)
    {
      int nsmear   = stringtonum(   Li[si*2+0]);
      double width = stringtodouble(Li[si*2+1]);
      print0("sn%03dsk%6.4f \n", nsmear, width);

      GaugeField gf1;
      //GaugeFieldT<Complexq> gf1;
      ///gf1.init(geo);
      set_left_expanded_gauge_field(gf1, gf);

      /////==============Normal smear
      //long Nvol = geo.local_volume();
      //qlat::vector<Complexq > gfE;gfE.resize(6*Nvol*9);
      //const int dir_limit = 3;
      //qacc_for(index,  geo.local_volume(),{
      //  for (int dir = -dir_limit; dir < dir_limit; ++dir) {
      //    const Coordinate xl = geo.coordinate_from_index(index);
      //    const ColorMatrixT<qlat::Complex > link =
      //        dir >= 0 ? gf1.get_elem(xl, dir)
      //                 : (ColorMatrixT<qlat::Complex >)matrix_adjoint(
      //                       gf1.get_elem(coordinate_shifts(xl, dir), -dir - 1));
      //    for(int ci=0; ci<9; ci++){gfE[(dir+3)*Nvol*9 + index*9 + ci] = link.p[ci];}
      //  }
      //
      //});
      //smear_propagator_gpu(propS, gfE, width, nsmear, 3);
      /////==============Normal smear

      smear_propagator_gpu4(propS, gf1, width, nsmear);
      qlat::WilsonMatrixT<Complexq >& v0 =  propS.get_elem(0);
      print0("check %.3e %.3e \n", v0(0,0).real(), v0(0,0).imag() );


      sprintf(namep, "%s.sn%03dsk%6.4f", prop_name, nsmear, width);
      print0("%s \n",namep);
      load_gwu_prop(namep, prop_s1);
      fflush_MPI();

      diff_prop(propS,prop_s1, 1e-7);
      fflush_MPI();
    

    }
  }

  fflush_MPI();
  qlat::Timer::display();

  qlat::end();
  return 0;
}



//const Coordinate expansion_left(1, 1, 1, 1);
//const Coordinate expansion_right(0, 0, 0, 0);
//const Geometry geo1 = geo_resize(gf.geo(), expansion_left, expansion_right);
//GaugeFieldT<Complexq> gfF;
//gfF.init(geo1);
////gfF.init(geo);
//qacc_for(isp, long(geo1.local_volume_expanded()), {
//  auto* res = gfF.get_elem(isp).p;
//  auto* src =  gf.get_elem(isp).p;
//  for(int m=0;m<3*3*4;m++)res[isp*3*3*4 + m] = src[isp*3*3*4 + m];
//});
//gfF = gf1;
//smear_propagator_gwu_convension(propS, width ,nsmear, zero, false);

/////mode 0, expanded
//smear_propagator_gwu_convension_cpu(propS, gf1, width ,nsmear, 0);

/////mode 1, normal



//smear_propagator1(propS, gfE, width ,nsmear);
//smear_propagator1(propS, gf1, width ,nsmear);

//CoordinateD zero(0,0,0,0);
//const double aw = 3.0*width*width/(2*nsmear);
////const double bw = width*width/(4.0*nsmear - 6.0*width*width);
////const double coef = aw - bw + aw*bw;
//const double coef = aw;
////const double fac = std::pow(aw*bw, nsmear);
//smear_propagator(propS, gf1, coef ,nsmear);

/////smear_propagator_gwu_convension(propS, gf1, width ,nsmear);
/////propS = propS*fac;
//#pragma omp parallel for
//for (long index = 0; index < geo.local_volume(); ++index) {
//  const Coordinate xl = geo.coordinate_from_index(index);
//  WilsonMatrix& wm = propS.get_elem(index);
//  wm = wm * fac;
//}

