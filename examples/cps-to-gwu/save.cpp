#include <qlat/qlat.h>
#include "io_gwu.h"
#include "meson_contra.h"

#include "../load-select-data/compute-check-prop.h"
#include "../load-select-data/compute-chvp.h"
#include "../load-select-data/compute-meson-chvp.h"
#include "../load-select-data/compute-meson-snk-src.h"
#include "../load-select-data/compute-meson-vv-meson.h"
#include "../load-select-data/compute-meson-vv.h"
#include "../load-select-data/compute-psel-fsel-distribution.h"
#include "../load-select-data/compute-three-point-func.h"
#include "../load-select-data/compute-two-point-func.h"
#include "../load-select-data/compute-wall-src-prop-norm-ratio.h"

namespace qlat{

inline void save_prop(const std::string& job_tag, const int traj){
  int a =0;

  //return 0;
}

}

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();

  inputpara in;
  //in.load_para("input.txt");
  in.load_para(argc,argv);

  //qmkdir_info(ssprintf("analysis"));
  //
  //test();
  //
  std::vector<std::string> job_tags;
  // SADJUST ME
  //job_tags.push_back("24D");
  //job_tags.push_back("32D");
  //job_tags.push_back("24DH");
  //job_tags.push_back("32Dfine");
  //job_tags.push_back("48I");
  //job_tags.push_back("64I");

  job_tags.push_back(in.lat);

  qmkdir_info(ssprintf("data/dwf_gwu/"));

  for (int k = 0; k < (int)job_tags.size(); ++k)
  {
    const std::string& job_tag = job_tags[k];
    qmkdir_info(ssprintf("data/dwf_gwu/%s",job_tag.c_str()));

    const std::vector<int> trajs = get_data_trajs(job_tag);
    for (int itraj = 0; itraj < (int)trajs.size(); ++itraj)
    //if(trajs[itraj]==in.conf)
    if(itraj==in.conf)
    {
    //for (int itraj = 0; itraj < 1; ++itraj)
    const int traj = trajs[itraj];
    setup(job_tag);
    TIMER_VERBOSE("compute_traj");

    if (check_prop_psrc_exact(job_tag, traj)) {
      setup(job_tag, traj);
      TIMER_VERBOSE("check_all_prop_psrc_exact");
      qassert(check_sparse_parameters(job_tag, traj));
      qassert(check_point_src_info(job_tag, traj));

      const std::vector<PointInfo>& pis = get_point_src_info(job_tag, traj);
      //
      qassert(not does_file_exist_sync_node(
          get_prop_psrc_exact_path(job_tag, traj), "CHECK-FILE"));
      //
      int countsave = 0;
      for (int i = 0; i < (int)pis.size(); ++i) {
        ////check_sigint();
        check_sigterm();
        check_time_limit();

        const PointInfo& pi = pis[i];
        const Coordinate& xg = pi.xg;
        const int type = pi.type;
        const int accuracy = pi.accuracy;
        if (accuracy == 2 and type == 0) {
          qassert(get_does_prop_psrc_exact_exist(job_tag, traj, xg, type));

          const SelProp& s_prop = get_prop_psrc_exact(job_tag, traj, xg, type);
          //int tslice_src = xg_ini[3];
          int tslice_src = xg[3];

          const PointSelection& psel = get_point_selection(job_tag, traj);
          const FieldSelection& fsel = get_field_selection(job_tag, traj);

          Geometry geo0 = s_prop.geo();
          Geometry geo;
          geo.init(geo0.total_site(),1);
          geo.multiplicity=1;

          int ionum = 8;
          io_gwu io_use(geo,ionum);

          qlat::FieldM<qlat::Complex,1> noi,grid;
          noi.init(geo);grid.init(geo);
          qlat::set_zero(noi);qlat::set_zero(grid);
          std::vector<qlat::FermionField4dT<qlat::Complex> > prop_qlat;
          prop_qlat.resize(12);for(int iv=0;iv<12;iv++){prop_qlat[iv].init(geo);qlat::set_zero(prop_qlat[iv]);}

          {
            /////long g_index = geo.g_index_from_g_coordinate(xg);
            const Coordinate& xl = geo.coordinate_l_from_g(xg);
            if(geo.is_on_node(xl)){
              const long index = geo.index_from_coordinate(xl);
              noi.get_elem(index) = qlat::Complex(1.0,0.0);
            }

          }

          ////const Coordinate xl = geo.coordinate_from_index(index);
          ////const Coordinate xg = geo.coordinate_g_from_l(xl);

          const Coordinate vg = geo.total_site();
          int nt = vg[3];double nvol= double(vg[0])*vg[1]*vg[2];
          int tini = xg[3];
          std::vector<double > normT;normT.resize(nt);
          for(int ti=0;ti<normT.size();ti++){normT[ti]=0.0;}
          std::vector<double > write;write.resize(2*nt);
          for(int ti=0;ti<write.size();ti++){write[ti]=0.0;}
          for (long idx = 0; idx < fsel.n_elems; ++idx){
            const long index = fsel.indices[idx];
            double* p = (double*)&(s_prop.get_elem(idx));
            qlat::Complex *pc = (qlat::Complex*)&(s_prop.get_elem(idx));

            qlat::Complex buf(0.0,0.0);
            for(int dc1=0;dc1<12;dc1++)
            for(int dc0=0;dc0<12;dc0++)
            {
              double* s = (double*)&(prop_qlat[dc0].get_elem(index));
              s[(dc1)*2+0] = p[dc1 * 24 + dc0*2 + 0];
              s[(dc1)*2+1] = p[dc1 * 24 + dc0*2 + 1];
              //(*prop_s.vec[dc0]).data[gwu].d[dc1/3].c[dc1%3].real = p[dc1 * 24 + dc0*2 + 0];
              //(*prop_s.vec[dc0]).data[gwu].d[dc1/3].c[dc1%3].imag = p[dc1 * 24 + dc0*2 + 1];
              buf += pc[dc1*12+dc0]*qlat::qconj(pc[dc1*12+dc0]);
            }
            grid.get_elem(index) = qlat::Complex(1.0,0.0);

            Coordinate xl0 = geo.coordinate_from_index(index);
            Coordinate xg0 = geo.coordinate_g_from_l(xl0);
            int t = xg0[3];
            int toff = ((t-tini+nt)%nt);
            normT[toff] += 1.0;
            write[toff*2 + 0 ] += buf.real();
            write[toff*2 + 1 ] += buf.imag();

          }
          sum_all_size((double*) &normT[0],nt);
          sum_all_size((double*) &write[0],2*nt);

          char filename[500];
          if(in.save_prop == 1){
            sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.S",job_tag.c_str(),job_tag.c_str(),traj,countsave);
            save_gwu_noi(filename,noi ,io_use);
            sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.G",job_tag.c_str(),job_tag.c_str(),traj,countsave);
            save_gwu_noi(filename,grid,io_use);
            sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.prop",job_tag.c_str(),job_tag.c_str(),traj,countsave);
            save_gwu_prop(filename,prop_qlat,io_use);
          }

          //get_corr_pion(prop_qlat,xg,write);

          for(int ti=0;ti<nt;ti++){
            write[ti*2+0] = write[ti*2+0]*nvol/normT[ti];write[ti*2+1] = write[ti*2+1]*nvol/normT[ti];
          }

          sprintf(filename,"data/dwf_gwu/%s/rbc.%s.%06d.N%06d.dat",job_tag.c_str(),job_tag.c_str(),traj,countsave);
          write_data(write,filename);

          countsave += 1;

          //write_gwu_prop(s_prop,job_tag,traj,xg);

        }
      }
    }
    }
  }

  end();
  return 0;
}


