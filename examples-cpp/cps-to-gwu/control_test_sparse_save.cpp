#include <sys/sysinfo.h>
//////#include <qlat/qcd.h>
#include <qlat/qlat.h>
#include <qlat/selected-field-io.h>

int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 2, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);


  int nl = 8;
  int nt = 16;
  int seed = 123;

  Coordinate total_site = Coordinate(nl, nl, nl, nt);
  Geometry geo;
  geo.init(total_site, 1); 

  PointSelection pconf;
  FieldSelection fsel;

  qlat::RngState rG(seed);
  qlat::RngState rs(seed + qlat::get_id_node());
  Coordinate tem;
  for(int iv=0;iv<10;iv++){
    for(int j=0;j<4;j++){tem[j] = int(qlat::u_rand_gen(rG)*total_site[j]);}
    pconf.push_back(tem);
  }


  long n_per_tslice = long(long(nl*nl*nl)/(16));

  //set_field_selection(fsel, total_site, n_per_tslice, rs);
  set_field_selection(fsel, total_site, n_per_tslice, rs, pconf);

  std::string namew("res/test.fsel");
  
  write_field_selection(fsel, namew);

  ////===save point selection info
  namew = std::string("res/test.psel");
  save_point_selection_info(pconf, namew);
  ////===save point selection info

  qlat::FieldM<Complex, 12*12 > propM;propM.init(geo);
  //Propagator4dT<Complex > propM;propM.init(geo);
  //Complex* data = (Complex*) &propM.get_elem(0);
  Complex* data = qlat::get_data(propM).data();
  for(long isp=0;isp<geo.local_volume();isp++)
  {
    double ini = qlat::u_rand_gen(rs);
    for(int j=0;j<12*12;j++){data[isp*12*12 + j] = Complex(ini, ini/2.0);}
  }

  qlat::SelectedField<Complex > sf;sf.init(fsel, 12*12);
  set_selected_field(sf, propM, fsel);
  Coordinate new_size_node = Coordinate(1, 1, 2, 4);
  const ShuffledBitSet sbs = mk_shuffled_bitset(fsel, new_size_node);


  {
  ShuffledFieldsWriter sfw("res/test.selfield", new_size_node);
  /////write_float_from_double(sfw, "fa.float.field", sf);


  write_float_from_double(sfw, "f.float.sfield", sf, sbs);
  //write_float_from_double(sfw, "f.float.sfield", propM, sbs);
  //write(sfw, "f.float.field", propM);
  //write(sfw, "f.float.sfield", propM, sbs);

  {
    const std::string path = "res/test.selfield";
    const std::vector<std::string> fns = list_fields(path);
    for (long i = 0; i < (long)fns.size(); ++i) {
      displayln_info(ssprintf(": %5d : '%s' from '%s'.", i, fns[i].c_str(), path.c_str()));
    }
  }


  }

  {
  ShuffledFieldsReader sfr("res/test.selfield");
  //read_double_from_float(sfr, "f.float.sfield", propM);
  sf.init();
  read_double_from_float(sfr, "f.float.sfield", sbs, sf);

  //read(sfr, "f.float.sfield", propM);
  }

  //////namew = std::string("res/test.selfield");
  ////write_selected_field(sf, namew, fsel );
  ////read_selected_field(sf, namew, fsel);


  qlat::Timer::display();

  qlat::end();
  return 0;
}

