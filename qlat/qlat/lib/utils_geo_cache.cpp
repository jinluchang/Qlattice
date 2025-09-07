// utils_geo_cache.cpp
// Gen Wang
// Sep. 2025

#include <qlat/vector_utils/utils_geo_cache.h>

namespace qlat
{

bool Compare_geo(const Geometry& g0, const Geometry& g1){
  int equal = 1;
  if(g0.initialized           != g1.initialized ){ return 0; }
  if(g0.eo                    != g1.eo ){ return 0; }
  if(g0.is_only_local         != g1.is_only_local    ){ return 0; }
  //
  if(g0.geon                  != g1.geon ){ return 0; }
  if(g0.node_site             != g1.node_site    ){ return 0; }
  if(g0.node_site_expanded    != g1.node_site_expanded    ){ return 0; }
  //
  if(g0.expansion_left        != g1.expansion_left  ){ return 0; }
  if(g0.expansion_right       != g1.expansion_right ){ return 0; }
  //
  //if(g0.total_site()    != g1.total_site()    ){ return 0; }
  //
  return equal;
}

bool Compare_geo_less(const Geometry& g0, const Geometry& g1){
  if(g0.total_site()    < g1.total_site()    ){  return true;}
  if(g1.total_site()    < g0.total_site()    ){  return false;}
  //
  if(g0.expansion_left  < g1.expansion_left  ){  return true;}
  if(g1.expansion_left  < g0.expansion_left  ){  return false;}
  //
  if(g0.expansion_right < g1.expansion_right ){  return true;}
  if(g1.expansion_right < g0.expansion_right ){  return false;}
  //
  return false;
}

void Get_geo_local(const qlat::Geometry& geo, Geometry& geo_l){
  Coordinate total_site;
  for(int i=0;i<4;i++){
    total_site[i] = geo.node_site[i] * geo.geon.size_node[i];
  }
  geo_l.init(total_site);
}

Geometry& get_geo(const Coordinate& tot, box<Geometry>& geo_BOX){
  Geometry geo_tmp;
  geo_tmp.init(tot);
  geo_BOX.set(geo_tmp);
  return geo_BOX();
}

Geometry& get_geo_resize(const Geometry& geo, const Coordinate& gl, const Coordinate& gr, box<Geometry>& geo_BOX){
  Geometry geo_ex = geo_resize(geo, gl, gr );
  geo_BOX.set(geo_ex);
  return geo_BOX();
}

Geometry& get_geo_local(const Geometry& geo, box<Geometry>& geo_BOX){
  Geometry geo_local;
  Get_geo_local(geo, geo_local);
  geo_BOX.set(geo_local);
  return geo_BOX();
}

void geo_to_nv(const Geometry& geo, std::vector<int >& nv, std::vector<int >& Nv, std::vector<int    >& mv){
  // read geo into vectors
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}

void geo_to_nv(const Geometry& geo, vector<int >& nv, vector<int >& Nv, vector<int >& mv){
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}

bool operator<(const Gbox_key& x, const Gbox_key& y){
  int sr = Compare_geo_less(x.geo(), y.geo());
  if (sr < 0) {
    return false;
  }
  if (sr > 0) {
    return true;
  }
  return false;
}

inline Cache<Gbox_key, box<Geometry> >& get_Gbox_key_cache()
{
  static Cache<Gbox_key, box<Geometry> > cache("Gbox_Key", QLAT_MAX_GEO_CACHE + 3);
  return cache;
}

Geometry& get_geo(const Geometry& geo)
{
  TIMERA("get_geo");
  Gbox_key ekey(geo);
  if (!get_Gbox_key_cache().has(ekey)) {
    TIMERA("get_geo construct");
    Qassert(get_Gbox_key_cache().size() < QLAT_MAX_GEO_CACHE);
    box<Geometry>& gbox = get_Gbox_key_cache()[ekey];
    gbox.set_mem_type(MemType::Uvm);
    gbox.set(geo);
  }
  Geometry& geo_buf = get_Gbox_key_cache()[ekey]();
  return geo_buf;
}

Geometry& get_geo(const Coordinate& Lat)
{
  Geometry geo;geo.init(Lat);
  return get_geo(geo);
}

Geometry& get_geo(const Geometry& geo, const Coordinate& gl, const Coordinate& gr)
{
  Geometry geo_ex = geo_resize(geo, gl, gr );
  return get_geo(geo_ex);
}

Geometry& get_geo_local(const Geometry& geo){
  Geometry geo_local;
  Get_geo_local(geo, geo_local);
  return get_geo(geo_local);
}

void clear_geo_cache()
{
  get_Gbox_key_cache().clear();
}


}

