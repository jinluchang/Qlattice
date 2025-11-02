// utils_geo_cache.cpp
// Gen Wang
// Sep. 2025

#include <qlat/vector_utils/utils_geo_cache.h>

namespace qlat
{

// compare in bytes may have issues when geo is have undefined values due to alignment
bool Compare_geo(const Geometry& g0, const Geometry& g1, const bool compare_bytes){
  if(compare_bytes){
    const void* c0 = (void*) &g0;
    const void* c1 = (void*) &g1;
    const size_t Nd  = sizeof(Geometry);
    //
    const Int res = std::memcmp(c0, c1, Nd);
    if(res == 0){return true;}
    else{return false;}
  }else{
    Int equal = 1;
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
    if(g0.total_site()          != g1.total_site()    ){ return 0; }
    //
    return equal;
  }
}

bool Compare_geo_less(const Geometry& g0, const Geometry& g1, const bool compare_bytes){
  // simple compare of char
  if(compare_bytes){
    const void* c0 = (void*) &g0;
    const void* c1 = (void*) &g1;
    const size_t Nd  = sizeof(Geometry);
    //
    const Int res = std::memcmp(c0, c1, Nd);
    if(res < 0){return true ;}
    if(res > 0){return false;}
    //
    //for(Long i=0;i<Nd;i++){
    //  if(c0[i] < c1[i]){return true ;}
    //  if(c1[i] < c0[i]){return false;}
    //}
    return false;
  }else{
    if(g0.initialized     < g1.initialized     ){ return true; }
    if(g1.initialized     < g0.initialized     ){ return false;}
    //
    if(g0.total_site()    < g1.total_site()    ){ return true;}
    if(g1.total_site()    < g0.total_site()    ){ return false;}
    //
    if(g0.is_only_local   < g1.is_only_local   ){ return true;}
    if(g1.is_only_local   < g0.is_only_local   ){ return false;}
    //
    if(g0.eo              < g1.eo              ){ return true;}
    if(g1.eo              < g0.eo              ){ return false;}
    //
    if(g0.geon.size_node  < g1.geon.size_node  ){ return true;}
    if(g1.geon.size_node  < g0.geon.size_node  ){ return false;}
    //
    if(g0.geon.coor_node  < g1.geon.coor_node  ){ return true;}
    if(g1.geon.coor_node  < g0.geon.coor_node  ){ return false;}
    //
    if(g0.node_site       < g1.node_site       ){ return true;}
    if(g1.node_site       < g0.node_site       ){ return false;}
    //
    if(g0.node_site_expanded < g1.node_site_expanded){ return true;}
    if(g1.node_site_expanded < g0.node_site_expanded){ return false;}
    //
    if(g0.expansion_left  < g1.expansion_left  ){  return true;}
    if(g1.expansion_left  < g0.expansion_left  ){  return false;}
    //
    if(g0.expansion_right < g1.expansion_right ){  return true;}
    if(g1.expansion_right < g0.expansion_right ){  return false;}
    //
    return false;
  }
}

void Get_geo_local(const qlat::Geometry& geo, Geometry& geo_l){
  Coordinate total_site;
  for(Int i=0;i<4;i++){
    total_site[i] = geo.node_site[i] * geo.geon.size_node[i];
  }
  geo_l.init(total_site);
}

Geometry& get_geo_cache(const Coordinate& tot, box<Geometry>& geo_BOX){
  Geometry geo_tmp;
  geo_tmp.init(tot);
  geo_BOX.set(geo_tmp);
  return geo_BOX();
}

Geometry& get_geo_cache_resize(const Geometry& geo, const Coordinate& gl, const Coordinate& gr, box<Geometry>& geo_BOX){
  Geometry geo_ex = geo_resize(geo, gl, gr );
  geo_BOX.set(geo_ex);
  return geo_BOX();
}

Geometry& get_geo_cache_local(const Geometry& geo, box<Geometry>& geo_BOX){
  Geometry geo_local;
  Get_geo_local(geo, geo_local);
  geo_BOX.set(geo_local);
  return geo_BOX();
}

void geo_to_nv(const Geometry& geo, std::vector<int >& nv, std::vector<int >& Nv, std::vector<int    >& mv){
  // read geo into vectors
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}

void geo_to_nv(const Geometry& geo, vector<int >& nv, vector<int >& Nv, vector<int >& mv){
  Nv.resize(4);nv.resize(4);mv.resize(4);
  for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
}

bool operator<(const Gbox_key& x, const Gbox_key& y){
  const bool sr = Compare_geo_less(x.geo, y.geo);
  return sr;
}

inline Cache<Gbox_key, box<Geometry> >& get_Gbox_key_cache()
{
  static Cache<Gbox_key, box<Geometry> > cache("Gbox_Key", QLAT_MAX_GEO_CACHE + 3);
  return cache;
}

Geometry& get_geo_cache(const Geometry& geo)
{
  TIMERA("get_geo_cache");
  Gbox_key ekey(geo);
  if (!get_Gbox_key_cache().has(ekey)) {
    TIMERA("get_geo_cache construct");
    //if(get_id_node() == 0)printf("==GGGGG get_geo cons \n");
    Qassert(get_Gbox_key_cache().size() < QLAT_MAX_GEO_CACHE);
    box<Geometry>& gbox = get_Gbox_key_cache()[ekey];
    gbox.set_mem_type(MemType::Uvm);
    gbox.set(geo);
  }
  Geometry& geo_buf = get_Gbox_key_cache()[ekey]();
  return geo_buf;
}

Geometry& get_geo_cache(const Coordinate& Lat)
{
  Geometry geo;geo.init(Lat);
  return get_geo_cache(geo);
}

Geometry& get_geo_cache(const Geometry& geo, const Coordinate& gl, const Coordinate& gr)
{
  Geometry geo_ex = geo_resize(geo, gl, gr );
  //if(get_id_node() == 0){
  //  printf("expand test \n%s\n\n", show(geo_ex ).c_str());
  //}
  return get_geo_cache(geo_ex);
}

Geometry& get_geo_local(const Geometry& geo){
  Geometry geo_local;
  Get_geo_local(geo, geo_local);
  return get_geo_cache(geo_local);
}

void clear_geo_cache()
{
  get_Gbox_key_cache().clear();
}


}

