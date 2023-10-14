#define QLAT_INSTANTIATE_CORE

#include <qlat/core.h>

namespace qlat
{  //

std::string show(const qlat::GeometryNode& geon)
{
  std::string s;
  s += ssprintf("{ initialized = %s\n", show(geon.initialized).c_str());
  s += ssprintf(", num_node    = %d\n", geon.num_node);
  s += ssprintf(", id_node     = %d\n", geon.id_node);
  s += ssprintf(", size_node   = %s\n", show(geon.size_node).c_str());
  s += ssprintf(", coor_node   = %s }", show(geon.coor_node).c_str());
  return s;
}

void Geometry::init(const Coordinate& total_site, const int multiplicity_)
{
  if (!initialized) {
    init();
    geon = get_geometry_node();
    multiplicity = multiplicity_;
    for (int i = 0; i < DIMN; ++i) {
      qassert(0 == total_site[i] % geon.size_node[i]);
      node_site[i] = total_site[i] / geon.size_node[i];
    }
    reset_node_site_expanded();
    initialized = true;
  }
}

std::string show(const qlat::Geometry& geo)
{
  std::string s;
  s += ssprintf("{ initialized   = %s\n", show(geo.initialized).c_str());
  s += ssprintf(", eo            = %s\n", show(geo.eo).c_str());
  s += ssprintf(", geon          =\n%s\n", show(geo.geon).c_str());
  s += ssprintf(", multiplicity  = %s\n", show(geo.multiplicity).c_str());
  s += ssprintf(", node_site     = %s\n", show(geo.node_site).c_str());
  s += ssprintf(", expan_left    = %s\n", show(geo.expansion_left).c_str());
  s += ssprintf(", expan_right   = %s\n", show(geo.expansion_right).c_str());
  s += ssprintf(", node_site_exp = %s\n", show(geo.node_site_expanded).c_str());
  s += ssprintf(", is_only_local = %s }", show(geo.is_only_local).c_str());
  return s;
}

void FieldSelection::init()
{
  f_rank.init();
  n_per_tslice = 0;
  prob = 0.0;
  f_local_idx.init();
  n_elems = 0;
  ranks.init();
  indices.init();
}

}  // namespace qlat
