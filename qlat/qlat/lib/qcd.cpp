#include <qlat/qcd.h>

namespace qlat
{  //

RealD gf_avg_spatial_plaq(const GaugeField& gf)
{
  return gf_avg_spatial_plaq<RealD>(gf);
}

RealD gf_avg_plaq(const GaugeField& gf) { return gf_avg_plaq<RealD>(gf); }

RealD gf_avg_link_trace(const GaugeField& gf)
{
  return gf_avg_link_trace<RealD>(gf);
}

}  // namespace qlat
