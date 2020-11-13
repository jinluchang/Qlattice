#pragma once

#include <qlat/qcd.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

struct GaugeMomentum : FieldM<ColorMatrix, 4> {
  virtual const std::string& cname()
  {
    static const std::string s = "GaugeMomentum";
    return s;
  }
};

inline bool metropolis_accept(double& accept_prob, const double delta_h,
                              const RngState& rs_)
// only compute at get_id_node() == 0
// broad_cast the result to all nodes
{
  TIMER_VERBOSE("metropolis_accept");
  double flag = 0.0;
  accept_prob = 0;
  if (get_id_node() == 0) {
    if (delta_h <= 0.0) {
      flag = 1.0;
      accept_prob = 1.0;
    } else {
      RngState rs = rs_;
      const double rand_num = u_rand_gen(rs, 0.0, 1.0);
      accept_prob = std::exp(-delta_h);
      if (rand_num <= accept_prob) {
        flag = 1.0;
      }
    }
  }
  bcast(get_data_one_elem(accept_prob));
  bcast(get_data_one_elem(flag));
  return flag > 0.5;
}

inline void set_rand_gauge_momentum(GaugeMomentum& gm, const double sigma,
                                    const RngState& rs)
//  Creates a field of antihermitian 3x3 complex matrices with each complex
//  element drawn at random from a gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  TIMER_VERBOSE("set_rand_gauge_momentum");
  set_g_rand_anti_hermitian_matrix_field(gm, rs, sigma);
}

}  // namespace qlat
