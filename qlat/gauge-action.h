#pragma once

#include <qlat/config.h>

namespace qlat
{  //

struct GaugeAction {
  bool initialized;
  double beta;
  double c1;
  //
  qacc void init()
  {
    initialized = false;
    beta = 5.5;
    c1 = -0.331;
  }
  //
  qacc GaugeAction() { init(); }
  qacc GaugeAction(const double beta_, const double c1_ = 0.0)
  {
    init();
    initialized = true;
    beta = beta_;
    c1 = c1_;
  }
};

}  // namespace qlat
