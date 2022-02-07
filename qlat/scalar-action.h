#pragma once

#include <qlat/config.h>

namespace qlat
{  //

struct ScalarAction {
  bool initialized;
  double m_sq;
  double lmbd;
  double alpha;
  //
  qacc void init()
  {
    initialized = false;
    m_sq = 1.0;
    lmbd = 1.0;
    alpha = 0.0;
  }
  //
  qacc ScalarAction() { init(); }
  qacc ScalarAction(const double m_sq_, const double lmbd_, const double alpha_)
  {
    init();
    initialized = true;
    m_sq = m_sq_;
    lmbd = lmbd_;
    alpha = alpha_;
  }
};

}  // namespace qlat
