#ifndef UTILS_TEST_UNIFIED
#define UTILS_TEST_UNIFIED

#pragma once
#include <qlat/qlat.h>

namespace qlat{

struct test_unified
{
  qlat::vector<qlat::Complex> alpha;

  test_unified(int a){
    alpha.resize(a);
    alpha.resize(0);
  }

  void initiallize_mass(int nmass,int Ns=12)
  {
    alpha.resize(2*nmass*Ns);set_zero(alpha);
  }
  void prop_L()
  {
    //TIMER("==prop_L");
    qacc_for(coff, long(alpha.size()),{
      alpha[coff] = 0.0;
    });
  

  }

  ~test_unified()
  {
    alpha.resize(0);
  }

};


}


#endif
