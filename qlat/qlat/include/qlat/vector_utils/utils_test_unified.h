#ifndef UTILS_TEST_UNIFIED
#define UTILS_TEST_UNIFIED

#pragma once
#include <qlat/qcd.h>

namespace qlat{

struct test_unified
{
  qlat::vector<qlat::vector<qlat::ComplexD > > alpha;

  test_unified(Int a){
    qlat::displayln_info(qlat::ssprintf("construct"));
    ////alpha.resize(a);
    ////alpha.resize(0);
  }

  void initiallize_mass(Int nmass,Int Ns=12)
  {
    alpha.resize(2*nmass*Ns);
    for(Int i=0;i<alpha.size();i++){alpha[i].resize(10);set_zero(alpha[i]);}
  }
  void prop_L()
  {
    //TIMER("==prop_L");
    qlat::vector<qlat::ComplexD > b;b.resize(20);
    qacc_for(coff, Long(b.size()),{
      b[coff] = 0.0;
    });
    qlat::displayln_info(qlat::ssprintf("b call"));
 
    alpha.resize(20);
    for(Int i=0;i<alpha.size();i++){alpha[i].resize(10);set_zero(alpha[i]);}

    for(Int coff=0;coff<int(alpha.size());coff++)
    {
      alpha[coff][0] = 0.0;
    }
    qlat::displayln_info(qlat::ssprintf("a for call"));

    test_unified& f0 = *this;

    Int n_vec = 10;
    for(Int i=0;i<alpha.size();i++){
      qlat::displayln_info(qlat::ssprintf("a qacc for size %5d \n", alpha[i].size()));
    }

    auto& alpha = this -> alpha;
    qacc_for(coff, Long(alpha.size()),{
      for(Int i=0;i<n_vec;i++)alpha[coff][i] = 0.0;
      for(Int i=0;i<10;i++)alpha[coff][i] = 0.0;
      for(Int i=0;i<alpha[coff].size();i++)alpha[coff][i] = 0.0;
    });
    qlat::displayln_info(qlat::ssprintf("a qacc for call"));
  
    //{
    //TIMER("Global sum");
    //sum_all_size(reinterpret_cast<Ftype* > (&alpha[0]),2*alpha.size());
    //}
  }

  //~test_unified()
  //{
  //  alpha.resize(0);
  //}

};


}


#endif
