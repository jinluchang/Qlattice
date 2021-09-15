#ifndef UTILS_TEST_UNIFIED
#define UTILS_TEST_UNIFIED

#pragma once
#include <qlat/qcd.h>

namespace qlat{

struct test_unified
{
  qlat::vector<qlat::vector<qlat::Complex > > alpha;

  test_unified(int a){
    qlat::displayln_info(qlat::ssprintf("construct"));
    ////alpha.resize(a);
    ////alpha.resize(0);
  }

  void initiallize_mass(int nmass,int Ns=12)
  {
    alpha.resize(2*nmass*Ns);
    for(int i=0;i<alpha.size();i++){alpha[i].resize(10);set_zero(alpha[i]);}
  }
  void prop_L()
  {
    //TIMER("==prop_L");
    qlat::vector<qlat::Complex > b;b.resize(20);
    qacc_for(coff, long(b.size()),{
      b[coff] = 0.0;
    });
    qlat::displayln_info(qlat::ssprintf("b call"));
 
    alpha.resize(20);
    for(int i=0;i<alpha.size();i++){alpha[i].resize(10);set_zero(alpha[i]);}

    for(int coff=0;coff<int(alpha.size());coff++)
    {
      alpha[coff][0] = 0.0;
    }
    qlat::displayln_info(qlat::ssprintf("a for call"));

    test_unified& f0 = *this;

    int n_vec = 10;
    for(int i=0;i<alpha.size();i++){
      qlat::displayln_info(qlat::ssprintf("a qacc for size %5d \n", alpha[i].size()));
    }

    auto& alpha = this -> alpha;
    qacc_for(coff, long(alpha.size()),{
      for(int i=0;i<n_vec;i++)alpha[coff][i] = 0.0;
      for(int i=0;i<10;i++)alpha[coff][i] = 0.0;
      for(int i=0;i<alpha[coff].size();i++)alpha[coff][i] = 0.0;
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
