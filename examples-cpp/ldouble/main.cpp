#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <qlat/vector_utils/general_funs.h>

using namespace qlat;

#if !defined(__QLAT_NO_FLOAT128__)
#define CFloat __float128
#else
#define CFloat double
#endif

void simple_tests()
{
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site);

  qlat::RngState rs(1235232);
  {
    double a = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));
    double b = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));

    CFloat a128 = a;
    CFloat b128 = b;
    displayln_info(ssprintf("CHECK: a   %+.8e, b %+.8e \n", a, b));
    const int Ncase = 8;
    qlat::vector<double > ld;ld.resize(4);

    for(int casei = 0; casei < Ncase; casei++)
    {
      const Long V = 16;
      qacc_for(c, V, {
        if(c == 0){
          //ld[0] = a;
          //ld[1] = b;
          RealDD la = a;
          RealDD lb = b;
          RealDD lz = 0.0;

          // add check
          if(casei == 0)
          {
            lz = la + lb;
          }

          // sub check
          if(casei == 1)
          {
            lz = la - lb;
          }

          // * check
          if(casei == 2)
          {
            lz = la * lb;
          }

          // / check
          if(casei == 3)
          {
            lz = la / lb;
          }

          // sqrt check
          if(casei == 4)
          {
            lz = sqrtT(fabsT(la));
          }
          // sin  check
          if(casei == 5)
          {
            lz = qsin(la + lb);
          }
          // cos  check
          if(casei == 6)
          {
            lz = qcos(la + lb);
          }
          // acos check
          if(casei == 7)
          {
            lz = qacos(qsin(la + lb));
          }

          ld[0] = lz.Y();
          ld[1] = lz.X();

          //ld[0] = a;
          //ld[1] = b;
          ////ld[1] = z0;
          ////ld[3] = lz.X();
        }
      });

      CFloat z128;
      std::string ss = "CHECK: ";
      std::string ssi = "CHECK: ";
      if(casei == 0){ss += "add ";z128 = a128 + b128;}
      if(casei == 1){ss += "sub ";z128 = a128 - b128;}
      if(casei == 2){ss += "mul ";z128 = a128 * b128;}
      if(casei == 3){ss += "div ";z128 = a128 / b128;}
      #if !defined(__QLAT_NO_FLOAT128__)
      if(casei == 4){ss += "qrt ";z128 = sqrtq(fabsq(a128));}
      #else
      if(casei == 4){ss += "qrt ";z128 = std::sqrt(std::fabs(a128));}
      #endif
      if(casei == 5){ss += "sin ";z128 = qsin(double(a128 + b128));}
      if(casei == 6){ss += "cos ";z128 = qcos(double(a128 + b128));}
      if(casei == 7){ss += "acos";z128 = qacos(qsin(double(a128 + b128)));}
      double rety = (double) z128;
      double retx = (double)(z128-(CFloat)rety);
      ld[2] = rety;
      ld[3] = retx;
      qlat::vector<double >& l = ld;
      if(qlat::get_id_node() == 0){
        printf("%+.8e, %+.8e, float128 %+.8e, %+.8e, diff %+.8e, %+.8e \n", l[0], l[1], l[2], l[3], l[0]-l[2], l[1]-l[3]);
      }
      ssi += ssprintf("%+.8e %+.8e", l[0], l[1]);
      ss += ssprintf("double %+.5e", l[2]);
      displayln_info(ss);
      displayln_info(ssi);
    }
  }

  {
    double a0 = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));
    double a1 = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));
    double b0 = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));
    double b1 = (1 - 2 * int(qlat::u_rand_gen(rs)*2)) * qlat::u_rand_gen(rs) * (std::pow(10, int(qlat::u_rand_gen(rs)*10)));

    std::complex<CFloat > aC(a0, a1);//(a0, b1);
    std::complex<CFloat > bC(b0, b1);//(a0, b1);
    ComplexT<RealDD > adC(a0, a1);//(a0, b1);
    ComplexT<RealDD > bdC(b0, b1);//(a0, b1);

    ComplexT<RealDD > zdC;
    //std::complex<CFloat > zC;
    //zC  = aC  / bC;
    zdC = adC / bdC;
    zdC = adC * ComplexT<RealDD>(1.0, 0.0);

    displayln_info(ssprintf("CHECK: a   %+.8e %+.8e, b %+.8e %+.8e \n", a0, a1, b0, b1));
    const int Ncase = 5;
    qlat::vector<double > ld;ld.resize(8);

    for(int casei = 0; casei < Ncase; casei++)
    {
      const Long V = 16;
      qacc_for(c, V, {
        if(c == 0){
          //ld[0] = a;
          //ld[1] = b;
          ComplexT<RealDD > la(a0, a1);//(a0, b1);
          ComplexT<RealDD > lb(b0, b1);//(a0, b1);
          ComplexT<RealDD > lz(0.0, 0.0);

          // add check
          if(casei == 0)
          {
            lz = la + lb;
          }

          // sub check
          if(casei == 1)
          {
            lz = la - lb;
          }

          // * check
          if(casei == 2)
          {
            lz = la * lb;
          }

          // / check
          if(casei == 3)
          {
            lz = la / lb;
          }

          // / check
          if(casei == 4)
          {
            lz = ComplexT<RealDD>(qlat::qnorm(la), 0.0);
          }

          ld[0] = lz.real().Y();
          ld[1] = lz.real().X();
          ld[2] = lz.imag().Y();
          ld[3] = lz.imag().X();
        }
      });

      std::complex<CFloat > zC;
      std::string ss = "CHECK: ";
      std::string ssi = "CHECK: ";
      if(casei == 0){ss += "add ";zC = aC + bC;}
      if(casei == 1){ss += "sub ";zC = aC - bC;}
      if(casei == 2){ss += "mul ";zC = aC * bC;}
      if(casei == 3){ss += "div ";zC = aC / bC;}
      if(casei == 4){ss += "nor ";zC = aC.real() * aC.real() + aC.imag() * aC.imag();}
      CFloat zt;
      zt = zC.real();
      double rety = (double) zt;
      double retx = (double)(zt-(CFloat)rety);
      ld[4] = rety;
      ld[5] = retx;

      zt = zC.imag();
             rety = (double) zt;
             retx = (double)(zt-(CFloat)rety);
      ld[6] = rety;
      ld[7] = retx;
      qlat::vector<double >& l = ld;
      if(qlat::get_id_node() == 0){
        printf("%+.8e, %+.8e, float128 %+.8e, %+.8e, diff %+.8e, %+.8e \n"    , l[0], l[1], l[2], l[3], l[0]-l[4], l[1]-l[5]);
        printf("    %+.8e, %+.8e, float128 %+.8e, %+.8e, diff %+.8e, %+.8e \n", l[4], l[5], l[6], l[7], l[2]-l[6], l[3]-l[7]);
      }
      //ss += ssprintf("%+.5e %+.5e, %+.5e %+.5e; "    , l[0], l[1], l[2], l[3]);
      //ss += ssprintf("float128 %+.5e %+.5e, %+.5e %+.5e "    , l[4], l[5], l[6], l[7]);
      ssi += ssprintf("%+.5e %+.5e, %+.5e %+.5e"    , l[0], l[1], l[2], l[3]);
      ss += ssprintf("double %+.5e, %+.5e"    , l[4], l[6]);
      displayln_info(ss);
      displayln_info(ssi);
    }
  }

  {
    ComplexT<RealDD > a(1.0, 0.0);
    ComplexT<double  > b(0.6, 0.7);
    ComplexT<double  > c(0.0, 0.0);
    c = b;
    a = b;
    b = a;

    c -= b;
    displayln_info(ssprintf("copy %+.8e %+.8e \n", c.real(), c.real()));

  }

}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qlat-basic-structure");
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
