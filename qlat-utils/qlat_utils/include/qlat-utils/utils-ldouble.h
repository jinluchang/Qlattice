/* Based on: Andrew Thall, Extended-Precision Floating-Point Numbers for GPU

   Computation. Retrieved from http://andrewthall.org/papers/df64_qf128.pdf

   on 7/12/2011.
*/

#ifndef UTILS_LDOUBLE_H
#define UTILS_LDOUBLE_H

#include <qlat-utils/config.h>
#include <qlat-utils/show.h>
#include <qlat-utils/qacc.h>

//#include "utils_float_type.h"
#if !defined(__QLAT_NO_FLOAT128__)
#include <quadmath.h>
#endif
//#include <float.h>
//#include <fenv.h>
//#include <math.h>

#include <cfloat>
#include <cmath>

namespace qlat
{

//y to double, x remaining stuff
class API RealDD{
  public:
  qacc double& X(){return x;}
  qacc double& Y(){return y;}

  qacc const double& X() const {return x;}
  qacc const double& Y() const {return y;}

  qacc RealDD() = default;
  //qacc RealDD(){x = 0;y = 0;}
  qacc RealDD(const Int&    a){y = RealD(a);x = 0;}
  qacc RealDD(const RealD& a){y = a;        x = 0;}
  qacc RealDD(const RealF&  a){y = a;        x = 0;}
  qacc RealDD(const RealDD& a)
  {
    y = a.y;
    x = a.x;
  }
#if !defined(__QLAT_NO_FLOAT128__)
  inline RealDD(const __float128& a)
  {
    RealD rety = (RealD)a;
    RealD retx = (RealD)(a - (__float128)rety);
    y = rety;
    x = retx;
  }
#endif

  qacc const RealDD& operator= (const RealDD& a)
  {
    y = a.Y();
    x = a.X();
    return *this;
  }

  //qacc RealDD& operator= (const RealDD& a)
  //{
  //  x = a.x;
  //  y = a.y;
  //  return *this;
  //}

  qacc RealDD& operator= (const RealD& a)
  {
    y = a;
    x = 0;
    return *this;
  }

  qacc RealDD& operator= (RealF& a)
  {
    y = a;
    x = 0;
    return *this;
  }

#if !defined(__QLAT_NO_FLOAT128__)
  inline RealDD& operator= (__float128& a)
  {
    RealD rety = (RealD) a;
    RealD retx = (RealD)(a-(__float128)rety);
    y = rety;
    x = retx;
    return *this;
  }
#endif

  qacc RealDD operator+=(RealDD b);
  qacc RealDD operator*=(RealDD b);
  qacc RealDD operator/=(RealDD b);

  qacc RealDD operator-()
  {
    RealDD am;
    am.y = -1.0 * y;
    am.x = -1.0 * x;
    return  am;
  }

  qacc RealDD operator-=(RealDD a)
  {
    RealDD am = a;
    am.y = -1.0 * am.y;
    am.x = -1.0 * am.x;
    *this += am;
    return  *this;
  }

  qacc operator RealD() const { return Y() + X(); }
  qacc operator  float() const { return Y() + X(); }
  //qacc operator double() const { return Y(); }
  //qacc operator  float() const { return Y(); }

  private:
  double y;
  double x;
};

#if !defined(__QLAT_NO_FLOAT128__)
inline __float128 copy_to_float128(const RealDD& a)
{
  __float128 a128 = (__float128) a.Y() + (__float128) a.X();
  return a128;
}
#endif

qacc RealDD minus(const RealDD& a)
{
  RealDD z;
  z.Y() = -1.0 * a.Y();
  z.X() = -1.0 * a.X();
  return z;
}

qacc RealDD fabsT(RealDD& a)
{
  RealDD z;
  z.Y() = std::fabs(a.Y());
  z.X() = std::fabs(a.X());
  return z;
}

qacc bool operator<(RealDD a, RealDD b)
{
  if(a.Y() < b.Y()){return true ;}
  if(a.Y() > b.Y()){return false;}
  if(a.X() < b.X()){return true ;}
  if(a.X() > b.X()){return false;}
  return false;
}

qacc RealDD operator+(RealDD a, RealDD b)
{
  RealDD z = a;
  z += b;
  return z;
}

qacc RealDD operator-(RealDD a, RealDD b)
{
  RealDD z = a;
  z -= b;
  return  z;
}

/* Take full advantage of FMA. Only 8 DP operations */
qacc RealD __Fma_rn(RealD x, RealD y, RealD z)
{
  #if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__)
  return std::fma(x , y , z);
  #else
  return __fma_rn(x, y, z);
  #endif
}

qacc RealD __Dmul_rn(RealD x, RealD y)
{
#if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__)
  return x * y;
  //return fmul(x, y);
#else
  return __dmul_rn(x, y);
#endif
}

qacc RealDD RealDD::operator+=(RealDD b)
{
  RealDD& a = *this;
  //RealDD z;
#if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__) && !defined(__QLAT_NO_FLOAT128__)
  //on CPU with __float128
  __float128 a128 = copy_to_float128(a);
  __float128 b128 = copy_to_float128(b);
  __float128 z128;
  z128 = a128 + b128;
  *this = z128;
  return *this;
#else

  RealD t1, t2, t3, t4, t5, e;

  t1 = a.y + b.y;

  t2 = t1 - a.y;

  t3 = (a.y - (t1 - t2)) + (b.y - t2);

  t4 = a.x + b.x;

  t2 = t4 - a.x;

  t5 = (a.x - (t4 - t2)) + (b.x - t2);

  t3 = t3 + t4;

  t4 = t1 + t3;

  t3 = t3 + (t1 - t4);

  t3 = t3 + t5;

  a.y = e = t4 + t3;

  a.x = t3 + (e - t4);

  return *this;

  //*this = *this + a;
  //return  *this;
#endif
}

qacc RealDD RealDD::operator*=(RealDD b)
{
  RealDD& a = *this;
#if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__) && !defined(__QLAT_NO_FLOAT128__)
  //on CPU with __float128
  __float128 a128 = copy_to_float128(a);
  __float128 b128 = copy_to_float128(b);
  __float128 z128;
  z128 = a128 * b128;
  *this = z128;
  return *this;
#else
  RealD e;
  RealDD t;

  t.y = __Dmul_rn (a.y, b.y);     /* prevent FMA-merging */

  t.x = __Fma_rn  (a.y, b.y, -t.y);

  t.x = __Fma_rn  (a.x, b.x,  t.x); /* remove this for even higher performance */

  t.x = __Fma_rn  (a.y, b.x,  t.x);

  t.x = __Fma_rn  (a.x, b.y,  t.x);

  a.y = e = t.y + t.x;

  a.x = (t.y - e) + t.x;

  return a;
#endif
}

qacc RealDD  operator* (RealDD a, RealDD b)
{
  RealDD z = a;
  z *= b;
  return z;
}

/* Take full advantage of the FMA (fused-multiply add) capability of the GPU.

   Straighforward Long-hand division that treats each double as a "digit".

*/

qacc RealDD RealDD::operator/= (RealDD b)
{
  RealDD& a = *this;
#if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__) && !defined(__QLAT_NO_FLOAT128__)
  //on CPU with __float128
  __float128 a128 = copy_to_float128(a);
  __float128 b128 = copy_to_float128(b);
  __float128 z128;
  z128 = a128 / b128;
  *this = z128;
  return *this;
#else

  //RealDD z;

   RealD c, cc, up;

  c = a.y / b.y;     /* compute most significant quotient "digit" */

  cc = fma (b.y, -c, a.y);

  cc = fma (b.x, -c, cc);

  cc = cc + a.x;     /* remainder after most significant quotient "digit" */

  cc = cc / b.y;     /* compute least significant quotient "digit" */

  a.y = up = c + cc;

  a.x = (c - up) + cc;

  return *this;
#endif
}


qacc RealDD operator/ (RealDD a, RealDD b)
{
  RealDD z = a;
  z /= b;
  return z;
}

/* Based on the binomial theorem. x = (a+b)^2 = a*a + 2*a*2b + b*b. If |b| is

   much smaller than |a|, then x ~= a*a + 2*a*b, thus b = (x - a*a) / (2*a).

*/

qacc RealDD sqrtT (const RealDD a)
{
#if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__) && !defined(__QLAT_NO_FLOAT128__)
  //on CPU with __float128
  __float128 a128 = copy_to_float128(a);
  a128 = sqrtq(a128);
  RealDD z(a128);
  return z;
#else

  RealD e, y, x;

  RealDD z;
  if(a.Y() == 0.0){
    z.Y() = std::sqrt(a.X());
    z.X() = 0.0;
    return z;
  }

  y     = std::sqrt (a.Y());

  x     = __Fma_rn (y    , -y    , a.Y());

  x     = x     + a.X();

  x     = x     / (y     + y    );

  z.Y() = e = y     + x  ;

  z.X() = (y     - e) + x    ;

  return z;
#endif
}

qacc bool qisnan(const RealDD &a)
{
  return std::isnan(a.X()) or std::isnan(a.Y());
}


}

#endif
