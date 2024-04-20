/* Based on: Andrew Thall, Extended-Precision Floating-Point Numbers for GPU 

   Computation. Retrieved from http://andrewthall.org/papers/df64_qf128.pdf

   on 7/12/2011.
*/

#ifndef UTILS_LDOUBLE_H
#define UTILS_LDOUBLE_H

//#include "utils_float_type.h"
#include <quadmath.h>
//#include <float.h>
//#include <fenv.h>
//#include <math.h>

#include <float.h>
#include <math.h>

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
  qacc RealDD(const int&    a){y = double(a);x = 0;}
  qacc RealDD(const double& a){y = a;        x = 0;}
  qacc RealDD(const float&  a){y = a;        x = 0;}
  qacc RealDD(const RealDD&  a){x = a.x;y = a.y;}
  #if !defined(__QLAT_NO_FLOAT128__)
  inline RealDD(const __float128&  a){
    double rety = (double) a;
    double retx = (double)(a-(__float128)rety);
    y = rety;
    x = retx;
  }
  #endif

  qacc const RealDD& operator= (const RealDD& a)
  {
    x = a.x;
    y = a.y;
    return *this;
  }

  //qacc RealDD& operator= (const RealDD& a)
  //{
  //  x = a.x;
  //  y = a.y;
  //  return *this;
  //}

  qacc RealDD& operator= (const double& a)
  {
    x = 0;
    y = a;
    return *this;
  }
  
  qacc RealDD& operator= (float& a)
  {
    x = 0;
    y = a;
    return *this;
  }

  inline RealDD& operator= (__float128& a)
  {
    double rety = (double) a;
    double retx = (double)(a-(__float128)rety);
    y = rety;
    x = retx;
    return *this;
  }

  qacc RealDD operator+=(RealDD b);
  //{
  //  RealDD& a = *this;
  //  //RealDD z;

  //  double t1, t2, t3, t4, t5, e;

  //  t1 = a.y + b.y;

  //  t2 = t1 - a.y;

  //  t3 = (a.y - (t1 - t2)) + (b.y - t2);

  //  t4 = a.x + b.x;

  //  t2 = t4 - a.x;

  //  t5 = (a.x - (t4 - t2)) + (b.x - t2);

  //  t3 = t3 + t4;

  //  t4 = t1 + t3;

  //  t3 = t3 + (t1 - t4);

  //  t3 = t3 + t5;

  //  y = e = t4 + t3;

  //  x = t3 + (e - t4);

  //  return *this;

  //  //*this = *this + a;
  //  //return  *this;
  //}

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

  operator double() const { return Y() + X(); }
  operator  float() const { return Y() + X(); }

  private:
  double y;
  double x;

};

inline __float128 copy_to_float128(const RealDD& a)
{
  __float128 a128 = (__float128) a.X() + (__float128) a.Y();
  return a128;
}

qacc RealDD minus(const RealDD& a)
{
  RealDD z;
  z.Y() = -1 * a.Y();
  z.X() = -1 * a.X();
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

  //double t1, t2, t3, t4, t5, e;

  //t1 = a.y + b.y;

  //t2 = t1 - a.y;

  //t3 = (a.y - (t1 - t2)) + (b.y - t2);

  //t4 = a.x + b.x;

  //t2 = t4 - a.x;

  //t5 = (a.x - (t4 - t2)) + (b.x - t2);

  //t3 = t3 + t4;

  //t4 = t1 + t3;

  //t3 = t3 + (t1 - t4);

  //t3 = t3 + t5;

  //z.y = e = t4 + t3;

  //z.x = t3 + (e - t4);

  //return z;
}

qacc RealDD operator-(RealDD a, RealDD b)
{
  RealDD z = a;
  z -= b;
  //RealDD bm = minus(b);
  //RealDD z = a + bm;
  return  z;
}

/* Take full advantage of FMA. Only 8 DP operations */
qacc double __Fma_rn(double x, double y, double z)
{
  #if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__)
  return fma(x , y , z);
  #else
  return __fma_rn(x, y, z);
  #endif
}

qacc double __Dmul_rn(double x, double y)
{
  #if !defined(__CUDA_ARCH__) && !defined(__HIP_ARCH__)
  return fmul(x , y);
  #else
  return __dmul_rn(x, y);
  #endif
}

qacc RealDD RealDD::operator+=(RealDD b)
{
  RealDD& a = *this;
  //RealDD z;

  double t1, t2, t3, t4, t5, e;

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
  double e;
  RealDD t;

  t.y = __Dmul_rn (a.y, b.y);     /* prevent FMA-merging */

  t.x = __Fma_rn  (a.y, b.y, -t.y);

  t.x = __Fma_rn  (a.x, b.x,  t.x); /* remove this for even higher performance */

  t.x = __Fma_rn  (a.y, b.x,  t.x);

  t.x = __Fma_rn  (a.x, b.y,  t.x);

  a.y = e = t.y + t.x;

  a.x = (t.y - e) + t.x;

  return *this;
  #endif
}

qacc RealDD  operator* (RealDD x, RealDD y)
{
  RealDD z = x;
  z *= y;
  return z;
  
  //#ifndef __CUDA_ARCH__
  ////on CPU with __float128
  //__float128 a128 = copy_to_float128(x);
  //__float128 b128 = copy_to_float128(y);
  //__float128 z128;
  //z128 = a128 * b128;
  //RealDD z(z128);
  //return z;
  //#else
  //double e;
  //RealDD t, z;

  //t.y = __Dmul_rn (x.y, y.y);     /* prevent FMA-merging */

  //t.x = __Fma_rn (x.y, y.y, -t.y);

  //t.x = __Fma_rn (x.x, y.x, t.x); /* remove this for even higher performance */

  //t.x = __Fma_rn (x.y, y.x, t.x);

  //t.x = __Fma_rn (x.x, y.y, t.x);

  //z.y = e = t.y + t.x;

  //z.x = (t.y - e) + t.x;

  //return z;
  //#endif
}

/* Take full advantage of the FMA (fused-multiply add) capability of the GPU.

   Straighforward long-hand division that treats each double as a "digit".

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
  a = z128;
  return *this;
  #else

  //RealDD z;

  double c, cc, up;

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
 
  //#ifndef __CUDA_ARCH__
  ////on CPU with __float128
  //__float128 a128 = copy_to_float128(a);
  //__float128 b128 = copy_to_float128(b);
  //__float128 z128;
  //z128 = a128 / b128;
  //RealDD z(z128);
  //return z;
  //#else

  //RealDD z;

  //double c, cc, up;

  //c = a.y / b.y;     /* compute most significant quotient "digit" */

  //cc = fma (b.y, -c, a.y); 

  //cc = fma (b.x, -c, cc);

  //cc = cc + a.x;     /* remainder after most significant quotient "digit" */

  //cc = cc / b.y;     /* compute least significant quotient "digit" */

  //z.y = up = c + cc;

  //z.x = (c - up) + cc;

  //return z;
  //#endif
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

  double e;

  RealDD t, z;

  t.Y() = std::sqrt (a.Y());

  t.X() = __Fma_rn (t.Y(), -t.Y(), a.Y());

  t.X() = t.X() + a.X();

  t.X() = t.X() / (t.Y() + t.Y());

  z.Y() = e = t.Y() + t.X();

  z.X() = (t.Y() - e) + t.X();

  return z;
  #endif
}

}

#endif
