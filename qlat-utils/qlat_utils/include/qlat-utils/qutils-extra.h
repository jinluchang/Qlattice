#pragma once

#include <qlat-utils/core.h>

#include <vector>
#include <map>
#include <set>

namespace qlat
{  //

qacc int8_t& operator*=(int8_t& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

qacc int& operator*=(int& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

qacc Long& operator*=(Long& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

qacc double& operator*=(double& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

qacc float& operator*=(float& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

qacc char& operator*=(char& x, const ComplexD& factor)
{
  (void)x;
  (void)factor;
  assert(false);
  return x;
}

template <class M>
int get_type_precision()
// 0: Complex
// 1: double
// 2: ComplexF
// 3: float
// 4: long
// 5: int
// 6: char
{
  return 0;
}

template <>
inline int get_type_precision<ComplexF>()
{
  return 2;
}

template <>
inline int get_type_precision<double>()
{
  return 1;
}

template <>
inline int get_type_precision<float>()
{
  return 3;
}

template <class Vec>
bool is_equal_vec(const Vec& v1, const Vec& v2)
{
  const bool b = v1.size() == v2.size();
  if (not b) {
    return false;
  } else {
    const Long s = v1.size();
    for (Long i = 0; i < s; ++i) {
      if (not(v1[i] == v2[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class M>
bool operator==(const std::vector<M>& v1, const std::vector<M>& v2)
{
  return is_equal_vec(v1, v2);
}

template <class M, int N>
bool operator==(const array<M, N>& v1, const array<M, N>& v2)
{
  return is_equal_vec(v1, v2);
}

template <typename F>
double simpson(const F& f, const double a, const double b)
{
  return 1.0 / 6.0 * (f(a) + 4 * f(0.5 * (a + b)) + f(b)) * (b - a);
}

API inline int adaptive_simpson_min_level()
{
  static int level = 6;
  return level;
}

API inline int adaptive_simpson_max_level()
{
  static int level = 20;
  return level;
}

template <typename F>
double adaptive_simpson_level(const F& f, const double a, const double b,
                              const double eps, const int level)
{
  const double w = simpson(f, a, b);
  const double l = simpson(f, a, 0.5 * (a + b));
  const double r = simpson(f, 0.5 * (a + b), b);
  const double error = 1.0 / 15.0 * (l + r - w);
  const double result = l + r + error;
  if ((level >= adaptive_simpson_min_level() and std::abs(error) <= eps) or
      level >= adaptive_simpson_max_level()) {
    return result;
  } else {
    return adaptive_simpson_level(f, a, 0.5 * (a + b), eps / 2.0, level + 1) +
           adaptive_simpson_level(f, 0.5 * (a + b), b, eps / 2.0, level + 1);
  }
  return result;
}

template <typename F>
struct AdaptiveSimpsonToInf {
  ConstHandle<F> f;
  double start;
  //
  double operator()(const double x) const
  {
    if (x == 1.0) {
      return 0.0;
    } else {
      return f()(start + x / (1.0 - x)) / sqr(1.0 - x);
    }
  }
};

template <typename F>
struct AdaptiveSimpsonFromInf {
  ConstHandle<F> f;
  double end;
  //
  double operator()(const double x) const
  {
    if (x == 1.0) {
      return 0.0;
    } else {
      return f()(end - x / (1.0 - x)) / sqr(1.0 - x);
    }
  }
};

template <typename F>
double adaptive_simpson(const F& f, const double a, const double b,
                        const double eps)
{
  if (b < a) {
    return -adaptive_simpson(f, b, a, eps);
  } else if (a == b) {
    return 0.0;
  } else {
    const double inf = 1.0 / 0.0;
    if (a == -inf and b == inf) {
      return adaptive_simpson(f, a, 0.0, eps / 2.0) +
             adaptive_simpson(f, 0.0, b, eps / 2.0);
    } else if (b == inf) {
      AdaptiveSimpsonToInf<F> ff;
      ff.f.init(f);
      ff.start = a;
      return adaptive_simpson_level(ff, 0.0, 1.0, eps, 0);
    } else if (a == -inf) {
      AdaptiveSimpsonFromInf<F> ff;
      ff.f.init(f);
      ff.end = b;
      return adaptive_simpson_level(ff, 0.0, 1.0, eps, 0);
    } else {
      return adaptive_simpson_level(f, a, b, eps, 0);
    }
  }
}

qacc void split_work(Long& start, Long& size, const Long total,
                     const Long num_worker, const Long id_worker)
{
  const Long size_max = (total - 1) / num_worker + 1;
  start = std::min(id_worker * size_max, total);
  const Long stop = std::min(start + size_max, total);
  size = stop - start;
}

qacc Long find_worker(const Long idx, const Long total, const Long num_worker)
{
  const Long size_max = (total - 1) / num_worker + 1;
  return idx / size_max;
}

template <typename T>
std::vector<T> vector_drop(const std::vector<T>& vs, const Long n)
{
  if (n <= 0) {
    return vs;
  }
  const Long len = vs.size() - n;
  if (len <= 0) {
    return std::vector<T>();
  }
  std::vector<T> ret(len);
  for (Long i = 0; i < len; ++i) {
    ret[i] = vs[i + n];
  }
  return ret;
}

template <typename T>
std::vector<T> vector_block(const std::vector<T>& vs, const Long n_block)
// need =, +=, *=
{
  qassert(n_block >= 1);
  const Long size = vs.size();
  if (n_block > size) {
    return vs;
  }
  const Long block_size = size / n_block;
  qassert(block_size >= 1);
  const Long reminder = size - block_size * n_block;
  std::vector<T> ret(n_block);
  Long cur = 0;
  for (int i = 0; i < n_block; ++i) {
    Long count = 0;
    qassert(cur < size);
    ret[i] = vs[cur];
    cur += 1;
    count += 1;
    for (int j = 1; j < block_size + (i < reminder ? 1 : 0); ++j) {
      qassert(cur < size);
      ret[i] += vs[cur];
      cur += 1;
      count += 1;
    }
    ret[i] *= 1.0 / (double)count;
  }
  return ret;
}

template <typename T>
T average(const std::vector<T>& vs)
// need =, +=, *=
{
  const Long size = vs.size();
  qassert(size >= 1);
  T val;
  val = vs[0];
  for (Long i = 1; i < size; ++i) {
    val += vs[i];
  }
  val *= 1.0 / (double)size;
  return val;
}

template <typename T>
std::vector<T> jackknife(const std::vector<T>& vs)
// need =, +=, *=
{
  const Long size = vs.size();
  qassert(size >= 1);
  std::vector<T> ret(size + 1);
  ret[0] = average(vs);
  if (size == 1) {
    ret[1] = ret[0];
    return ret;
  }
  for (Long i = 0; i < size; ++i) {
    ret[i + 1] = vs[i];
    ret[i + 1] *= -1.0 / (double)size;
    ret[i + 1] += ret[0];
    ret[i + 1] *= (double)size / (double)(size - 1);
  }
  return ret;
}

template <typename T>
T jackknife_sigma(const std::vector<T>& vs)
// need =, *=, -, *, sqrt
{
  const Long size = vs.size();
  qassert(size >= 2);
  T val_sub, val_diff, val_sum, val2_sum;
  val_sub = vs[0];
  val_diff = vs[0] - val_sub;
  val_sum = val_diff;
  val2_sum = sqr(val_diff);
  for (Long i = 1; i < size; ++i) {
    val_diff = vs[i] - val_sub;
    val_sum += val_diff;
    val2_sum += sqr(val_diff);
  }
  val_sum *= 1.0 / (double)size;
  val2_sum *= 1.0 / (double)size;
  return std::sqrt((double)size * std::abs(val2_sum - sqr(val_sum)));
}

template <typename T>
T identity(const T& x)
{
  return x;
}

}  // namespace qlat
