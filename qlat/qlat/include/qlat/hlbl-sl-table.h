#pragma once

#include <qlat/qlat.h>

namespace qlat
{  //

struct SlTable {
  Int s_limit;
  Int l_limit;
  std::vector<ComplexD> table;
  //
  void init()
  {
    s_limit = 0;
    l_limit = 0;
    clear(table);
  }
  void init(const Int s, const Int l)
  {
    s_limit = s;
    l_limit = l;
    table = std::vector<Complex>(s * l, 0.0);
  }
  void init(const Coordinate& total_site)
  {
    const Int s = std::min(total_site[0] / 2 + 10, total_site[0] + 2);
    const Int l =
        std::ceil(sqrt(distance_sq_relative_coordinate_g(total_site / 2))) + 2;
    qassert(s > 0);
    qassert(l > 0);
    init(s, l);
  }
};

inline void acc_sl_table(SlTable& t)
{
  TIMER("acc_sl_table");
  qassert((int)t.table.size() == t.s_limit * t.l_limit);
  std::vector<Complex> l_table(t.l_limit, 0.0);
  for (Int i = 0; i < t.s_limit; ++i) {
    Complex sum = 0.0;
    for (Int j = 0; j < t.l_limit; ++j) {
      sum += t.table[i * t.l_limit + j];
      l_table[j] = l_table[j] + sum;
      t.table[i * t.l_limit + j] = l_table[j];
    }
  }
}

inline std::string show_sl_table(const SlTable& t)
{
  TIMER("show_sl_table");
  std::ostringstream out;
  out << ssprintf("# %24.17E %24.17E", t.table.back().real(),
                  t.table.back().imag())
      << std::endl;
  for (Int i = 0; i < t.s_limit; ++i) {
    for (Int j = 0; j < t.l_limit; ++j) {
      if (j > 0) {
        out << " ";
      }
      out << show(t.table[i * t.l_limit + j].real());
    }
    out << std::endl;
  }
  return out.str();
}

inline void add_to_sl_table(SlTable& t, const Complex& val,
                            const Int dis_sq_min, const Int dis_sq_max)
{
  const Int l_len = std::min(t.l_limit - 1, (int)std::ceil(std::sqrt(dis_sq_max)));
  const Int s_len = std::min(t.s_limit - 1, (int)std::ceil(std::sqrt(dis_sq_min)));
  qassert(s_len < t.s_limit && l_len < t.l_limit);
  t.table[s_len * t.l_limit + l_len] += val;
}

inline void add_to_sl_table(SlTable& t, const Complex& val, const Coordinate& x,
                            const Coordinate& y, const Coordinate& z,
                            const Coordinate& total_site)
{
  const long dis_sq_xy = sqr(smod(x - y, total_site));
  const long dis_sq_xz = sqr(smod(x - z, total_site));
  const long dis_sq_yz = sqr(smod(y - z, total_site));
  const long dis_sq_min = std::min(dis_sq_xy, std::min(dis_sq_xz, dis_sq_yz));
  const long dis_sq_max = std::max(dis_sq_xy, std::max(dis_sq_xz, dis_sq_yz));
  add_to_sl_table(t, val, dis_sq_min, dis_sq_max);
}

inline const SlTable& operator*=(SlTable& t, const Complex& coef)
{
  for (Int i = 0; i < (int)t.table.size(); ++i) {
    t.table[i] *= coef;
  }
  return t;
}

inline const SlTable& operator+=(SlTable& t, const SlTable& t1)
{
  qassert(t.s_limit == t1.s_limit);
  qassert(t.l_limit == t1.l_limit);
  qassert(t.table.size() == t1.table.size());
  for (Int i = 0; i < (int)t.table.size(); ++i) {
    t.table[i] += t1.table[i];
  }
  return t;
}

inline const SlTable& operator-=(SlTable& t, const SlTable& t1)
{
  qassert(t.s_limit == t1.s_limit);
  qassert(t.l_limit == t1.l_limit);
  qassert(t.table.size() == t1.table.size());
  for (Int i = 0; i < (int)t.table.size(); ++i) {
    t.table[i] -= t1.table[i];
  }
  return t;
}

}  // namespace qlat
