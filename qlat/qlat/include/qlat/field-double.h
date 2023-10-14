#pragma once

#include <qlat/field.h>

namespace qlat
{  //

template <class M>
void set_checkers_double(Field<M>& f)
{
  TIMER("set_checkers");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    Vector<M> v = f.get_elems(xl);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      if ((xg[0] + xg[1] + xg[2] + xg[3]) % 2 == 0)
        dv[m] = 1.0;
      else
        dv[m] = -1.0;
    }
  }
}

template <class M>
void set_complex_from_double(Field<M>& cf, const Field<double>& sf)
{
  TIMER("set_complex_from_double");
  const Geometry geo = sf.geo();
  // cf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = cf.get_elems(xl);
    Vector<Complex> cf_v((Complex*)v.data(), v.data_size() / sizeof(Complex));
    int N = cf_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < N; ++m) {
      cf_v[m] = Complex(sf.get_elem(xl, m));
    }
  });
}

template <class M>
void set_double_from_complex(Field<M>& sf, const Field<Complex>& cf)
{
  TIMER("set_double_from_complex");
  const Geometry geo = cf.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = sf.get_elems(xl);
    Vector<double> sf_v((double*)v.data(), v.data_size() / sizeof(double));
    int N = sf_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < N; ++m) {
      sf_v[m] = cf.get_elem(xl, m).real();
    }
  });
}

template <class M>
void set_abs_from_complex(Field<M>& sf, const Field<Complex>& cf)
{
  TIMER("set_mod_sq_from_complex");
  const Geometry geo = cf.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = sf.get_elems(xl);
    Vector<double> sf_v((double*)v.data(), v.data_size() / sizeof(double));
    int N = sf_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < N; ++m) {
      double r = cf.get_elem(xl, m).real();
      double i = cf.get_elem(xl, m).imag();
      sf_v[m] = std::pow(r * r + i * i, 0.5);
    }
  });
}

template <class M>
void set_ratio_double(Field<M>& sf, const Field<double>& sf1,
                             const Field<double>& sf2)
{
  TIMER("set_ratio_double");
  const Geometry geo = sf.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = sf.get_elems(xl);
    Vector<double> sf_v((double*)v.data(), v.data_size() / sizeof(double));
    int N = sf_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < N; ++m) {
      sf_v[m] = sf1.get_elem(xl, m) / sf2.get_elem(xl, m);
    }
  });
}

template <class M>
void less_than_double(Field<M>& sf1, const Field<double>& sf2,
                             Field<double>& mask)
{
  TIMER("less_than");
  const Geometry geo = sf1.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    const Vector<M> v = sf1.get_elems(xl);
    Vector<double> mask_v = mask.get_elems(xl);
    const Vector<double> sf1_v((double*)v.data(),
                               v.data_size() / sizeof(double));
    int N = sf1_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < N; ++m) {
      mask_v[m] = sf1_v[m] < sf2.get_elem(xl, m);
    }
  });
}

template <class M>
void invert_double(Field<M>& sf)
{
  TIMER("invert");
  const Geometry geo = sf.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = sf.get_elems(xl);
    Vector<double> sf_v((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < geo.multiplicity; ++m) {
      sf_v[m] = 1 / sf_v[m];
    }
  });
}

template <class M>
void multiply_double(Field<M>& sf, const Field<double>& factor)
{
  TIMER("invert");
  const Geometry geo = factor.geo();
  // sf.init(geo);
  qacc_for(index, geo.local_volume(), {
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> v = sf.get_elems(xl);
    Vector<double> sf_v((double*)v.data(), v.data_size() / sizeof(double));
    int N = sf_v.size();
    qassert(N == geo.multiplicity);
    for (int m = 0; m < geo.multiplicity; ++m) {
      sf_v[m] *= factor.get_elem(xl, m);
    }
  });
}

}  // namespace qlat
