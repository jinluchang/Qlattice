// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2022 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <qlat-utils/qassert.h>
#include <qlat-utils/complex.h>
#include <qlat-utils/config.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/qcd-setting.h>
#include <qlat-utils/timer.h>

namespace qlat
{

template <class T>
struct CommaLoader {
  Vector<T> m;
  Int i;
  //
  CommaLoader(Vector<T> m_, Int i_)
  {
    m = m_;
    i = i_;
  }
  //
  CommaLoader operator,(const T& x)
  {
    m[i] = x;
    return CommaLoader(m, i + 1);
  }
};

// --------------------

template <Int DIMN, class T>
struct API QLAT_ALIGN(sizeof(T) * DIMN * DIMN) MatrixT {
  T p[DIMN * DIMN];
  //
  qacc T* data() { return p; }
  qacc const T* data() const { return p; }
  //
  // convert to RealD array
  qacc RealD* d() { return (RealD*)p; }
  qacc const RealD* d() const { return (const RealD*)p; }
  //
  qacc T& operator()(Int i, Int j)
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  qacc const T& operator()(Int i, Int j) const
  {
    qassert(0 <= i && i < DIMN);
    qassert(0 <= j && j < DIMN);
    return p[i * DIMN + j];
  }
  //
  qacc const MatrixT& operator+=(const MatrixT& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  qacc const MatrixT& operator-=(const MatrixT& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  qacc const MatrixT& operator*=(const MatrixT& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  qacc const MatrixT& operator*=(const T& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  qacc const MatrixT& operator/=(const T& x)
  {
    *this = *this / x;
    return *this;
  }
  //
  CommaLoader<T> operator<<(const T& x)
  {
    Vector<T> m(p, DIMN * DIMN);
    m[0] = x;
    return CommaLoader<T>(m, 1);
  }
};

template <class T>
struct API ColorMatrixT : MatrixT<NUM_COLOR, ComplexT<T>> {
  qacc ColorMatrixT() {}
  qacc ColorMatrixT(const MatrixT<NUM_COLOR, ComplexT<T>>& m) { *this = m; }
  //
  qacc ColorMatrixT& operator=(const MatrixT<NUM_COLOR, ComplexT<T>>& m)
  {
    *this = (const ColorMatrixT&)m;
    return *this;
  }
};

template <class T>
struct API WilsonMatrixT : MatrixT<4 * NUM_COLOR, ComplexT<T>> {
  qacc WilsonMatrixT() {}
  qacc WilsonMatrixT(const MatrixT<4 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = m;
  }
  //
  qacc WilsonMatrixT& operator=(const MatrixT<4 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = (const WilsonMatrixT&)m;
    return *this;
  }
};

template <class T>
struct API SpinMatrixT : MatrixT<4, ComplexT<T>> {
  qacc SpinMatrixT() {}
  qacc SpinMatrixT(const MatrixT<4, ComplexT<T>>& m) { *this = m; }
  //
  qacc SpinMatrixT& operator=(const MatrixT<4, ComplexT<T>>& m)
  {
    *this = (const SpinMatrixT&)m;
    return *this;
  }
};

template <class T>
struct API NonRelWilsonMatrixT : MatrixT<2 * NUM_COLOR, ComplexT<T>> {
  qacc NonRelWilsonMatrixT() {}
  qacc NonRelWilsonMatrixT(const MatrixT<2 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = m;
  }
  //
  qacc NonRelWilsonMatrixT& operator=(
      const MatrixT<2 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = (const NonRelWilsonMatrixT&)m;
    return *this;
  }
};

template <class T>
struct API IsospinMatrixT : MatrixT<2, ComplexT<T>> {
  qacc IsospinMatrixT() {}
  qacc IsospinMatrixT(const MatrixT<2, ComplexT<T>>& m) { *this = m; }
  //
  qacc IsospinMatrixT& operator=(const MatrixT<2, ComplexT<T>>& m)
  {
    *this = (const IsospinMatrixT&)m;
    return *this;
  }
};

template <class T, QLAT_ENABLE_IF(is_real<T>())>
struct API AdjointColorMatrixT : MatrixT<8, T> {
  qacc AdjointColorMatrixT() {}
  qacc AdjointColorMatrixT(const MatrixT<8, T>& m) { *this = m; }
  //
  qacc AdjointColorMatrixT& operator=(const MatrixT<8, T>& m)
  {
    *this = (const AdjointColorMatrixT&)m;
    return *this;
  }
};

using ColorMatrix = ColorMatrixT<Real>;

using WilsonMatrix = WilsonMatrixT<Real>;

using SpinMatrix = SpinMatrixT<Real>;

using NonRelWilsonMatrix = NonRelWilsonMatrixT<Real>;

using IsospinMatrix = IsospinMatrixT<Real>;

using AdjointColorMatrix = AdjointColorMatrixT<Real>;

using ColorMatrixD = ColorMatrixT<RealD>;

using WilsonMatrixD = WilsonMatrixT<RealD>;

using SpinMatrixD = SpinMatrixT<RealD>;

using NonRelWilsonMatrixD = NonRelWilsonMatrixT<RealD>;

using IsospinMatrixD = IsospinMatrixT<RealD>;

using AdjointColorMatrixD = AdjointColorMatrixT<RealD>;

using ColorMatrixF = ColorMatrixT<RealF>;

using WilsonMatrixF = WilsonMatrixT<RealF>;

using SpinMatrixF = SpinMatrixT<RealF>;

using NonRelWilsonMatrixF = NonRelWilsonMatrixT<RealF>;

using IsospinMatrixF = IsospinMatrixT<RealF>;

using AdjointColorMatrixF = AdjointColorMatrixT<RealF>;

// --------------------

template <Int DIMN, class T>
struct API QLAT_ALIGN(sizeof(T) * DIMN) MvectorT {
  T p[DIMN];
  //
  qacc T* data() { return p; }
  qacc const T* data() const { return p; }
  //
  // convert to RealD array
  qacc RealD* d() { return (RealD*)p; }
  qacc const RealD* d() const { return (const RealD*)p; }
  //
  qacc T& operator()(Int i)
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  qacc const T& operator()(Int i) const
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  //
  qacc const MvectorT& operator+=(const MvectorT& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  qacc const MvectorT& operator-=(const MvectorT& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  qacc const MvectorT& operator*=(const T& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  qacc const MvectorT& operator/=(const T& x)
  {
    *this = *this / x;
    return *this;
  }
  //
  CommaLoader<T> operator<<(const T& x)
  {
    Vector<T> m(p, DIMN);
    m[0] = x;
    return CommaLoader<T>(m, 1);
  }
};

template <class T>
struct API WilsonVectorT : MvectorT<4 * NUM_COLOR, ComplexT<T>> {
  qacc WilsonVectorT() {}
  qacc WilsonVectorT(const MvectorT<4 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = m;
  }
  //
  qacc WilsonVectorT& operator=(const MvectorT<4 * NUM_COLOR, ComplexT<T>>& m)
  {
    *this = (const WilsonVectorT&)m;
    return *this;
  }
};

template <class T>
struct API SpinVectorT : MvectorT<4, ComplexT<T>> {
  qacc SpinVectorT() {}
  qacc SpinVectorT(const MvectorT<4, ComplexT<T>>& m) { *this = m; }
  //
  qacc SpinVectorT& operator=(const MvectorT<4, ComplexT<T>>& m)
  {
    *this = (const SpinVectorT&)m;
    return *this;
  }
};

using WilsonVector = WilsonVectorT<Real>;

using SpinVector = SpinVectorT<Real>;

using WilsonVectorD = WilsonVectorT<RealD>;

using SpinVectorD = SpinVectorT<RealD>;

using WilsonVectorF = WilsonVectorT<RealF>;

using SpinVectorF = SpinVectorT<RealF>;

// --------------------

}  // namespace qlat
