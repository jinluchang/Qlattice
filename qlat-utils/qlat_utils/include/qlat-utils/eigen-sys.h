#pragma once

#pragma GCC system_header
#pragma clang system_header

#ifdef QLAT_USE_GRID_EIGEN

#include <Grid/Eigen/Eigen>

namespace qlat
{

inline std::string get_eigen_type() { return "grid"; }

}  // namespace qlat

#else

#include <Eigen/Eigen>

namespace qlat
{

inline std::string get_eigen_type() { return "system"; }

}  // namespace qlat

#endif

// -------------------------------------------------------------------------------------
