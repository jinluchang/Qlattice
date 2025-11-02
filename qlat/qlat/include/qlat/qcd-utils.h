#pragma once

#include <qlat-utils/coordinate-d.h>
#include <qlat/qcd.h>

#include <map>
#include <set>

namespace qlat
{  //

void gf_wilson_line_no_comm(Field<ColorMatrix>& wilson_line_field,
                            const Int wilson_line_field_m,
                            const GaugeField& gf_ext,
                            const std::vector<Int>& path);

void gf_wilson_line_no_comm(Field<ColorMatrix>& wilson_line_field,
                            const Int wilson_line_field_m,
                            const GaugeField& gf_ext,
                            const std::vector<Int>& path,
                            const std::vector<Int>& path_n);

struct WilsonLinePathStop {
  Coordinate x;
  std::vector<std::vector<Int>> paths;
  Int num_origins;
  //
  WilsonLinePathStop() { num_origins = 0; }
};

struct WilsonLinePathSegment {
  Coordinate target;
  std::map<Coordinate, WilsonLinePathStop> stops;
};

struct WilsonLinePath {
  std::vector<WilsonLinePathSegment> ps;
};

void set_g_rand_anti_hermitian_matrix_field(Field<ColorMatrix>& fc,
                                            const RngState& rs,
                                            const double sigma);

void set_g_rand_color_matrix_field(Field<ColorMatrix>& fc, const RngState& rs,
                                   const double sigma, const Int n_step = 1);

void set_local_current_from_props(FieldM<WilsonMatrix, 4>& cf,
                                  const Propagator4d& prop1,
                                  const Propagator4d& prop2);

double coordinate_distance_from_wilson_line(
    const Coordinate& x, const Coordinate& target_wilson_line);

std::vector<Int> find_next_dirs(const Coordinate& loc,
                                const Coordinate& target_wilson_line);

void acc_wilson_line_path_segment(WilsonLinePathSegment& path);

WilsonLinePathSegment make_wilson_line_path_segment(const Coordinate& target);

void set_multiply_simple_wilson_line_field_partial_comm(
    FieldM<ColorMatrix, 1>& wlf, FieldM<ColorMatrix, 1>& wlf1,
    const GaugeField& gf1, const std::vector<Int>& path);

void set_multiply_wilson_line_field_partial_comm(
    FieldM<ColorMatrix, 1>& wlf, FieldM<ColorMatrix, 1>& wlf1,
    const GaugeField& gf1, const WilsonLinePathSegment& path);

void set_left_expanded_gauge_field(GaugeField& gf1, const GaugeField& gf);

ColorMatrix gf_avg_wilson_line(const GaugeField& gf,
                               const WilsonLinePath& path);

WilsonLinePath make_wilson_loop_path(const Coordinate& target_l, const Int t);

ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const Int l, const Int t);

std::vector<Coordinate> spatial_permute_direction(const Coordinate& l);

ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const Coordinate& l,
                               const Int t);

void gf_show_info(const GaugeField& gf, const Int level = 0);

}  // namespace qlat
