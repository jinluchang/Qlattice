#include <qlat/qlat.h>

namespace qlat
{  //

void test()
{
  displayln_info(ssprintf("sizeof(ColorMatrix) = %d", sizeof(ColorMatrix)));
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  RngState rs;
  set_g_rand_color_matrix_field(gf, rs, 0.1);
  gf_show_info(gf);
  qmkdir("results");
  save_gauge_field(gf, "results/ckpoint.0");
  gf.init();
  load_gauge_field(gf, "results/ckpoint.0");
  gf_show_info(gf);
  FieldM<double, 1> fsum;
  fsum.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> fv = gf.get_elems_const(xl);
    double sum = 0;
    for (int m = 0; m < (int)fv.size(); ++m) {
      sum += matrix_trace(fv[m]).real();
    }
    fsum.get_elem(xl) = sum;
  }
  double gf_trace_sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    gf_trace_sum += fsum.get_elem(xl);
  }
  glb_sum(gf_trace_sum);
  displayln_info(ssprintf("INFO: gf_trace_sum = %24.17E", gf_trace_sum));
  displayln_info(ssprintf("INFO: gf_trace_avg = %24.17E",
                          gf_trace_sum / geo.total_volume() / 4.0 / 3.0));
  displayln_info(ssprintf("INFO: gf_avg_plaq = %24.17E", gf_avg_plaq(gf)));
}

}  // namespace qlat

int main(int argc, char** argv)
{
  using namespace qlat;
  begin(&argc, &argv);
  test();
  displayln_info("CHECK: finished successfully.");
  end();
  return 0;
}
