#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

int main(int argc, char* argv[])
{
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  begin(&argc, &argv, size_node_list);
  if (argc != 4) {
    displayln_info("usage: ./field-compare float/double/long fn1 fn2");
    exit(-1);
  }
  const std::string dtype = argv[1];
  displayln_info(ssprintf("dtype = '%s'", dtype.c_str()));
  const std::string fn1 = remove_trailing_slashes(argv[2]);
  displayln_info(ssprintf("fn1 = '%s'", fn1.c_str()));
  const std::string fn2 = remove_trailing_slashes(argv[3]);
  displayln_info(ssprintf("fn2 = '%s'", fn2.c_str()));
  if (dtype == "float") {
    Field<double> f1, f2;
    read_field_double_from_float(f1, fn1);
    const double qnorm1 = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1) = %24.17E", qnorm1));
    read_field_double_from_float(f2, fn2);
    const double qnorm2 = qnorm(f2);
    displayln_info(ssprintf("qnorm(f2) = %24.17E", qnorm2));
    f1 -= f2;
    const double qnorm_diff = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1-f2) = %24.17E", qnorm_diff));
  } else if (dtype == "double") {
    Field<double> f1, f2;
    read_field_double(f1, fn1);
    const double qnorm1 = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1) = %24.17E", qnorm1));
    read_field_double(f2, fn2);
    const double qnorm2 = qnorm(f2);
    displayln_info(ssprintf("qnorm(f2) = %24.17E", qnorm2));
    f1 -= f2;
    const double qnorm_diff = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1-f2) = %24.17E", qnorm_diff));
  } else if (dtype == "long") {
    Field<long> f1, f2;
    read_field_64(f1, fn1);
    const double qnorm1 = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1) = %24.17E", qnorm1));
    read_field_64(f2, fn2);
    const double qnorm2 = qnorm(f2);
    displayln_info(ssprintf("qnorm(f2) = %24.17E", qnorm2));
    f1 -= f2;
    const double qnorm_diff = qnorm(f1);
    displayln_info(ssprintf("qnorm(f1-f2) = %24.17E", qnorm_diff));
  }
  Timer::display();
  end();
  return 0;
}
