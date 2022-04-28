#include <qlat/qlat.h>

#include <complex>
#include <iostream>

using namespace std;
using namespace qlat;

void test1()
{
  TIMER("test1");
  RngState rs(0);
  const int limit = 1024 * 1024;
  double sum = 0.0;
  for (int i = 0; i < limit; ++i) {
    double x = u_rand_gen(rs);
    double y = std::sin(PI / x);
    sum += y;
  }
  // Integrate[Sin[Pi/x],{x,0,1}] = -Pi Ci(Pi) = -0.231435
  cout << "expected: " << -0.231435 << endl;
  cout << sum / limit << endl;
}

void test2()
{
  TIMER("test2");
  RngState rs(0);
  const int limit = 1024 * 1024;
  int count = 0;
  for (int i = 0; i < limit; ++i) {
    double x = g_rand_gen(rs);
    if (std::abs(x) <= 1.0) {
      count += 1;
    }
  }
  cout << "expected: " << 0.682689492137 << endl;
  cout << (double)count / limit << endl;
}

void test2a()
{
  TIMER("test2a");
  RngState rs(0);
  const int limit = 1024 * 1024;
  int count = 0;
  for (int i = 0; i < limit; ++i) {
    double x = g_rand_gen(rs, 1.0, 2.0);
    if (std::abs(x - 1.0) <= 2.0) {
      count += 1;
    }
  }
  cout << "expected: " << 0.682689492137 << endl;
  cout << (double)count / limit << endl;
}

inline double sqr(double x) { return x * x; }

void test3()
{
  TIMER("test3");
  RngState rs(123123);
  const int Nb = 128;
  const int Ni = 8;
  const int Ndrop = 1024;
  const int Ntake = 8;
  double sum = 0;
  double sigma2 = 0;
  for (int block = 0; block < Nb; block++) {
    Complex a = 0;
    for (int id = 0; id < Ni; id++) {
      int index = block * Ni + id;
      reset(rs, index);
      for (int i = 0; i < Ndrop; i++) {
        u_rand_gen(rs);
      }
      for (int i = 0; i < Ntake; i++) {
        a += polar(1.0, u_rand_gen(rs, PI, -PI));
      }
    }
    sum += qnorm(a);
    sigma2 += sqr(qnorm(a));
  }
  cout << "Expected : " << Ni * Ntake << endl;
  cout << "Mean     : " << sum / Nb << endl;
  cout << "Var      : " << sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt(Nb - 1)
       << endl;
}

void test4()
{
  TIMER("test4");
  RngState rs(get_global_rng_state(), "test4");
  const int Nb = 128;
  const int Ni = 8;
  const int Ndrop = 1024;
  const int Ntake = 8;
  double sum = 0;
  double sigma2 = 0;
  for (int block = 0; block < Nb; block++) {
    Complex a = 0;
    for (int id = 0; id < Ni; id++) {
      int index = block * Ni + id;
      RngState rsi = rs.newtype(index);
      for (int i = 0; i < Ndrop; i++) {
        u_rand_gen(rsi);
      }
      for (int i = 0; i < Ntake; i++) {
        a += polar(1.0, u_rand_gen(rsi, PI, -PI));
      }
    }
    sum += qnorm(a);
    sigma2 += sqr(qnorm(a));
  }
  cout << "Expected : " << Ni * Ntake << endl;
  cout << "Mean     : " << sum / Nb << endl;
  cout << "Var      : " << sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt(Nb - 1)
       << endl;
}

void profile()
{
  TIMER_FLOPS("profile");
  RngState rs(0);
  const int size = 1000 * 1000;
  double sum = 0.0;
  for (int i = 0; i < size; ++i) {
    sum += u_rand_gen(rs);
  }
  cout << "sum: " << sum << endl;
  cout << "avg: " << sum / size << endl;
  timer.flops += size;
}

void profileOmp()
{
  TIMER_FLOPS("profileOmp");
  const int num_threads = 64;
  const int size = 16 * 1024;
  double sums[num_threads];
#pragma omp parallel for
  for (int i = 0; i < num_threads; ++i) {
    RngState rs(i);
    double psum = 0.0;
    for (int k = 0; k < size; ++k) {
      psum += u_rand_gen(rs);
    }
    sums[i] = psum;
  }
  double sum = 0.0;
  for (int i = 0; i < num_threads; ++i) {
    sum += sums[i];
  }
  cout << "sum: " << sum << endl;
  cout << "avg: " << sum / (num_threads * size) << endl;
  timer.flops += num_threads * size;
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  reset(get_global_rng_state(), "rng-state-tests");
  if (get_id_node() == 0) {
    cout << "sizeof RngState = " << sizeof(RngState) << endl;
    test1();
    test2();
    test2a();
    test3();
    test4();
    profileOmp();
    profile();
  }
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
