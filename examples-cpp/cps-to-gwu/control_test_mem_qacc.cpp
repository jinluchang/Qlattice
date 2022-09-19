#include <sys/sysinfo.h>
#include <unistd.h>

#include <qlat/qcd.h>

#define Cy qlat::Complex

int main(int argc, char* argv[])
{
  using namespace qlat;
  //MPI_Init(&argc, &argv);
  //int id_node;
  //MPI_Comm_rank(MPI_COMM_WORLD, &id_node);

  const int seed = 131;
  qlat::RngState rs(qlat::get_id_node() + 1 + seed);

  const long MAX = 100;

  double ini = qlat::u_rand_gen(rs);
  {
  qlat::vector_acc<Cy > va;va.resize((MAX)*9);

  qacc_for(isp, MAX, {
    for(int ic=0;ic<9;ic++){va[isp*9+ic] = Cy(std::cos((ini+isp + ic)*0.5) , (5.0/(isp + ic+1))*ini*0.1);}
  })
  //va.resize(0);
  }

  qlat::end();
  printf("aaaaaaaa");

  return 0;
}

