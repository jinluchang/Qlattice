#include <qlat/qcd.h>
#include <sys/sysinfo.h>



int main(int argc, char* argv[])
{
  using namespace qlat;

  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1,  1));
  size_node_list.push_back(Coordinate(1, 1, 1,  2));
  size_node_list.push_back(Coordinate(1, 1, 1,  4));
  size_node_list.push_back(Coordinate(1, 1, 1,  8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(1, 1, 1, 64));
  size_node_list.push_back(Coordinate(4, 4, 8, 16));
  size_node_list.push_back(Coordinate(4, 8, 8, 16));

  begin(&argc, &argv, size_node_list);
  {
  #ifdef QLAT_USE_ACC
  int num_node;MPI_Comm_size(get_comm(), &num_node);
  int id_node;MPI_Comm_rank(get_comm(), &id_node);

  int num_gpus = 0; 
  qlat_GPU_GetDeviceCount(&num_gpus);
  ////qlat_GPU_DeviceReset();
  qlat_GPU_SetDevice(id_node % num_gpus);
  int gpu_id = -1;  
  qlat_GPU_GetDevice(&gpu_id);
  //printf("CPU node %d (of %d) uses CUDA device %d\n", id_node, num_node, gpu_id);
  fflush(stdout);
  MPI_Barrier(get_comm());
  #endif
  }

  Long na = 64;

  qlat::vector<ComplexD > src;src.resize(na);
  qlat::vector<ComplexD > s0;s0.resize(na);
  qlat::vector<ComplexD > s1;s1.resize(na);
  ComplexD* k = (ComplexD*) qlat::get_data(src).data();
  qacc_barrier(dummy);
  qacc_for(isp , na, {
    //src[isp] = src[isp] + 1;
    k[isp] = ComplexD(1.0, 0.0);
    s1[isp] = isp;
    s0[isp] = s1[isp] * s1[isp];
  });
  qacc_barrier(dummy);
  for(int isp=0;isp<src.size();isp++)
  {
    qlat::displayln_info(qlat::ssprintf("t %5d %.8e ", isp, s0[isp].real()));
  }

  qlat::Timer::display();
  qlat::end();
  return 0;
}

