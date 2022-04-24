#include <sys/sysinfo.h>
#include <unistd.h>

#include <qlat/qcd.h>

#define Cy qlat::Complex

__global__ void eigen_test(Cy* va, Cy* vb, Cy* vc, const long m)
{
  unsigned long index =  threadIdx.y*blockDim.x + threadIdx.x;
  if(index < m){
    Cy buf[9];
    for(int ci=0;ci<9;ci++){buf[ci] = va[index*9 +  ci];}
    Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>*) buf);

    Cy* bP = &vb[index*3*2];
    Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>&     bE = *((Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>*) bP);
    Cy* cP = &vc[index*3*2];
    Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>&     cE = *((Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>*) cP);

    ////bE = bE * lE; 
    cE = bE * lE;
    //////bE = cE;
  }
}

int main(int argc, char* argv[])
{
  using namespace qlat;
  MPI_Init(&argc, &argv);
  int id_node;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_node);

  const int seed = 131;
  qlat::RngState rs(qlat::get_id_node() + 1 + seed);

  const long MAX = 100;

  double ini = qlat::u_rand_gen(rs);
  //qlat::vector_acc<Cy > qva;qva.resize((MAX)*9);
  //qlat::vector_acc<Cy > qvb;qvb.resize((MAX)*2*3);
  //qlat::vector_acc<Cy > qvc;qvc.resize((MAX)*2*3);
  //qlat::vector_acc<Cy > qvd;qvd.resize((MAX)*2*3);
  //qlat::vector_acc<Cy > qve;qve.resize((MAX)*2*3);
  //Cy* va = (Cy*) qlat::get_data(qva).data();
  //Cy* vb = (Cy*) qlat::get_data(qvb).data();
  //Cy* vc = (Cy*) qlat::get_data(qvc).data();
  //Cy* vd = (Cy*) qlat::get_data(qvd).data();
  //Cy* ve = (Cy*) qlat::get_data(qve).data();

  Cy* va = NULL;cudaMalloc(&va, MAX*3*3*sizeof(Cy));
  Cy* vb = NULL;cudaMalloc(&vb, MAX*3*2*sizeof(Cy));
  Cy* vc = NULL;cudaMalloc(&vc, MAX*3*2*sizeof(Cy));
  Cy* vd = NULL;cudaMalloc(&vd, MAX*3*2*sizeof(Cy));
  Cy* ve = NULL;cudaMalloc(&ve, MAX*3*2*sizeof(Cy));

  Cy* vaH = NULL;vaH = (Cy *)malloc(MAX*3*3*sizeof(Cy));
  Cy* vbH = NULL;vbH = (Cy *)malloc(MAX*3*2*sizeof(Cy));
  Cy* vcH = NULL;vcH = (Cy *)malloc(MAX*3*2*sizeof(Cy));
  Cy* vdH = NULL;vdH = (Cy *)malloc(MAX*3*2*sizeof(Cy));
  Cy* veH = NULL;veH = (Cy *)malloc(MAX*3*2*sizeof(Cy));

  //qlat::vector_acc<Cy > va;va.resize((MAX)*9);
  //qlat::vector_acc<Cy > vb;vb.resize((MAX)*2*3);
  //qlat::vector_acc<Cy > vc;vc.resize((MAX)*2*3);

  qacc_for(isp, MAX, {
    for(int ic=0;ic<9;ic++){va[isp*9+ic] = Cy(std::cos((ini+isp + ic)*0.5) , (5.0/(isp + ic+1))*ini*0.1);}
    for(int ic=0;ic<6;ic++){
      vb[isp*6+ic] = Cy(std::cos((ini+isp + ic)*0.8) , (8.0/(isp + ic+1))*ini*0.2);
      //vc[isp*6+ic] = vb[isp*6+ic];
      ve[isp*6+ic] = vb[isp*6+ic];
    }
  })
  cudaMemcpy(vaH, va, MAX*3*3*sizeof(Cy), cudaMemcpyDeviceToHost);
  cudaMemcpy(vbH, vb, MAX*3*2*sizeof(Cy), cudaMemcpyDeviceToHost);
  cudaMemcpy(vcH, vc, MAX*3*2*sizeof(Cy), cudaMemcpyDeviceToHost);

  const long nc =  (MAX + 32-1)/32;
  dim3 dimGrid( nc, 1, 1); 
  dim3 dimBlock(32, 1, 1); 
  eigen_test<<< dimGrid, dimBlock >>>(va, vb, vd, MAX);
  cudaDeviceSynchronize();

  qacc_for(index, MAX, {
    Cy buf[9];
    for(int ci=0;ci<9;ci++){buf[ci] = va[index*9 +  ci];}
    Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>*) buf);

    Cy* bP = &vb[index*3*2];
    Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>&     bE = *((Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>*) bP);
    Cy* eP = &vb[index*3*2];
    Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>&     eE = *((Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>*) eP);
    eE = bE * lE;
    bE = eE;
    //bE *= lE; 
  }); 

  for(long index=0;index<MAX;index++ )
  {
    Cy buf[9];
    for(int ci=0;ci<9;ci++){buf[ci] = vaH[index*9 +  ci];}
    Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Eigen::RowMajor>*) buf);

    Cy* bP = &vcH[index*3*2];
    Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>&     bE = *((Eigen::Matrix<Cy, 2, 3, Eigen::ColMajor>*) bP);
    bE *= lE; 
  }

  cudaMemcpy(vbH, vb , MAX*3*2*sizeof(Cy), cudaMemcpyDeviceToHost);
  cudaMemcpy(vdH, vd , MAX*3*2*sizeof(Cy), cudaMemcpyDeviceToHost);
  cudaMemcpy(veH, ve , MAX*3*2*sizeof(Cy), cudaMemcpyDeviceToHost);

  std::vector<double >  diff(3);for(int i=0;i<diff.size();i++){diff[i] = 0;}
  for(long index=0;index<MAX*6;index++ )
  {
    diff[0] += qlat::qnorm(vcH[index] - vbH[index]);
    diff[1] += qlat::qnorm(vcH[index] - vdH[index]);
    diff[2] += qlat::qnorm(veH[index]);
  }

  printf("rank %d, diff %.3e, %.3e %.3e \n", id_node, diff[0], diff[1], diff[2]);
  MPI_Barrier(MPI_COMM_WORLD);
  fflush(stdout);

  //qva.resize(0);
  //qvb.resize(0);
  //qvc.resize(0);
  //qvd.resize(0);
  //qve.resize(0);
  cudaFree(va);
  cudaFree(vb);
  cudaFree(vc);
  cudaFree(vd);
  cudaFree(ve);
  free(vaH);
  free(vbH);
  free(vcH);
  free(vdH);
  free(veH);

  return 0;
}

