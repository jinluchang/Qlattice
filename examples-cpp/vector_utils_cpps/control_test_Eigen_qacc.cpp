#include <sys/sysinfo.h>
#include <unistd.h>

#include <qlat/qcd.h>

#define Cy qlat::Complex
///#define Et Eigen::RowMajor
#define Et Eigen::ColMajor
#define D 3

__global__ void eigen_test(Cy* va, Cy* vb, Cy* vc, const Long m)
{
  unsigned long index =  threadIdx.y*blockDim.x + threadIdx.x;
  if(index < m){
    //Cy buf[9];
    //for(int ci=0;ci<9;ci++){buf[ci] = va[index*9 +  ci];}
    //Eigen::Matrix<Cy, 3   , 3, Et>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Et>*) &va[index*9 + 0]);
    Eigen::Map< const Eigen::Matrix<Cy, 3   , 3, Et>> lE(&va[index*9 + 0]);

    Cy* bP = &vb[index*3*D];
    //Eigen::Matrix<Cy, D, 3, Et>&     bE = *((Eigen::Matrix<Cy, D, 3, Et>*) bP);
    Eigen::Map< const Eigen::Matrix<Cy, D   , 3, Et>> bE(bP);
    Cy* cP = &vc[index*3*D];
    //Eigen::Matrix<Cy, D, 3, Et>&     cE = *((Eigen::Matrix<Cy, D, 3, Et>*) cP);
    Eigen::Map< Eigen::Matrix<Cy, D   , 3, Et>> cE(cP);

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

  const Long MAX = 100;

  double ini = qlat::u_rand_gen(rs);
  //qlat::vector<Cy > qva;qva.resize((MAX)*9);
  //qlat::vector<Cy > qvb;qvb.resize((MAX)*D*3);
  //qlat::vector<Cy > qvc;qvc.resize((MAX)*D*3);
  //qlat::vector<Cy > qvd;qvd.resize((MAX)*D*3);
  //qlat::vector<Cy > qve;qve.resize((MAX)*D*3);
  //Cy* va = (Cy*) qlat::get_data(qva).data();
  //Cy* vb = (Cy*) qlat::get_data(qvb).data();
  //Cy* vc = (Cy*) qlat::get_data(qvc).data();
  //Cy* vd = (Cy*) qlat::get_data(qvd).data();
  //Cy* ve = (Cy*) qlat::get_data(qve).data();

  Cy* va = NULL;qlat_GPU_Malloc(&va, MAX*3*3*sizeof(Cy));
  Cy* vb = NULL;qlat_GPU_Malloc(&vb, MAX*3*D*sizeof(Cy));
  Cy* vc = NULL;qlat_GPU_Malloc(&vc, MAX*3*D*sizeof(Cy));
  Cy* vd = NULL;qlat_GPU_Malloc(&vd, MAX*3*D*sizeof(Cy));
  Cy* ve = NULL;qlat_GPU_Malloc(&ve, MAX*3*D*sizeof(Cy));

  Cy* vaH = NULL;vaH = (Cy *)malloc(MAX*3*3*sizeof(Cy));
  Cy* vbH = NULL;vbH = (Cy *)malloc(MAX*3*D*sizeof(Cy));
  Cy* vcH = NULL;vcH = (Cy *)malloc(MAX*3*D*sizeof(Cy));
  Cy* vdH = NULL;vdH = (Cy *)malloc(MAX*3*D*sizeof(Cy));
  Cy* veH = NULL;veH = (Cy *)malloc(MAX*3*D*sizeof(Cy));

  //qlat::vector<Cy > va;va.resize((MAX)*9);
  //qlat::vector<Cy > vb;vb.resize((MAX)*D*3);
  //qlat::vector<Cy > vc;vc.resize((MAX)*D*3);

  qacc_for(isp, MAX, {
    for(int ic=0;ic<9;ic++){va[isp*9+ic] = Cy(std::cos((ini+isp + ic)*0.5) , (5.0/(isp + ic+1))*ini*0.1);}
    for(int ic=0;ic<3*D;ic++){
      vb[isp*3*D+ic] = Cy(std::cos((ini+isp + ic)*0.8) , (8.0/(isp + ic+1))*ini*0.2);
      vc[isp*3*D+ic] = vb[isp*3*D+ic];
      ve[isp*3*D+ic] = vb[isp*3*D+ic];
    }
  })
  qlat_GPU_Memcpy(vaH, va, MAX*3*3*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);
  qlat_GPU_Memcpy(vbH, vb, MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);
  qlat_GPU_Memcpy(vcH, vc, MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);
  qlat_GPU_DeviceSynchronize();

  const Long nc =  (MAX + 32-1)/32;
  dim3 dimGrid( nc, 1, 1); 
  dim3 dimBlock(32, 1, 1); 
  eigen_test<<< dimGrid, dimBlock >>>(va, vb, vd, MAX);
  qlat_GPU_DeviceSynchronize();

  //for(Long index=0;index<MAX;index++ ){
  //  Cy buf[9];
  //  for(int ci=0;ci<9;ci++){buf[ci] = vaH[index*9 +  ci];}
  //  Eigen::Matrix<Cy, 3   , 3, Et>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Et>*) buf);

  //  Cy* bP = &vbH[index*3*D];
  //  Eigen::Matrix<Cy, D, 3, Et>&     bE = *((Eigen::Matrix<Cy, D, 3, Et>*) bP);
  //  Cy* eP = &veH[index*3*D];
  //  Eigen::Matrix<Cy, D, 3, Et>&     eE = *((Eigen::Matrix<Cy, D, 3, Et>*) eP);
  //  bE *= lE;
  //  //eE = bE * lE;
  //  //bE = eE;
  //  //bE *= lE; 
  //}
  //qlat_GPU_DeviceSynchronize();

  //qthread_for(index, MAX, {
  //  Cy buf[9];
  //  for(int ci=0;ci<9;ci++){buf[ci] = vaH[index*9 +  ci];}
  //  Eigen::Matrix<Cy, 3   , 3, Et>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Et>*) buf);

  //  Cy* bP = &vbH[index*3*D];
  //  Eigen::Matrix<Cy, D, 3, Et>&     bE = *((Eigen::Matrix<Cy, D, 3, Et>*) bP);
  //  Cy* eP = &veH[index*3*D];
  //  Eigen::Matrix<Cy, D, 3, Et>&     eE = *((Eigen::Matrix<Cy, D, 3, Et>*) eP);
  //  eE = bE * lE;
  //  bE = eE;
  //  //bE *= lE; 
  //}); 
  //qlat_GPU_Memcpy(vb, vbH, MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyHostToDevice);

  qacc_for(index, MAX, {
    //Cy buf[9];
    //for(int ci=0;ci<9;ci++){buf[ci] = va[index*9 +  ci];}
    Eigen::Matrix<Cy, 3   , 3, Et>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Et>*) &va[index*9 + 0]);

    Cy* bP = &vb[index*3*D];
    Eigen::Matrix<Cy, D, 3, Et>&     bE = *((Eigen::Matrix<Cy, D, 3, Et>*) bP);
    Cy* eP = &ve[index*3*D];
    Eigen::Matrix<Cy, D, 3, Et>&     eE = *((Eigen::Matrix<Cy, D, 3, Et>*) eP);
    eE = bE * lE;
    bE = eE;
    //bE *= lE; 
  }); 


  for(Long index=0;index<MAX;index++ )
  {
    //Cy buf[9];
    //for(int ci=0;ci<9;ci++){buf[ci] = vaH[index*9 +  ci];}
    Eigen::Matrix<Cy, 3   , 3, Et>&   lE = *((Eigen::Matrix<Cy, 3   , 3, Et>*) &vaH[index*9 + 0]);

    Cy* bP = &vcH[index*3*D];
    Eigen::Matrix<Cy, D, 3, Et>&     bE = *((Eigen::Matrix<Cy, D, 3, Et>*) bP);
    bE *= lE; 
  }

  qlat_GPU_Memcpy(vbH, vb , MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);
  qlat_GPU_Memcpy(vdH, vd , MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);
  qlat_GPU_Memcpy(veH, ve , MAX*3*D*sizeof(Cy), qlat_GPU_MemcpyDeviceToHost);

  std::vector<double >  diff(3);for(int i=0;i<diff.size();i++){diff[i] = 0;}
  for(Long index=0;index<MAX*3*D;index++ )
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
  qlat_GPU_Free(va);
  qlat_GPU_Free(vb);
  qlat_GPU_Free(vc);
  qlat_GPU_Free(vd);
  qlat_GPU_Free(ve);
  free(vaH);
  free(vbH);
  free(vcH);
  free(vdH);
  free(veH);

  return 0;
}

