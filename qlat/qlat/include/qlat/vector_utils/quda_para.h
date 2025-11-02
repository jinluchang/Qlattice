#ifndef QUDA_PARA_H
#define QUDA_PARA_H

#pragma once

#include <quda.h>
#include <invert_quda.h>

#include <cstdlib>
#include "utils_float_type.h"
#include "utils_stagger_contractions.h"
#include "utils_gauge_field.h"

namespace quda
{

//quda::Complex& operator+=(const thrust::complex<double>& rhs){
//  double a = real + rhs.real();
//  double b = imag + rhs.imag();
//  *this = quda::Complex(a, b);
//  //real() += rhs.real();
//  //imag() += rhs.imag();
//  return *this;
//}
//quda::Complex& operator+=(const thrust::complex<float>& rhs){
//  real() += rhs.real();
//  imag() += rhs.imag();
//  return *this;
//}
//inline qlat::ComplexT<double> operator+(const quda::Complex &a, const qlat::ComplexT<double> &b) {
//    return qlat::ComplexT<double>(a.real() + b.real(), a.imag() + b.imag());
//}
//inline quda::Complex operator+(const quda::Complex &a, const qlat::ComplexT<double> &b) {
//    return quda::Complex<double>(a.real() + b.real(), a.imag() + b.imag());
//}
//inline qlat::ComplexT<double> operator*(const quda::Complex &a, const qlat::ComplexT<double> &b) {
//    return qlat::ComplexT<double>(a.real() + b.real(), a.imag() + b.imag());
//}


void massRescale(ColorSpinorField &b, QudaInvertParam &param, bool for_multishift)
{
  double kappa5 = (0.5/(5.0 + param.m5));
  double kappa = (param.dslash_type == QUDA_DOMAIN_WALL_DSLASH || param.dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH
                  || param.dslash_type == QUDA_MOBIUS_DWF_DSLASH || param.dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH) ?
    kappa5 :
    param.kappa;

  logQuda(QUDA_DEBUG_VERBOSE, "Mass rescale: Kappa is: %g\n", kappa);
  logQuda(QUDA_DEBUG_VERBOSE, "Mass rescale: mass normalization: %d\n", param.mass_normalization);
  logQuda(QUDA_DEBUG_VERBOSE, "Mass rescale: norm of source in = %g\n", blas::norm2(b));

  // staggered dslash uses mass normalization internally
  if (param.dslash_type == QUDA_ASQTAD_DSLASH || param.dslash_type == QUDA_STAGGERED_DSLASH) {
    switch (param.solution_type) {
      case QUDA_MAT_SOLUTION:
      case QUDA_MATPC_SOLUTION:
        if (param.mass_normalization == QUDA_KAPPA_NORMALIZATION) blas::ax(2.0*param.mass, b);
        break;
      case QUDA_MATDAG_MAT_SOLUTION:
      case QUDA_MATPCDAG_MATPC_SOLUTION:
        if (param.mass_normalization == QUDA_KAPPA_NORMALIZATION) blas::ax(4.0*param.mass*param.mass, b);
        break;
      default:
        errorQuda("Not implemented");
    }
    return;
  }

  // multiply the source to compensate for normalization of the Dirac operator, if necessary
  // you are responsible for restoring what's in param.offset
  switch (param.solution_type) {
    case QUDA_MAT_SOLUTION:
      if (param.mass_normalization == QUDA_MASS_NORMALIZATION ||
          param.mass_normalization == QUDA_ASYMMETRIC_MASS_NORMALIZATION) {
        blas::ax(2.0*kappa, b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 2.0 * kappa;
      }
      break;
    case QUDA_MATDAG_MAT_SOLUTION:
      if (param.mass_normalization == QUDA_MASS_NORMALIZATION ||
          param.mass_normalization == QUDA_ASYMMETRIC_MASS_NORMALIZATION) {
        blas::ax(4.0*kappa*kappa, b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 4.0 * kappa * kappa;
      }
      break;
    case QUDA_MATPC_SOLUTION:
      if (param.mass_normalization == QUDA_MASS_NORMALIZATION) {
        blas::ax(4.0*kappa*kappa, b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 4.0 * kappa * kappa;
      } else if (param.mass_normalization == QUDA_ASYMMETRIC_MASS_NORMALIZATION) {
        blas::ax(2.0*kappa, b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 2.0 * kappa;
      }
      break;
    case QUDA_MATPCDAG_MATPC_SOLUTION:
      if (param.mass_normalization == QUDA_MASS_NORMALIZATION) {
        blas::ax(16.0*std::pow(kappa,4), b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 16.0 * std::pow(kappa, 4);
      } else if (param.mass_normalization == QUDA_ASYMMETRIC_MASS_NORMALIZATION) {
        blas::ax(4.0*kappa*kappa, b);
        if (for_multishift)
          for (Int i = 0; i < param.num_offset; i++) param.offset[i] *= 4.0 * kappa * kappa;
      }
      break;
    default:
      errorQuda("Solution type %d not supported", param.solution_type);
  }

  logQuda(QUDA_DEBUG_VERBOSE, "Mass rescale: norm of source out = %g\n", blas::norm2(b));
}

}

namespace qlat
{  //

void begin_Qlat_with_quda(bool t = false)
{
  qlat::Coordinate nodeD;int dims[4];int coords[4];
  for(Int d=0;d<4;d++){
    nodeD[d] = quda::comm_dim(d);dims[d] = quda::comm_dim(d);
    coords[d] = quda::comm_coord(d);
  }

  Int rank = 0;
  if(t == true ){
    rank = coords[0];
    for (Int i = 1; i <= 3; i++) {
      rank = dims[i] * rank + coords[i];
    }
  }
  if(t == false){
    rank = coords[3];
    for (Int i = 2; i >= 0; i--) {
      rank = dims[i] * rank + coords[i];
    }
  }

  ////printf("rank %d, nodeD, %d %d %d %d \n", id_node, nodeD[0], nodeD[1], nodeD[2], nodeD[3]);
  qlat::begin(rank, nodeD);
  //printf("rank %d check %d, ", quda::comm_rank(), qlat::get_id_node());
  /////Qassert(quda::comm_rank() == qlat::get_id_node());

  for (Int d = 0; d < 4; d++) {
    Qassert(quda::comm_coord(d) == qlat::get_coor_node()[d]);
    //printf("%d  %d, ", quda::comm_coord(d),  qlat::get_coor_node()[d]);
  }
  //printf("\n");
}

static Int mpi_rank_from_coords_x(const Int* coords, void* fdata)
{
  Int* dims = reinterpret_cast<int*>(fdata);
  //
  Int rank;
  rank = coords[3];
  for (Int i = 2; i >= 0; i--) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
static Int mpi_rank_from_coords_t(const Int* coords, void* fdata)
{
  Int* dims = reinterpret_cast<int*>(fdata);
  //
  Int rank;
  rank = coords[0];
  for (Int i = 1; i <= 3; i++) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
inline void quda_begin_internal(Int mpi_layout[4], Int quda_rankx = 1)
{
  using namespace quda;
  // The following sets the MPI comm stuff.
  MPI_Comm comm = get_comm();
  //qudaSetCommHandle((void*) &comm);
  setMPICommHandleQuda((void*) &comm);
  ////default t = false, x running the fast

  ///int t = 0;
  //Coordinate node =  get_size_node();
  //Coordinate cor  = qlat::get_coor_node();
  //Coordinate max  =  Coordinate(0,0,0, node[3]-1);
  //if(cor == max)
  //{ 
  //  if(qlat::get_id_node() == cor[3]){
  //    t = 1;
  //  }
  //}
  //sum_all_size(&t, 1, 0);

  if(quda_rankx >= 1)
  {
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_x,
                      reinterpret_cast<void*>(mpi_layout));
  }else{
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_t,
                      reinterpret_cast<void*>(mpi_layout));
  }

  // comm_set_gridsize(mpi_layout);
  //int gpu_id = -1;qacc_ErrCheck(qacc_GetDevice(&gpu_id));
  //initQuda(gpu_id);
  initQuda(-1000);
  //initQuda(-1);

  //printf(
  //    "initialized on quda rank #%03d (%03d,%03d,%03d,%03d), qlat rank #%03d "
  //    "(%03d,%03d,%03d,%03d).\n",
  //    comm_rank(), comm_coord(0), comm_coord(1), comm_coord(2), comm_coord(3),
  //    get_id_node(), get_coor_node()[0], get_coor_node()[1], get_coor_node()[2],
  //    get_coor_node()[3]);
  // Make sure there is no mismatch

  Qassert(comm_rank() == get_id_node());
  for (Int d = 0; d < 4; d++) {
    Qassert(comm_coord(d) == get_coor_node()[d]);
  }
  qacc_ErrCheck(qacc_DeviceSetCacheConfig(qacc_FuncCachePreferNone ));// for smearing settings
}

/*
  need to begin quda before any GPU memory allocation
  Assume Qlattice does not assign GPUs
*/
inline void begin_quda_with_gpu()
{
  using namespace quda;

  Int mpi_layout[4]={0,0,0,0};
  qlat::GeometryNode geon = qlat::get_geometry_node();for(Int i=0;i<4;i++){mpi_layout[i] = geon.size_node[i];}

  Int rank = qlat::get_id_node();
  qlat::Coordinate coords  = qlat::get_coor_node();
  Int rankx = coords[3];
  for (Int i = 2; i >= 0; i--) { 
    rankx = mpi_layout[i] * rankx + coords[i];
  }
  
  Int rankt = coords[0];
  for (Int i = 1; i <= 3; i++) { 
    rankt = mpi_layout[i] * rankt + coords[i];
  }
  
  Int quda_rankx = 0;
  if(rankt != rankx){
    //if(rank == rankt){
    //  printf("T rank! %3d %3d \n", rank, rankt);               
    //}
    if(rank == rankx){
      //printf("X rank! %3d %3d \n", rank, rankx);
      quda_rankx = 1;
    }
    if(rank != rankx and rank != rankt){printf("UNKNOWN rank! \n");}
  }
  qlat::sum_all_size(&quda_rankx, 1);
  //qmessage("Rank X %d \n", quda_rankx);

  qlat::quda_begin_internal(mpi_layout, quda_rankx);

}

inline void begin_Qlat_quda(int* argc, char** argv[], inputpara& in, Int read_Lat = 0)
{
  begin_Lat(argc, argv, in, 0, read_Lat);
  begin_quda_with_gpu();
}

inline void check_quda_layout_eo(const Geometry& geo)
{
  ////qlat::GeometryNode geon = qlat::get_geometry_node();for(Int i=0;i<4;i++){mpi_layout[i] = geon.size_node[i];}
  ////Nv.resize(4);nv.resize(4);mv.resize(4);
  ////for(Int i=0;i<4;i++){Nv[i]=geo.node_site[i];nv[i] = geo.node_site[i] * geo.geon.size_node[i];}
  ////for(Int i=0;i<4;i++){mv[i] = nv[i]/Nv[i];}
  for(Int i=0;i<4;i++)
  {
    if(geo.node_site[i] % 2 != 0)
    {
      qmessage("dir %5d, size %5d wrong! \n", i, geo.node_site[i]);
      abort_r();
    }
  }
}

inline void quda_end()
{
  using namespace quda;
  endQuda();
}

template <class Ta, class Ty>
void quda_convert_gauge(Ta* qgf, Ty* gf, const Geometry& geo, Int dir = 0, Int GPU = 0)
{
  TIMER("quda_convert_gauge(qgf,gf)");
  //Qassert(geo.multiplicity == 1);

  //ColorMatrix* quda_pt = reinterpret_cast<ColorMatrix*>(qgf.data());
  const Long V = geo.local_volume();
  const Long Vh = V / 2;
  qGPU_for(qlat_idx, V, GPU, {
    Coordinate xl = geo.coordinate_from_index(qlat_idx);
    //const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
    //Vector<ColorMatrix> ms = gf.get_elems(xl);
    Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    Long quda_idx = (qlat_idx / 2 + eo * Vh) * 4 + 0;
    Ty* ms  =  &gf[su3_n(geo, xl, 0)];
    Ta*  q  = &qgf[quda_idx*9 + 0];
    if(dir == 0){for(Int i=0;i< 4 * 9;i++){q[i]  = ms[i];}}
    if(dir == 1){for(Int i=0;i< 4 * 9;i++){ms[i] =  q[i];}}
  });
}

template <class Ta, class Td>
void quda_convert_gauge(qlat::vector<qlat::ComplexT<Ta > >& qgf, GaugeFieldT<Td >& gf, Int dir = 0, Int GPU = 0)
{
  TIMER("quda_convert_gauge(qgf,gf)");
  const Geometry& geo = gf.geo();
  //Qassert(geo.multiplicity == 4);
  //geo.multiplicity = 1;
  qlat::ComplexT<Ta >* res = (qlat::ComplexT<Ta >*) qgf.data();
  qlat::ComplexT<Td >* src = (qlat::ComplexT<Td >*) qlat::get_data(gf).data();
  quda_convert_gauge(res, src, geo, dir, GPU);

  //qlat::ComplexT<T >* res = (qlat::ComplexT<T >*) qgf.data();
  //qlat::ComplexT<T >* src = (qlat::ComplexT<T >*) qlat::get_data(gf).data();
  //const Long V = geo.local_volume();
  //const Long Vh = V / 2;

  //qGPU_for(qlat_idx, V, GPU, {
  //  const Coordinate xl = geo.coordinate_from_index(qlat_idx);
  //  const Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
  //  const Long quda_idx = qlat_idx / 2 + eo * Vh;
  //  qlat::ComplexT<T >* s0 = &src[(qlat_idx*4 + 0) * 9 ];
  //  qlat::ComplexT<T >* s1 = &res[(quda_idx*4 + 0) * 9 ];
  //  if(dir == 0)for(Int i=0;i<4*9;i++){s1[i] = s0[i];}
  //  if(dir == 1)for(Int i=0;i<4*9;i++){s0[i] = s1[i];}
  //});
  //ColorMatrix* quda_pt = reinterpret_cast<ColorMatrix*>(qgf.data());
  ////Qassert(geo.multiplicity == 4);
  //const Long V = geo.local_volume();
  //const Long Vh = V / 2;
  //qGPU_for(qlat_idx, V, GPU, {
  //  Coordinate xl = geo.coordinate_from_index(qlat_idx);
  //  //const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
  //  //Vector<ColorMatrix> ms = gf.get_elems(xl);
  //  Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
  //  for (Int mu = 0; mu < 4; mu++) {
  //    //ColorMatrixT<T>& ms = gf.get_elem_offset(qlat_idx*gf.geo().multiplicity+mu);
  //    //ColorMatrixT<Td >& ms = Td* gf.get_elem(xl, mu).p;
  //    Long quda_idx = (qlat_idx / 2 + eo * Vh) * 4 + mu;
  //    Td* ms  = (Td*) gf.get_elem(xl, mu).p;
  //    Ta*  q  = (Ta*)  &quda_pt[quda_idx];
  //    if(dir == 0){for(Int i=0;i<18;i++){q[i]  = ms[i];}}
  //    if(dir == 1){for(Int i=0;i<18;i++){ms[i] = q[i]; }}
  //    //if(dir == 0){quda_pt[quda_idx] = ms;}
  //    //if(dir == 1){ms = quda_pt[quda_idx];}
  //  }
  //});
}

QudaPrecision get_quda_precision(Int byte)
{
  switch (byte) {
    case 8:
      return QUDA_DOUBLE_PRECISION;
      break;
    case 4:
      return QUDA_SINGLE_PRECISION; 
      break;
    case 2:
      return QUDA_HALF_PRECISION; 
      break;
    default:
      Qassert(false);
      return QUDA_INVALID_PRECISION;
  }
}


//template <class T>
//void quda_convert_fermion_copy(FermionField5d& ff, const std::vector<T>& qff)
//{
//  TIMER("quda_convert_fermion(ff,qff)");
//  const Geometry& geo = ff.geo();
//  const WilsonVector* quda_pt =
//      reinterpret_cast<const WilsonVector*>(qff.data());
//  Int Ls = geo.multiplicity;
//  Qassert(Ls > 0);
//  Long V = geo.local_volume();
//  Long Vh = V / 2;
//
//  //
//  #pragma omp parallel for
//  for (Long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
//    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
//    Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
//    Vector<WilsonVector> wvs = ff.get_elems(xl);
//    for (Int s = 0; s < Ls; s++) {
//      Int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
//      wvs[s] = quda_pt[quda_idx];
//    }
//  }
//
//}

//template <class T>
//void quda_convert_fermion_copy(std::vector<T>& qff, const FermionField5d& ff)
//{
//  TIMER("quda_convert_fermion(qff,ff)");
//  const Geometry& geo = ff.geo();
//  WilsonVector* quda_pt = reinterpret_cast<WilsonVector*>(qff.data());
//  Int Ls = geo.multiplicity;
//  Qassert(Ls > 0);
//  Long V = geo.local_volume();
//  Long Vh = V / 2;
////
//#pragma omp parallel for
//  for (Long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
//    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
//    Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
//    const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
//    for (Int s = 0; s < Ls; s++) {
//      Int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
//      quda_pt[quda_idx] = wvs[s];
//    }
//  }
//}


// Helper functions
//------------------------------------------------------
inline bool isPCSolution(QudaSolutionType solution_type)
{
  return (solution_type == QUDA_MATPC_SOLUTION || solution_type == QUDA_MATPC_DAG_SOLUTION
          || solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
}
//------------------------------------------------------

void constructWilsonTestSpinorParam(quda::ColorSpinorParam *cs_param, const QudaInvertParam *inv_param,
                                    const QudaGaugeParam *gauge_param)
{
  // Lattice vector spacetime/colour/spin/parity properties
  cs_param->nColor = 3;
  cs_param->nSpin = 4;
  if (inv_param->dslash_type == QUDA_DOMAIN_WALL_DSLASH || inv_param->dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH
      || inv_param->dslash_type == QUDA_MOBIUS_DWF_DSLASH || inv_param->dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH) {
    cs_param->nDim = 5;
    cs_param->x[4] = inv_param->Ls;
  } else if ((inv_param->dslash_type == QUDA_TWISTED_MASS_DSLASH || inv_param->dslash_type == QUDA_TWISTED_CLOVER_DSLASH)
             && (inv_param->twist_flavor == QUDA_TWIST_NONDEG_DOUBLET)) {
    cs_param->nDim = 5;
    cs_param->x[4] = 2;
  } else {
    cs_param->nDim = 4;
  }
  cs_param->pc_type = inv_param->dslash_type == QUDA_DOMAIN_WALL_DSLASH ? QUDA_5D_PC : QUDA_4D_PC;
  for (Int d = 0; d < 4; d++) cs_param->x[d] = gauge_param->X[d];
  bool pc = isPCSolution(inv_param->solution_type);
  if (pc) cs_param->x[0] /= 2;
  cs_param->siteSubset = pc ? QUDA_PARITY_SITE_SUBSET : QUDA_FULL_SITE_SUBSET;

  // Lattice vector data properties
  cs_param->setPrecision(inv_param->cpu_prec);
  cs_param->pad = 0;
  cs_param->siteOrder = QUDA_EVEN_ODD_SITE_ORDER;
  cs_param->fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
  cs_param->gammaBasis = inv_param->gamma_basis;
  cs_param->create = QUDA_ZERO_FIELD_CREATE;
  cs_param->location = QUDA_CPU_FIELD_LOCATION;
}

//// data reordering routines
//template <typename Out, typename In>
//void reorderQDPtoMILC(Out *milc_out, In *qdp_in, Long V, Int siteSize)
//{
//  qthread_for(i, V,{
//    for (Int dir = 0; dir < 4; dir++) {
//      for (Int j = 0; j < siteSize; j++) {
//        milc_out[(i * 4 + dir) * siteSize + j] = static_cast<Out>(qdp_in[dir][i * siteSize + j]);
//      }
//    }
//  });
//}

bool last_node_in_t()
{
  using namespace quda;
  // only apply T-boundary at edge nodes
  return commCoords(3) == commDim(3) - 1; 
}

Int fullLatticeIndex(Int dim[4], Int index, Int oddBit)
{

  Int za = index / (dim[0] >> 1);
  Int zb = za / dim[1];
  Int x2 = za - zb * dim[1];
  Int x4 = zb / dim[2];
  Int x3 = zb - x4 * dim[2];

  return 2 * index + ((x2 + x3 + x4 + oddBit) & 1);
}

// given a "half index" i into either an even or odd half lattice (corresponding
// to oddBit = {0, 1}), returns the corresponding full lattice index.
// TODO check Z is the dimension
Int fullLatticeIndex(Int i, Int oddBit, Int Z[4])
{
  /*
    Int boundaryCrossings = i/(Z[0]/2) + i/(Z[1]*Z[0]/2) + i/(Z[2]*Z[1]*Z[0]/2);
    return 2*i + (boundaryCrossings + oddBit) % 2;
  */

  Int X1 = Z[0];
  Int X2 = Z[1];
  Int X3 = Z[2];
  // Int X4 = Z[3];
  Int X1h = X1 / 2;

  Int sid = i;
  Int za = sid / X1h;
  // Int x1h = sid - za*X1h;
  Int zb = za / X2;
  Int x2 = za - zb * X2;
  Int x4 = zb / X3;
  Int x3 = zb - x4 * X3;
  Int x1odd = (x2 + x3 + x4 + oddBit) & 1;
  // Int x1 = 2*x1h + x1odd;
  Int X = 2 * sid + x1odd;

  return X;
}


template <typename Ty>
void applyGaugeFieldScaling_long(Ty *gauge, Long Vh, QudaGaugeParam *param, QudaDslashType dslash_type = QUDA_STAGGERED_DSLASH)
{
  TIMER("applyGaugeFieldScaling_long");
  Int X1h = param->X[0] / 2;
  Int X1 = param->X[0];
  Int X2 = param->X[1];
  Int X3 = param->X[2];
  Int X4 = param->X[3];

  const Long V = Vh*2;

  // rescale Long links by the appropriate coefficient
  if (dslash_type == QUDA_ASQTAD_DSLASH) {
    #pragma omp parallel for
    for(Long isp=0;isp < V; isp ++)
    for (Int d = 0; d < 4 * 9; d++) {
      {
        gauge[isp*4*9 + d] /= (-24 * param->tadpole_coeff * param->tadpole_coeff);
      }
    }
  }

  // apply the staggered phases
  for (Int d = 0; d < 3; d++) {

    // even
    #pragma omp parallel for
    for (Int i = 0; i < Vh; i++) {

      Int index = fullLatticeIndex(i, 0, param->X);
      Int i4 = index / (X3 * X2 * X1);
      Int i3 = (index - i4 * (X3 * X2 * X1)) / (X2 * X1);
      Int i2 = (index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1)) / X1;
      Int i1 = index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1) - i2 * X1;
      Int sign = 1;

      if (d == 0) {
        if (i4 % 2 == 1) { sign = -1; }
      }

      if (d == 1) {
        if ((i4 + i1) % 2 == 1) { sign = -1; }
      }
      if (d == 2) {
        if ((i4 + i1 + i2) % 2 == 1) { sign = -1; }
      }

      Long indexe = ((0 * Vh + i)*4 + d)*9 ;
      for (Int j = 0; j < 9; j++) { gauge[indexe + j] *= sign; }
      //if(i1 == 1 and i2 == 1 and i3 == 1 and i4 == 1){qmessage("==even  %d \n", int(sign));}
    }
    // odd
    #pragma omp parallel for
    for (Int i = 0; i < Vh; i++) {
      Int index = fullLatticeIndex(i, 1, param->X);
      Int i4 = index / (X3 * X2 * X1);
      Int i3 = (index - i4 * (X3 * X2 * X1)) / (X2 * X1);
      Int i2 = (index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1)) / X1;
      Int i1 = index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1) - i2 * X1;
      Int sign = 1;

      if (d == 0) {
        if (i4 % 2 == 1) { sign = -1; }
      }

      if (d == 1) {
        if ((i4 + i1) % 2 == 1) { sign = -1; }
      }
      if (d == 2) {
        if ((i4 + i1 + i2) % 2 == 1) { sign = -1; }
      }

      Long indexo = ((1 * Vh + i)*4 + d)*9 ;
      ////for (Int j = 0; j < 18; j++) { gauge[d][(Vh + i) * gauge_site_size + j] *= sign; }
      for (Int j = 0; j < 9; j++) { gauge[indexo + j] *= sign; }
      //if(i1 == 1 and i2 == 1 and i3 == 1 and i4 == 0){qmessage("==odd  %d \n", int(sign));}
    }
  }

  // Apply boundary conditions to temporal links
  if (param->t_boundary == QUDA_ANTI_PERIODIC_T && last_node_in_t()) {
    #pragma omp parallel for
    for (Int j = 0; j < Vh; j++) {
      Int sign = 1;
      if (dslash_type == QUDA_ASQTAD_DSLASH) {
        if (j >= (X4 - 3) * X1h * X2 * X3) { sign = -1; }
      } else {
        if (j >= (X4 - 1) * X1h * X2 * X3) { sign = -1; }
      }

      Long indexe = ((0 * Vh + j)*4 + 3)*9 ;
      Long indexo = ((1 * Vh + j)*4 + 3)*9 ;
      for (Int i = 0; i < 9; i++) {
        gauge[indexe + i] *= sign;
        gauge[indexo + i] *= sign;
      }
    }
  }
}


template <class Td, class T1>
void Ffield4d_to_quda_ff(T1*  quda_ff, qlat::FermionField4dT<Td>& ff, Int dir = 1)
{
  TIMER("Ffield4d_to_quda_ff(ff, quda_ff)");
  Qassert(ff.initialized);
  const Geometry& geo = ff.geo();
  //const WilsonVector* quda_pt =
  //    reinterpret_cast<const WilsonVector*>(qff.data());
  Long V = geo.local_volume();
  Long Vh = V / 2;

  qlat::ComplexT<Td>* src = (qlat::ComplexT<Td>*) qlat::get_data(ff).data();

  //
  #pragma omp parallel for
  for (Long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (Int dc = 0; dc < 12; dc++)
    {
      Long quda_idx = eo * Vh + qlat_idx_4d / 2;
      if(dir == 1){quda_ff[quda_idx*12 + dc] = src[qlat_idx_4d*12 + dc];}
      if(dir == 0){src[qlat_idx_4d*12 + dc] = quda_ff[quda_idx*12 + dc];}
    }
  }

}

template <class Td, class T1>
void quda_ff_to_Ffield4d(qlat::FermionField4dT<Td>& ff, T1* quda_ff)
{
  Ffield4d_to_quda_ff(quda_ff, ff, 0);
}

template <class Ty, class T1>
void qlat_cf_to_quda_cf_test(T1*  quda_cf, Ty* src, const Geometry& geo, Int Dim, Int dir = 1)
{
  TIMER("qlat_cf_to_quda_cf_test(qlat_cf, quda_cf)");
  Qassert(sizeof(T1) == sizeof(ComplexT<double>));// quda can only use double interface
  Long V = geo.local_volume();
  Long Vh = V / 2;
  #pragma omp parallel for
  for (Long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (Int dc = 0; dc < Dim; dc++)
    {
      Long quda_idx = eo * Vh + qlat_idx_4d / 2;
      if(dir == 1){quda_cf[quda_idx*Dim + dc] = src[qlat_idx_4d*Dim + dc];}
      if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[quda_idx*Dim + dc];}
    }
  }
}

template <class Ty, class T1>
void quda_cf_to_qlat_cf_test(Ty* qlat_cf, T1* quda_cf, const Geometry& geo, Int Dim)
{
  qlat_cf_to_quda_cf_test(quda_cf, qlat_cf, geo, Dim, 0);
}

template <class Ty, class T1>
void qlat_cf_to_quda_cf_test(T1*  quda_cf, colorFT& qlat_cf, Int dir = 1)
{
  Qassert(qlat_cf.initialized);
  const Geometry& geo = qlat_cf.geo();
  Long V = geo.local_volume();
  //Long Vh = V / 2;
  //Qassert(geo.multiplicity == 3);
  //const Long Dim = geo.multiplicity;
  const Long Dim = 3;

  Ty* src = (Ty*) qlat::get_data(qlat_cf).data();
  qlat_cf_to_quda_cf_test(quda_cf, src, geo, Dim, dir);
}

template <class Ty, class T1>
void quda_cf_to_qlat_cf_test(colorFT& qlat_cf, T1* quda_cf)
{
  qlat_cf_to_quda_cf_test(quda_cf, qlat_cf, 0);
}

//template <class Ty>
//void copy_color_prop(qlat::vector_gpu<Ty >& res, std::vector<colorFT >& src, Int dir = 1)
//{
//  Qassert(src.size() == 3);
//  qlat::vector<Ty* > srcP;srcP.resize(3);
//  for(Int ic=0;ic<3;ic++){
//    Qassert(src[ic].initialized);
//    srcP[ic] = (Ty*) qlat::get_data(src[ic]).data();
//  }
//
//  const qlat::Geometry& geo = src[0].geo();
//  const Long V = geo.local_volume();
//  Qassert(res.size() == V*9);
//  Ty* resP = res.data();
//
//  if(dir == 1){
//  qacc_for(isp, V, {
//    for(Int c0=0;c0<3;c0++)
//    for(Int c1=0;c1<3;c1++)
//    {
//      resP[isp*9 + c1*3 + c0 ] = srcP[c0][isp*3 + c1];
//    }
//  });}
//
//  if(dir == 0){
//  qacc_for(isp, V, {
//    for(Int c0=0;c0<3;c0++)
//    for(Int c1=0;c1<3;c1++)
//    {
//      srcP[c0][isp*3 + c1] = resP[isp*9 + c1*3 + c0];
//    }
//  });}
//}
//
//template <class Ty>
//void copy_to_color_prop(std::vector<colorFT >& res, qlat::vector_gpu<Ty >& src)
//{
//  copy_color_prop(src, res, 0);
//}


//inline void get_index_mappings(qlat::vector<Long >& map, const Geometry& geo)
//{
//  const Long V = geo.local_volume();
//  const Long Vh = V / 2;
//
//  if(map.size() == V){return ;}
//  else{map.resize(V);}
//  qacc_for(qlat_idx_4d, V , {
//    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
//    const Int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
//    const Long quda_idx = eo * Vh + qlat_idx_4d / 2;
//    map[qlat_idx_4d] = quda_idx;
//  });
//}

///////GPU order with color to be outside even odd
//template <class T1, class Ty, Int dir>
//void qlat_cf_to_quda_cfT(T1*  quda_cf, Ty* src, const Int Dim, const Geometry& geo, qlat::vector<Long >& map_)
//{
//  TIMER("qlat_cf_to_quda_cf");
//  const Long V = geo.local_volume();
//  Long Vh = V / 2;
//  if(map_.size() != V){get_index_mappings(map_, geo);}
//  qlat::vector<Long >& map = map_;
//  qacc_for(qlat_idx_4d, V, {
//    const Long quda_idx = map[qlat_idx_4d];
//    const Long eo = quda_idx / Vh;
//    const Long qi = quda_idx % Vh;
//    for(Int dc = 0; dc < Dim; dc++)
//    {
//      //if(dir == 1){quda_cf[ quda_idx*Dim + dc] = src[qlat_idx_4d*Dim + dc];}
//      //if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[quda_idx*Dim + dc];}
//      if(dir == 1){quda_cf[(eo*Dim + dc)*Vh + qi] = src[qlat_idx_4d*Dim + dc];}
//      if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[(eo*Dim + dc)*Vh + qi];}
//    }
//  });
//}

//template <class T1, class Ty>
//void qlat_cf_to_quda_cf(T1*  quda_cf, Ty* src, const Int Dim, const Geometry& geo, qlat::vector<Long >& map)
//{
//  qlat_cf_to_quda_cfT<T1, Ty, 1>(quda_cf, src, Dim, geo, map);
//}
//
//template <class T1, class Ty>
//void quda_cf_to_qlat_cf(Ty* res, T1*  quda_cf, const Int Dim, const Geometry& geo, qlat::vector<Long >& map)
//{
//  qlat_cf_to_quda_cfT<T1, Ty, 0>(quda_cf, res, Dim, geo, map);
//}

/*
  Always copy to double precions Quda fields first
  x and src may share the same pointers
*/
template <class Ty, Int dir>
void qlat_cf_to_quda_cfT(quda::ColorSpinorField& x, Ty* src, const Geometry& geo, qlat::vector<Long >& map)
{
  //quda::ColorSpinorParam cs_tmp(x);
  /////ColorSpinorField& x
  const Long V = geo.local_volume();
  const Int Ndata = x.Nspin() * x.Ncolor();//Ncolor()
  const size_t Vl = size_t(V) * Ndata;
  //const size_t Vb = Vl * 2;
  quda::ColorSpinorParam param(x);
  //Qassert( x.TotalBytes() % Vb == 0);
  QudaPrecision Dtype = param.Precision();
  Qassert(Dtype == QUDA_DOUBLE_PRECISION or Dtype == QUDA_SINGLE_PRECISION or Dtype == QUDA_HALF_PRECISION);
  //
  quda::ColorSpinorField *g = NULL;
  qlat::vector_gpu<int8_t >& quda_buf = get_vector_gpu_plan<int8_t >(0, "quda_field_copy_buffers", 1);
  //if(Dtype == QUDA_HALF_PRECISION or Dtype == QUDA_SINGLE_PRECISION)
  //if(Dtype != QUDA_DOUBLE_PRECISION)
  {
    quda_buf.resizeL(Vl * sizeof(RealD) * 2 / sizeof(int8_t));
    //quda::ColorSpinorParam param(x);
    param.setPrecision(QUDA_DOUBLE_PRECISION, QUDA_DOUBLE_PRECISION, true);
    param.is_composite  = false;
    param.is_component  = false;
    param.composite_dim = 1;
    param.create = QUDA_REFERENCE_FIELD_CREATE;
    param.v = (void*) quda_buf.data();
    g = quda::ColorSpinorField::Create(param);
  }
  if(x.IsComposite()){
    for(Int di=0;di<x.CompositeDim();di++)
    {
      if(dir == 0)*g = x.Component(di);
      qlat_cf_to_quda_cfT<qlat::ComplexT<double >, Ty, dir>((qlat::ComplexT<double >* ) g->data(), &src[di*Vl], Ndata, geo, map);
      if(dir == 1)x.Component(di) = *g;

      //if(Dtype == QUDA_DOUBLE_PRECISION){
      //  qlat_cf_to_quda_cfT<qlat::ComplexT<double>, Ty, dir>((qlat::ComplexT<double>* ) x.Component(di).data(), &src[di*Vl], Ndata, geo, map);
      //}else{
      //  if(dir == 0)*g = x.Component(di);
      //  qlat_cf_to_quda_cfT<qlat::ComplexT<double >, Ty, dir>((qlat::ComplexT<double >* ) g->data(), &src[di*Vl], Ndata, geo, map);
      //  if(dir == 1)x.Component(di) = *g;
      //}
      //if(Dtype == QUDA_SINGLE_PRECISION)
      //qlat_cf_to_quda_cfT<qlat::ComplexT<float >, Ty, dir>((qlat::ComplexT<float >* ) x.Component(di).data(), &src[di*Vl], Ndata, geo, map);
      //if(Dtype == QUDA_HALF_PRECISION)
      //{
      //if(dir == 0)*g = x.Component(di);
      //qlat_cf_to_quda_cfT<qlat::ComplexT<float >, Ty, dir>((qlat::ComplexT<float >* ) g->data(), &src[di*Vl], Ndata, geo, map);
      //if(dir == 1)x.Component(di) = *g;
      //}
    }
  }else{
    if(dir == 0)*g = x;
    qlat_cf_to_quda_cfT<qlat::ComplexT<double >, Ty, dir>((qlat::ComplexT<double >* ) g->data(), src, Ndata, geo, map);
    if(dir == 1)x = *g;
    //if(Dtype == QUDA_DOUBLE_PRECISION){
    //  qlat_cf_to_quda_cfT<qlat::ComplexT<double>, Ty, dir>((qlat::ComplexT<double>* ) x.data(), src, Ndata, geo, map);
    //}else{
    //  if(dir == 0)*g = x;
    //  qlat_cf_to_quda_cfT<qlat::ComplexT<double >, Ty, dir>((qlat::ComplexT<double >* ) g->data(), src, Ndata, geo, map);
    //  if(dir == 1)x = *g;
    //}
    //if(Dtype == QUDA_SINGLE_PRECISION)
    //qlat_cf_to_quda_cfT<qlat::ComplexT<float >, Ty, dir>((qlat::ComplexT<float >* ) x.data(), src, Ndata, geo, map);
    //if(Dtype == QUDA_HALF_PRECISION)
    //{
    //if(dir == 0)*g = x;
    //qlat_cf_to_quda_cfT<qlat::ComplexT<float >, Ty, dir>((qlat::ComplexT<float >* ) g->data(), src, Ndata, geo, map);
    //if(dir == 1)x = *g;
    //}
  }
  if(g != NULL){delete g;g=NULL;}
}

template <class Ty>
void qlat_cf_to_quda_cf(quda::ColorSpinorField& x, Ty* src, const Geometry& geo, qlat::vector<Long >& map)
{
  qlat_cf_to_quda_cfT<Ty, 1>(x, src, geo, map);
}

template <class Ty>
void quda_cf_to_qlat_cf(Ty* src, quda::ColorSpinorField& x, const Geometry& geo, qlat::vector<Long >& map)
{
  qlat_cf_to_quda_cfT<Ty, 0>(x, src, geo, map);
}


}  // namespace qlat

#endif


