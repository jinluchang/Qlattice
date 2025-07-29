#pragma once

/*
 * The eigen system is stored in compressed format invented by Christoph Lehner
 * and Chulwoo Jung.
 *
 * Coded initially from https://github.com/lehner/eigen-comp by Christoph Lehner
 *
 * Added to Qlattice by Luchang Jin.
 *
 * Add the functionality to read from a different geometry.
 *
 */

#include <qlat-utils/cache.h>
#include <qlat-utils/matrix.h>
#include <qlat-utils/mvector.h>
#include <qlat/core.h>
#include <qlat/fermion-action.h>
#include <qlat/field-dist-io.h>
#include <qlat/field-expand.h>
#include <qlat/field-fft.h>
#include <qlat/field-serial-io.h>
#include <qlat/field.h>
#include <qlat/mpi.h>
#include <qlat/qcd-gauge-transformation.h>
#include <qlat/qcd-smear.h>
#include <qlat/qcd-topology.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>
#include <qlat/utils-coordinate.h>
#include <qlat/utils-io.h>

namespace qlat
{  //

inline int fp_map(const float in, const float min, const float max, const int N)
{
  // Idea:
  //
  // min=-6
  // max=6
  //
  // N=1
  // [-6,0] -> 0, [0,6] -> 1;  reconstruct 0 -> -3, 1-> 3
  //
  // N=2
  // [-6,-2] -> 0, [-2,2] -> 1, [2,6] -> 2;  reconstruct 0 -> -4, 1->0, 2->4
  int ret = (int)((float)(N + 1) * ((in - min) / (max - min)));
  if (ret == N + 1) {
    ret = N;
  }
  return ret;
}

inline float fp_unmap(const int val, const float min, const float max,
                      const int N)
{
  return min + (float)(val + 0.5) * (max - min) / (float)(N + 1);
}

inline float unmap_fp16_exp(unsigned short e)
{
  const double base = 1.4142135623730950488;
  const float de = (float)((int)e - USHRT_MAX / 2);
  return pow(base, de);
}

inline Long fp_16_size(const Long f_size, const Long nsc)
{
  return (f_size + f_size / nsc) * sizeof(uint16_t);
}

void read_floats(Vector<float> out, const Vector<uint8_t> fp_data);

void read_floats_fp16(float* out, const uint8_t* ptr, const int64_t n,
                      const int nsc);

void read_floats_fp16(Vector<float> out, const Vector<uint8_t> fp_data,
                      const int nsc);

struct CompressedEigenSystemInfo {
  bool initialized;
  // int s[5]; // local vol size
  // int b[5]; // block size
  // // derived
  // int nb[5];
  // int blocks;
  // //
  // int index;
  //
  int nkeep_fp16;          // nkeep - nkeep_single
  Coordinate total_node;   // number of node in each direction (total)
  Coordinate total_block;  // number of block in each direction (total)
  Coordinate node_block;   // number of block in a node in each direction (node)
  //
  Coordinate total_site;  // number of site in each direction (total)
  Coordinate node_site;   // number of site in a node in each direction (node)
  Coordinate
      block_site;  // number of site in each block in each direction (block)
  int ls;
  int neig;
  int nkeep;         // number base
  int nkeep_single;  // number base stored as single precision
  int FP16_COEF_EXP_SHARE_FLOATS;
  //
  std::vector<crc32_t> crcs;
  //
  CompressedEigenSystemInfo() { init(); }
  //
  void init() { initialized = false; }
};

struct CompressedEigenSystemDenseInfo {
  Coordinate total_site;
  Coordinate node_site;
  Coordinate block_site;
  int ls;
  int neig;
  int nkeep;
  int nkeep_single;
  int FP16_COEF_EXP_SHARE_FLOATS;
};

CompressedEigenSystemInfo populate_eigen_system_info(
    const CompressedEigenSystemDenseInfo& cesdi,
    const std::vector<crc32_t>& crcs);

CompressedEigenSystemInfo read_compressed_eigen_system_info(
    const std::string& root);

void write_compressed_eigen_system_info(const CompressedEigenSystemInfo& cesi,
                                        const std::string& root);

CompressedEigenSystemInfo resize_compressed_eigen_system_info(
    const CompressedEigenSystemInfo& cesi, const Coordinate& new_size_node);

struct HalfVector : Field<ComplexF> {
  static const int c_size = 12;  // number of complex number per wilson vector
  int ls;
  // geo.multiplicity = ls * c_size;
};

void init_half_vector(HalfVector& hv, const Geometry& geo, const int ls);

struct CompressedEigenSystemData : Field<uint8_t> {
  CompressedEigenSystemInfo cesi;
  Long block_vol_eo;
  int ls;
  Long basis_c_size;
  Long coef_c_size;
  Long basis_size_single;
  Long basis_size_fp16;
  Long coef_size_single;
  Long coef_size_fp16;
  Long coef_size;
  Long bases_size_single;
  Long bases_size_fp16;
  Long bases_offset_single;
  Long bases_offset_fp16;
  Long coefs_offset;
  Long end_offset;
};

Geometry get_geo_from_cesi(const CompressedEigenSystemInfo& cesi,
                           const int id_node,
                           const Coordinate& new_size_node = Coordinate());

Geometry block_geometry(const Geometry& geo_full, const Coordinate& block_site);

void init_compressed_eigen_system_data(
    CompressedEigenSystemData& cesd, const CompressedEigenSystemInfo& cesi,
    const int id_node, const Coordinate& new_size_node = Coordinate());

struct CompressedEigenSystemBases : Field<ComplexF> {
  int n_basis;
  Long block_vol_eo;  // even odd precondition
                      // (block_vol_eo is half of the number of site within the
                      // block n_t * n_z * n_y * n_x/2) n_t slowest varying and
                      // n_x fastest varying
  int ls;             // ls varying faster than n_x
  int c_size_vec;     // c_size_vec = block_vol_eo * ls * HalfVector::c_size;
  // geo.multiplicity = n_basis * c_size_vec;
  Geometry geo_full;
  Coordinate block_site;
};

struct CompressedEigenSystemCoefs : Field<ComplexF> {
  int n_vec;
  int n_basis;
  int ls;
  int c_size_vec;  // c_size_vec = n_basis;
  // geo.multiplicity = n_vec * c_size_vec;
  Geometry geo_full;
  Coordinate block_site;
};

void init_compressed_eigen_system_bases(CompressedEigenSystemBases& cesb,
                                        const int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const int ls);

void init_compressed_eigen_system_bases(
    CompressedEigenSystemBases& cesb, const CompressedEigenSystemInfo& cesi,
    const int id_node, const Coordinate& new_size_node = Coordinate());

void init_compressed_eigen_system_coefs(CompressedEigenSystemCoefs& cesc,
                                        const int n_vec, const int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const int ls);

void init_compressed_eigen_system_coefs(
    CompressedEigenSystemCoefs& cesc, const CompressedEigenSystemInfo& cesi,
    const int id_node, const Coordinate& new_size_node = Coordinate());

API inline Long& get_vbfile_buffer_limit()
{
  static Long buffer_limit = 4L * 1024L * 1024L * 1024L;
  return buffer_limit;
}

struct VBFile
// virtual bufferred IO (real IO performed during vbflush)
{
  std::string fn;
  std::string mode;
  QFile fp;
  Long buffer_limit;
  std::vector<uint8_t> buffer;
  Long entry_total_size;
  std::vector<Vector<uint8_t>> entries;
  //
  VBFile()
  {
    fp.init();
    buffer_limit = get_vbfile_buffer_limit();
    entry_total_size = 0;
  }
  //
  ~VBFile() { assert(fp.null()); }
};

VBFile vbopen(const std::string& fn, const std::string& mode);

void vbflush(VBFile& fp);

void vbclose(VBFile& fp);

void vbread_data(const Vector<uint8_t>& v, VBFile& fp);

void vbwrite_data(const Vector<uint8_t>& v, VBFile& fp);

int vbseek(VBFile& fp, const Long offset, const int whence);

struct VFile
// virtual IO (real IO performed during vbclose)
{
  std::string fn;
  std::string mode;
  Long pos;
  Long size;
  std::multimap<Long, Vector<uint8_t>> read_entries;
  std::multimap<Long, Vector<uint8_t>> write_entries;
};

void set_vfile_size(VFile& fp);

VFile vopen(const std::string& fn, const std::string& mode);

void vclose(VFile& fp);

template <class M>
Long vread_data(const Vector<M>& v, VFile& fp)
{
  qassert(fp.mode == "r");
  Vector<uint8_t> dv((uint8_t*)v.data(), v.data_size());
  fp.read_entries.insert(std::pair<Long, Vector<uint8_t>>(fp.pos, dv));
  fp.pos += dv.size();
  qassert(0 <= fp.pos);
  return dv.size();
}

template <class M>
Long vwrite_data(const Vector<M>& v, VFile& fp)
{
  qassert(fp.mode == "w");
  Vector<uint8_t> dv((uint8_t*)v.data(), v.data_size());
  fp.write_entries.insert(std::pair<Long, Vector<uint8_t>>(fp.pos, dv));
  fp.pos += dv.size();
  return dv.size();
}

int vseek(VFile& fp, const Long offset, const int whence);

Long vtell(const VFile& fp);

void load_block_data(CompressedEigenSystemData& cesd, const Coordinate& xl,
                     const CompressedEigenSystemInfo& cesi, VFile& fp,
                     const Coordinate& xl_file);

void save_block_data(const CompressedEigenSystemData& cesd,
                     const Coordinate& xl,
                     const CompressedEigenSystemInfo& cesi, VFile& fp);

crc32_t block_data_crc(const CompressedEigenSystemData& cesd,
                       const Coordinate& xl,
                       const CompressedEigenSystemInfo& cesi,
                       const Coordinate& xl_file);

std::vector<crc32_t> load_node_data(CompressedEigenSystemData& cesd,
                                    const CompressedEigenSystemInfo& cesi,
                                    const std::string& path);

crc32_t save_node_data(const CompressedEigenSystemData& cesd,
                       const CompressedEigenSystemInfo& cesi,
                       const std::string& path);

void load_block(CompressedEigenSystemBases& cesb,
                CompressedEigenSystemCoefs& cesc,
                const CompressedEigenSystemData& cesd, const Coordinate& xl,
                const CompressedEigenSystemInfo& cesi);

std::vector<crc32_t> load_node(CompressedEigenSystemBases& cesb,
                               CompressedEigenSystemCoefs& cesc,
                               const CompressedEigenSystemInfo& cesi,
                               const std::string& path);

std::vector<double> read_eigen_values(const std::string& path);

struct BlockedHalfVector : Field<ComplexF> {
  Long block_vol_eo;  // even odd precondition (block_vol is half of the number
                      // of site within the block)
  int ls;
  // geo.multiplicity = block_vol_eo * ls * HalfVector::c_size;
  Geometry geo_full;
  Coordinate block_site;
};

void init_blocked_half_vector(BlockedHalfVector& bhv, const Geometry& geo_full,
                              const Coordinate& block_site, const int ls);

void decompress_eigen_system(std::vector<BlockedHalfVector>& bhvs,
                             const CompressedEigenSystemBases& cesb,
                             const CompressedEigenSystemCoefs& cesc);

void convert_half_vector(BlockedHalfVector& bhv, const HalfVector& hv,
                         const Coordinate& block_site);

void convert_half_vector(HalfVector& hv, const BlockedHalfVector& bhv);

void convert_half_vectors(std::vector<HalfVector>& hvs,
                          std::vector<BlockedHalfVector>& bhvs);

void convert_half_vector_bfm_format(Vector<ComplexF> bfm_data,
                                    const HalfVector& hv);

Long load_compressed_eigen_vectors(vector<double>& eigen_values,
                                   CompressedEigenSystemInfo& cesi,
                                   CompressedEigenSystemBases& cesb,
                                   CompressedEigenSystemCoefs& cesc,
                                   const std::string& path);

crc32_t save_half_vectors(const std::vector<HalfVector>& hvs,
                          const std::string& fn,
                          const bool is_saving_crc = false,
                          const bool is_bfm_format = false);

Long decompress_eigen_vectors_node(
    const std::string& old_path, const CompressedEigenSystemInfo& cesi,
    const std::string& new_path, const int idx,
    const Coordinate& new_size_node = Coordinate());

crc32_t resize_compressed_eigen_vectors_node(
    std::vector<crc32_t>& crcs_acc, const std::string& old_path,
    const CompressedEigenSystemInfo& cesi, const std::string& new_path,
    const int idx, const Coordinate& size_node);

void combine_crc32(const std::string& path, const int idx_size,
                   const CompressedEigenSystemInfo& cesi);

void decompress_eigen_vectors(const std::string& old_path,
                              const std::string& new_path,
                              const Coordinate& new_size_node = Coordinate());

bool check_compressed_eigen_vectors(const std::string& path);

bool resize_compressed_eigen_vectors(const std::string& old_path,
                                     const std::string& new_path,
                                     const Coordinate& size_node);

void decompressed_eigen_vectors_check_crc32(const std::string& path);

bool eigen_system_repartition(const Coordinate& new_size_node,
                              const std::string& path,
                              const std::string& new_path = "");

}  // namespace qlat
