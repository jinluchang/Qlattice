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

inline Int fp_map(const RealF in, const RealF min, const RealF max, const Int N)
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
  Int ret = (int)((RealF)(N + 1) * ((in - min) / (max - min)));
  if (ret == N + 1) {
    ret = N;
  }
  return ret;
}

inline RealF fp_unmap(const Int val, const RealF min, const RealF max,
                      const Int N)
{
  return min + (RealF)(val + 0.5) * (max - min) / (RealF)(N + 1);
}

inline RealF unmap_fp16_exp(unsigned short e)
{
  const RealD base = 1.4142135623730950488;
  const RealF de = (RealF)((int)e - USHRT_MAX / 2);
  return pow(base, de);
}

inline Long fp_16_size(const Long f_size, const Long nsc)
{
  return (f_size + f_size / nsc) * sizeof(uint16_t);
}

void read_floats(Vector<RealF> out, const Vector<uint8_t> fp_data);

void read_floats_fp16(RealF* out, const uint8_t* ptr, const int64_t n,
                      const Int nsc);

void read_floats_fp16(Vector<RealF> out, const Vector<uint8_t> fp_data,
                      const Int nsc);

struct CompressedEigenSystemInfo {
  bool initialized;
  // Int s[5]; // local vol size
  // Int b[5]; // block size
  // // derived
  // Int nb[5];
  // Int blocks;
  // //
  // Int index;
  //
  Int nkeep_fp16;          // nkeep - nkeep_single
  Coordinate total_node;   // number of node in each direction (total)
  Coordinate total_block;  // number of block in each direction (total)
  Coordinate node_block;   // number of block in a node in each direction (node)
  //
  Coordinate total_site;  // number of site in each direction (total)
  Coordinate node_site;   // number of site in a node in each direction (node)
  Coordinate
      block_site;  // number of site in each block in each direction (block)
  Int ls;
  Int neig;
  Int nkeep;         // number base
  Int nkeep_single;  // number base stored as single precision
  Int FP16_COEF_EXP_SHARE_FLOATS;
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
  Int ls;
  Int neig;
  Int nkeep;
  Int nkeep_single;
  Int FP16_COEF_EXP_SHARE_FLOATS;
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
  static const Int c_size = 12;  // number of complex number per wilson vector
  Int ls;
  // geo.multiplicity = ls * c_size;
};

void init_half_vector(HalfVector& hv, const Geometry& geo, const Int ls);

struct CompressedEigenSystemData : Field<uint8_t> {
  CompressedEigenSystemInfo cesi;
  Long block_vol_eo;
  Int ls;
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
                           const Int id_node,
                           const Coordinate& new_size_node = Coordinate());

Geometry block_geometry(const Geometry& geo_full, const Coordinate& block_site);

void init_compressed_eigen_system_data(
    CompressedEigenSystemData& cesd, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node = Coordinate());

struct CompressedEigenSystemBases : Field<ComplexF> {
  Int n_basis;
  Long block_vol_eo;  // even odd precondition
                      // (block_vol_eo is half of the number of site within the
                      // block n_t * n_z * n_y * n_x/2) n_t slowest varying and
                      // n_x fastest varying
  Int ls;             // ls varying faster than n_x
  Int c_size_vec;     // c_size_vec = block_vol_eo * ls * HalfVector::c_size;
  // geo.multiplicity = n_basis * c_size_vec;
  Geometry geo_full;
  Coordinate block_site;
};

struct CompressedEigenSystemCoefs : Field<ComplexF> {
  Int n_vec;
  Int n_basis;
  Int ls;
  Int c_size_vec;  // c_size_vec = n_basis;
  // geo.multiplicity = n_vec * c_size_vec;
  Geometry geo_full;
  Coordinate block_site;
};

void init_compressed_eigen_system_bases(CompressedEigenSystemBases& cesb,
                                        const Int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const Int ls);

void init_compressed_eigen_system_bases(
    CompressedEigenSystemBases& cesb, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node = Coordinate());

void init_compressed_eigen_system_coefs(CompressedEigenSystemCoefs& cesc,
                                        const Int n_vec, const Int n_basis,
                                        const Geometry& geo_full,
                                        const Coordinate& block_site,
                                        const Int ls);

void init_compressed_eigen_system_coefs(
    CompressedEigenSystemCoefs& cesc, const CompressedEigenSystemInfo& cesi,
    const Int id_node, const Coordinate& new_size_node = Coordinate());

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

Int vbseek(VBFile& fp, const Long offset, const Int whence);

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

Int vseek(VFile& fp, const Long offset, const Int whence);

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

std::vector<RealD> read_eigen_values(const std::string& path);

struct BlockedHalfVector : Field<ComplexF> {
  Long block_vol_eo;  // even odd precondition (block_vol is half of the number
                      // of site within the block)
  Int ls;
  // geo.multiplicity = block_vol_eo * ls * HalfVector::c_size;
  Geometry geo_full;
  Coordinate block_site;
};

void init_blocked_half_vector(BlockedHalfVector& bhv, const Geometry& geo_full,
                              const Coordinate& block_site, const Int ls);

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

Long load_compressed_eigen_vectors(vector<RealD>& eigen_values,
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
    const std::string& new_path, const Int idx,
    const Coordinate& new_size_node = Coordinate());

crc32_t resize_compressed_eigen_vectors_node(
    std::vector<crc32_t>& crcs_acc, const std::string& old_path,
    const CompressedEigenSystemInfo& cesi, const std::string& new_path,
    const Int idx, const Coordinate& size_node);

void combine_crc32(const std::string& path, const Int idx_size,
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
