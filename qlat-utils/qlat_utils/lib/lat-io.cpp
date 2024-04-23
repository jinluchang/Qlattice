#include <qlat-utils/lat-io.h>

namespace qlat
{  //

bool is_lat_info_complex(const LatInfo& info)
{
  if ((Long)info.size() < 1) {
    return false;
  }
  const LatDim& dim = info.back();
  if (dim.name != "re-im") {
    return false;
  } else {
    qassert(dim.size == 2);
    qassert(dim.indices.size() == 2);
    qassert(dim.indices[0] == "re");
    qassert(dim.indices[1] == "im");
    return true;
  }
}

LatDim lat_dim_re_im()
{
  return lat_dim_string("re-im", make_array<std::string>("re", "im"));
}

LatDim lat_dim_number(const std::string& name, const Long start, const Long end,
                      const Long inc)
{
  LatDim dim;
  dim.name = name;
  if (start == 0 and inc == 1) {
    dim.size = end + 1;
  } else {
    for (Long i = start; i <= end; i += inc) {
      dim.size += 1;
      dim.indices.push_back(show(i));
    }
  }
  return dim;
}

std::string show(const LatDim& dim)
{
  std::ostringstream out;
  out << ssprintf("\"%s\"[%ld]:", dim.name.c_str(), dim.size);
  for (Long i = 0; i < (Long)dim.indices.size(); ++i) {
    out << ssprintf(" \"%s\"", dim.indices[i].c_str());
  }
  return out.str();
}

std::string show(const LatInfo& info)
{
  std::ostringstream out;
  out << ssprintf("ndim: %ld\n", info.size());
  for (int i = 0; i < (int)info.size(); ++i) {
    out << show(info[i]) << "\n";
  }
  return out.str();
}

LatDim read_lat_dim(const std::string& str)
{
  LatDim dim;
  Long cur = 0;
  char c;
  if (!parse_string(dim.name, cur, str)) {
    qerr("read_lat_dim: dim.name");
  } else if (!parse_char(c, cur, str) or c != '[') {
    qerr("read_lat_dim: [");
  } else if (!parse_long(dim.size, cur, str)) {
    qerr("read_lat_dim: Long");
  } else if (!parse_char(c, cur, str) or c != ']') {
    qerr("read_lat_dim: ]");
  } else if (!parse_char(c, cur, str) or c != ':') {
    qerr("read_lat_dim: :");
  } else {
    while (parse_char(c, cur, str)) {
      if (c == '\n') {
        qassert(cur == (Long)str.size());
        break;
      }
      qassert(c == ' ');
      std::string index;
      if (!parse_string(index, cur, str)) {
        qerr("read_lat_dim: string");
      }
      dim.indices.push_back(index);
    }
  }
  return dim;
}

LatInfo read_lat_info(const std::string& str)
{
  LatInfo info;
  const std::vector<std::string> infos = split_into_lines(str);
  qassert(infos.size() >= 1);
  const std::string ndim_prop = "ndim: ";
  qassert(infos[0].compare(0, ndim_prop.size(), ndim_prop) == 0);
  const Long ndim = read_long(std::string(infos[0], ndim_prop.size()));
  qassert(ndim == (Long)infos.size() - 1);
  for (int i = 1; i < (int)infos.size(); ++i) {
    info.push_back(read_lat_dim(infos[i]));
  }
  return info;
}

// -----------------------

template <class T>
const std::string& get_lat_data_header();

template <class T>
void lat_data_load(LatDataT<T>& ld, QFile& qfile)
{
  qassert(not qfile.null());
  std::vector<char> check_line(get_lat_data_header<T>().size(), 0);
  const Long fread_check_len =
      qfread(check_line.data(), get_lat_data_header<T>().size(), 1, qfile);
  if (fread_check_len == 0) {
    ld.init();
    return;
  }
  qassert(fread_check_len == 1);
  qassert(std::string(check_line.data(), check_line.size()) ==
          get_lat_data_header<T>());
  std::vector<std::string> infos;
  infos.push_back(get_lat_data_header<T>());
  while (infos.back() != "END_HEADER\n" && infos.back() != "") {
    infos.push_back(qgetline(qfile));
  }
  std::ostringstream out;
  for (int i = 3; i < (int)infos.size() - 2; ++i) {
    out << infos[i];
  }
  const std::string info_str = out.str();
  ld.info = read_lat_info(info_str);
  const std::string& crc_str = infos[infos.size() - 2];
  const std::string crc_prop = "crc32: ";
  qassert(crc_str.compare(0, crc_prop.size(), crc_prop) == 0);
  const crc32_t crc = read_crc32(std::string(crc_str, crc_prop.size()));
  lat_data_alloc(ld);
  qassert((Long)ld.res.size() == lat_info_size(ld.info));
  qassert((Long)ld.res.size() * (Long)sizeof(T) == read_long(infos[2]));
  const Long fread_res_len =
      qfread(ld.res.data(), sizeof(T), ld.res.size(), qfile);
  qassert(fread_res_len == (Long)ld.res.size());
  const crc32_t crc_computed =
      crc32_par(ld.res.data(), ld.res.size() * sizeof(T));
  if (crc != crc_computed) {
    qerr(ssprintf("ERROR: crc do not match: file=%08X computed=%08X path='%s'.",
                  crc, crc_computed, qfile.path().c_str()));
  }
  to_from_little_endian(get_data(ld.res));
}

template <class T>
void lat_data_save(const LatDataT<T>& ld, QFile& qfile)
{
  qassert(not qfile.null());
  std::vector<T> res_copy;
  if (!is_little_endian()) {
    res_copy = ld.res;
    qassert(res_copy.size() == ld.res.size());
    to_from_little_endian(get_data(res_copy));
  }
  const std::string data_size =
      ssprintf("data_size\n%ld\n", ld.res.size() * sizeof(T));
  const std::string info_str = show(ld.info);
  const std::string checksum_str =
      ssprintf("crc32: %08X\n",
               crc32_par(is_little_endian() ? ld.res.data() : res_copy.data(),
                         ld.res.size() * sizeof(T)));
  const std::string end_header = "END_HEADER\n";
  qfwrite(get_lat_data_header<T>().data(), get_lat_data_header<T>().size(), 1, qfile);
  qfwrite(data_size.data(), data_size.size(), 1, qfile);
  qfwrite(info_str.data(), info_str.size(), 1, qfile);
  qfwrite(checksum_str.data(), checksum_str.size(), 1, qfile);
  qfwrite(end_header.data(), end_header.size(), 1, qfile);
  qfwrite(is_little_endian() ? ld.res.data() : res_copy.data(), sizeof(T),
          ld.res.size(), qfile);
}

// -----------------------

template <>
const std::string& get_lat_data_header<Int>()
{
  return lat_data_int_header;
}

template <>
void LatDataT<Int>::load(QFile& qfile)
{
  lat_data_load(*this, qfile);
}

template <>
void LatDataT<Int>::save(QFile& qfile) const
{
  lat_data_save(*this, qfile);
}

// -----------------------

template <>
const std::string& get_lat_data_header<Long>()
{
  return lat_data_long_header;
}

template <>
void LatDataT<Long>::load(QFile& qfile)
{
  lat_data_load(*this, qfile);
}

template <>
void LatDataT<Long>::save(QFile& qfile) const
{
  lat_data_save(*this, qfile);
}

// -----------------------

template <>
const std::string& get_lat_data_header<RealD>()
{
  return lat_data_header;
}

template <>
void LatDataT<RealD>::load(QFile& qfile)
{
  lat_data_load(*this, qfile);
}

template <>
void LatDataT<RealD>::save(QFile& qfile) const
{
  lat_data_save(*this, qfile);
}

// -----------------------

template <>
const std::string& get_lat_data_header<RealF>()
{
  return lat_data_real_f_header;
}

template <>
void LatDataT<RealF>::load(QFile& qfile)
{
  lat_data_load(*this, qfile);
}

template <>
void LatDataT<RealF>::save(QFile& qfile) const
{
  lat_data_save(*this, qfile);
}

// -----------------------

LatData::LatData(const LatDataRealF& ld) { *this = ld; }

LatData& LatData::operator=(const LatDataRealF& ld)
{
  info = ld.info;
  const Long size = ld.res.size();
  res.resize(size);
  for (Long i = 0; i < size; ++i) {
    res[i] = ld.res[i];
  }
  return *this;
}

LatDataRealF::LatDataRealF(const LatData& ld) { *this = ld; }

LatDataRealF& LatDataRealF::operator=(const LatData& ld)
{
  info = ld.info;
  const Long size = ld.res.size();
  res.resize(size);
  for (Long i = 0; i < size; ++i) {
    res[i] = ld.res[i];
  }
  return *this;
}

// -----------------------

std::string show_real(const LatData& ld)
{
  std::ostringstream out;
  const LatInfo& info = ld.info;
  out << "# ";
  for (int a = 0; a < (int)info.size(); ++a) {
    if (0 == a) {
      out << ssprintf("%10s ", info[a].name.c_str());
    } else {
      out << ssprintf("%12s ", info[a].name.c_str());
    }
  }
  out << ssprintf("%24s\n", "VALUE");
  std::vector<Long> idx(info.size(), 0);
  for (Long k = 0; k < lat_data_size(ld); ++k) {
    for (int a = 0; a < (int)info.size(); ++a) {
      out << ssprintf("%12s ", idx_name(info[a], idx[a]).c_str());
    }
    const Vector<RealD> ldv = lat_data_get_const(ld, idx);
    out << ssprintf("%24.17E\n", ldv[0]);
    idx[(int)info.size() - 1] += 1;
    for (int a = (int)info.size() - 1; a > 0; --a) {
      if (idx[a] == info[a].size) {
        idx[a] = 0;
        idx[a - 1] += 1;
      }
    }
  }
  return out.str();
}

std::string show_complex(const LatData& ld)
{
  std::ostringstream out;
  const LatInfo& info = ld.info;
  out << "# ";
  for (int a = 0; a < (int)info.size() - 1; ++a) {
    if (0 == a) {
      out << ssprintf("%10s ", info[a].name.c_str());
    } else {
      out << ssprintf("%12s ", info[a].name.c_str());
    }
  }
  out << ssprintf("%24s %24s\n", "RE-VALUE", "IM-VALUE");
  std::vector<Long> idx((int)info.size() - 1, 0);
  for (Long k = 0; k < lat_data_size(ld) / 2; ++k) {
    for (int a = 0; a < (int)info.size() - 1; ++a) {
      out << ssprintf("%12s ", idx_name(info[a], idx[a]).c_str());
    }
    const Vector<ComplexD> ldv = lat_data_cget_const(ld, idx);
    out << ssprintf("%24.17E %24.17E\n", ldv[0].real(), ldv[0].imag());
    if ((int)info.size() - 2 >= 0) {
      idx[(int)info.size() - 2] += 1;
    }
    for (int a = (int)info.size() - 2; a > 0; --a) {
      if (idx[a] == info[a].size) {
        idx[a] = 0;
        idx[a - 1] += 1;
      }
    }
  }
  return out.str();
}

std::string show(const LatData& ld)
{
  const LatInfo& info = ld.info;
  if (is_lat_info_complex(info)) {
    return show_complex(ld);
  } else {
    return show_real(ld);
  }
}

void print(const LatData& ld)
{
  const LatInfo& info = ld.info;
  display(ssprintf("%s", show(info).c_str()));
  std::vector<Long> idx(info.size(), 0);
  for (Long k = 0; k < lat_data_size(ld); ++k) {
    for (int a = 0; a < (int)info.size(); ++a) {
      display(ssprintf("%s[%8s] ", info[a].name.c_str(),
                       idx_name(info[a], idx[a]).c_str()));
    }
    display(ssprintf("%24.17E\n", ld.res[lat_data_offset(ld, idx)]));
    idx[info.size() - 1] += 1;
    for (int a = info.size() - 1; a > 0; --a) {
      if (idx[a] == info[a].size) {
        idx[a] = 0;
        idx[a - 1] += 1;
      }
    }
  }
}

const LatData& operator*=(LatData& ld, const RealD factor)
{
  Vector<RealD> v = lat_data_get(ld);
  for (Long i = 0; i < v.size(); ++i) {
    v[i] *= factor;
  }
  return ld;
}

const LatData& operator*=(LatData& ld, const ComplexD& factor)
{
  Vector<ComplexD> v = lat_data_cget(ld);
  for (Long i = 0; i < v.size(); ++i) {
    v[i] *= factor;
  }
  return ld;
}

LatData operator*(const LatData& ld, const RealD factor)
{
  LatData ret = ld;
  ret *= factor;
  return ret;
}

LatData operator*(const RealD factor, const LatData& ld)
{
  return ld * factor;
}

LatData operator*(const LatData& ld, const ComplexD& factor)
{
  LatData ret = ld;
  ret *= factor;
  return ret;
}

LatData operator*(const ComplexD& factor, const LatData& ld)
{
  return ld * factor;
}

const LatData& operator+=(LatData& ld, const LatData& ld1)
{
  if (is_zero(ld)) {
    ld = ld1;
  } else {
    if (not is_matching(ld, ld1)) {
      displayln("operator+=(ld,ld1): ld.info: " + show(ld.info));
      displayln("operator+=(ld,ld1): ld1.info: " + show(ld1.info));
      qassert(false);
    }
    for (Long i = 0; i < (Long)ld.res.size(); ++i) {
      ld.res[i] += ld1.res[i];
    }
  }
  return ld;
}

const LatData& operator-=(LatData& ld, const LatData& ld1)
{
  if (is_zero(ld)) {
    ld.info = ld1.info;
    lat_data_alloc(ld);
    for (Long i = 0; i < (Long)ld.res.size(); ++i) {
      ld.res[i] = -ld1.res[i];
    }
  } else {
    if (not is_matching(ld, ld1)) {
      displayln("operator-=(ld,ld1): ld.info: " + show(ld.info));
      displayln("operator-=(ld,ld1): ld1.info: " + show(ld1.info));
      qassert(false);
    }
    for (Long i = 0; i < (Long)ld.res.size(); ++i) {
      ld.res[i] -= ld1.res[i];
    }
  }
  return ld;
}

LatData operator+(const LatData& ld, const LatData& ld1)
{
  LatData ret = ld;
  ret += ld1;
  return ret;
}

LatData operator-(const LatData& ld, const LatData& ld1)
{
  LatData ret = ld;
  ret -= ld1;
  return ret;
}

LatData operator-(const LatData& ld)
{
  LatData ret;
  ret -= ld;
  return ret;
}

}  // namespace qlat
