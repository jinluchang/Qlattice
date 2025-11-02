#pragma once

#include <qlat/field-base-io.h>
#include <qlat/field-dist-io.h>

#include <stdio.h>
#include <ctime>

#include <fstream>
#include <iostream>

namespace qlat
{  //

inline Coordinate get_default_serial_new_size_node(const Geometry& geo, const Int max_num_ = 0)
{
  const Int num_node = geo.geon.num_node;
  const Int max_num =
      (max_num_ <= 0) or (max_num_ > num_node) ? num_node : max_num_;
  const Coordinate total_site = geo.total_site();
  Coordinate new_size_node = Coordinate(1, 1, 1, total_site[3]);
  while (max_num < new_size_node[3]) {
    if (new_size_node[3] % 2 == 0) {
      new_size_node[3] /= 2;
    } else if (new_size_node[3] % 3 == 0) {
      new_size_node[3] /= 3;
    } else if (new_size_node[3] % 5 == 0) {
      new_size_node[3] /= 5;
    } else if (new_size_node[3] % 7 == 0) {
      new_size_node[3] /= 7;
    } else if (new_size_node[3] % 11 == 0) {
      new_size_node[3] /= 11;
    } else if (new_size_node[3] % 13 == 0) {
      new_size_node[3] /= 13;
    } else {
      new_size_node[3] = 1;
    }
  }
  qassert(total_site % new_size_node == Coordinate());
  return new_size_node;
}

template <class M>
Long serial_write_field(const Field<M>& f, const std::string& path,
                        const Coordinate& new_size_node)
// will append to the file
// assume new_size_node is properly chosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_write_field");
  std::vector<Field<M> > fs;
  shuffle_field(fs, f, new_size_node);
  const Int mpi_tag = 6;
  if (get_id_node() == 0) {
    qassert(fs.size() > 0);
    Field<M> f;
    f = fs[0];
    Vector<M> v = get_data(f);
    const Int num_node = get_num_node();
    const Int new_num_node = product(new_size_node);
    QFile qfile = qfopen(path, "a");
    qassert(not qfile.null());
    for (Int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const Int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      if (0 == id_node) {
        assign(v, get_data(fs[new_id_node]));
      } else {
        mpi_recv(v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm());
      }
      qwrite_data(v, qfile);
    }
    qfclose(qfile);
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      const Vector<M> v = get_data(fs[i]);
      mpi_send((void*)v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag,
               get_comm());
    }
  }
  const Long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
Long serial_read_field(Field<M>& f, const std::string& path,
                       const Coordinate& new_size_node, const Long offset = 0,
                       const Int whence = SEEK_SET)
// will read from offset relative to whence
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_read_field");
  if (not does_file_exist_qar_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  std::vector<Field<M> > fs;
  const std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i], multiplicity);
  }
  const Int mpi_tag = 7;
  if (get_id_node() == 0) {
    qassert(fs.size() > 0);
    Field<M> f;
    f.init(fs[0].geo(), multiplicity);
    Vector<M> v = get_data(f);
    const Int num_node = get_num_node();
    const Int new_num_node = product(new_size_node);
    QFile qfile = qfopen(path, "r");
    qassert(not qfile.null());
    qfseek(qfile, offset, whence);
    for (Int new_id_node = 0; new_id_node < new_num_node; ++new_id_node) {
      const Int id_node =
          get_id_node_from_new_id_node(new_id_node, new_num_node, num_node);
      qread_data(v, qfile);
      if (0 == id_node) {
        assign(get_data(fs[new_id_node]), v);
      } else {
        mpi_send((void*)v.data(), v.data_size(), MPI_BYTE, id_node, mpi_tag,
                 get_comm());
      }
    }
    qfclose(qfile);
  } else {
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      mpi_recv(v.data(), v.data_size(), MPI_BYTE, 0, mpi_tag, get_comm());
    }
  }
  shuffle_field_back(f, fs, new_size_node);
  const Long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
Long serial_read_field_par(Field<M>& f, const std::string& path,
                           const Coordinate& new_size_node,
                           const Long offset = 0, const Int whence = SEEK_SET)
// will read from offset relative to whence
// assume new_size_node is properly choosen so that concatenate the new fields
// would be correct. eg. new_size_node = Coordinate(1,1,1,2)
{
  TIMER_VERBOSE_FLOPS("serial_read_field_par");
  if (not does_file_exist_qar_sync_node(path)) {
    displayln_info(fname +
                   ssprintf(": file does not exist: '%s'", path.c_str()));
    return 0;
  }
  const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  std::vector<Field<M> > fs;
  const std::vector<Geometry> new_geos =
      make_dist_io_geos(geo.total_site(), new_size_node);
  fs.resize(new_geos.size());
  for (size_t i = 0; i < fs.size(); ++i) {
    fs[i].init(new_geos[i], multiplicity);
  }
  if (fs.size() > 0) {
    QFile qfile = qfopen(path, "r");
    qassert(not qfile.null());
    qfseek(qfile,
           offset + fs[0].geo().geon.id_node * get_data(fs[0]).data_size(),
           whence);
    for (size_t i = 0; i < fs.size(); ++i) {
      Vector<M> v = get_data(fs[i]);
      qread_data(v, qfile);
    }
    qfclose(qfile);
  }
  shuffle_field_back(f, fs, new_size_node);
  SYNC_NODE();
  const Long file_size = get_data(f).data_size() * f.geo().geon.num_node;
  timer.flops += file_size;
  return file_size;
}

template <class M>
Long serial_write_field(const Field<M>& f, const std::string& path)
// interface_function
{
  return serial_write_field(
      f, path, get_default_serial_new_size_node(f.geo(), dist_write_par_limit()));
}

template <class M>
Long serial_read_field(Field<M>& f, const std::string& path,
                       const Long offset = 0, const Int whence = SEEK_SET)
// interface_function
{
  return serial_read_field(
      f, path, get_default_serial_new_size_node(f.geo(), dist_read_par_limit()),
      offset, whence);
}

template <class M>
Long serial_read_field_par(Field<M>& f, const std::string& path,
                           const Long offset = 0, const Int whence = SEEK_SET)
// interface_function
{
  return serial_read_field_par(
      f, path, get_default_serial_new_size_node(f.geo(), dist_read_par_limit()),
      offset, whence);
}

}  // namespace qlat
